use std::collections::BTreeMap;
use std::sync::mpsc::channel;
use std::sync::{Arc, Mutex};

use error_stack::ResultExt;
use log::info;
use rust_htslib::bam::{self, Read};
use rust_vc_utils::{get_region_segments, ChromList, GenomeRef, ProgressReporter};
use unwrap::unwrap;

use crate::chrom_sites::{ChromSiteData, SiteData};
pub use crate::meth_read_processor::{
    ChromSiteMethInfo, ChromSiteMethTrainingData, GenomeSiteMethTrainingData,
};
use crate::meth_read_processor::{MethReadsProcessor, MethReadsProcessorOptions};
use crate::tflite_model::{ModelSelection, TFLiteModelData};

/// Data passed through the bam scanner methods used to improve error reports
struct BamScannerDebugInfo {
    bam_filename: String,
    chrom_list: ChromList,
}

#[allow(clippy::too_many_arguments)]
fn scan_chromosome_segment(
    thread_data: &mut WorkerThreadData,
    chrom_index: usize,
    start: i64,
    end: i64,
    chrom_size: u64,
    mrp_options: &MethReadsProcessorOptions,
    chrom_training_sites: Option<&ChromSiteData>,
    chrom_ref: Option<&Vec<u8>>,
    debug_info: &BamScannerDebugInfo,
    progress_reporter: &ProgressReporter,
) -> (ChromSiteMethInfo, ChromSiteMethTrainingData) {
    // Note that ideally we'd like to put &mut thread_data.tflite into meth_reads_processor too, but
    // so far I haven't figured out how to annotate the lifetimes...
    //
    let mut meth_reads_processor =
        MethReadsProcessor::new(start, end, mrp_options, chrom_training_sites, chrom_ref);

    let bam_reader = &mut thread_data.bam_reader;

    bam_reader
        .fetch(bam::FetchDefinition::Region(
            chrom_index as i32,
            start,
            chrom_size as i64,
        ))
        .unwrap();

    let mut record = bam::Record::new();
    while let Some(r) = bam_reader.read(&mut record) {
        unwrap!(r, "Failed to parse alignment record");

        // Remove filter to make this compatible with python script:
        //
        //if filter_out_alignment_record(&record) {
        //    continue;
        //}

        let result = meth_reads_processor
            .process_bam_record(&record, &mut thread_data.tflite)
            .attach_printable_lazy(|| {
                let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
                format!(
                    "Failed to process read `{qname}' \
                 while scanning region {}:{}-{} \
                 from alignment file `{}`",
                    debug_info.chrom_list.data[chrom_index].label,
                    start + 1,
                    end,
                    debug_info.bam_filename
                )
            });

        match result {
            Ok(is_continue) => {
                if !is_continue {
                    break;
                }
            }
            Err(err) => {
                panic!("\nError: {err:?}");
            }
        };
    }

    let result = meth_reads_processor.complete_processing(&mut thread_data.tflite);

    let kb_completed = (end - start) as u64 / 1000;
    progress_reporter.inc(kb_completed);

    result
}

/// Scan chromosome and delegate out segments of the chromosome to additional worker threads
///
#[allow(clippy::too_many_arguments)]
fn scan_chromosome_segments(
    worker_thread_data: Arc<Vec<Mutex<WorkerThreadData>>>,
    chrom_index: usize,
    chrom_size: u64,
    mrp_options: &MethReadsProcessorOptions,
    chrom_training_sites: Option<&ChromSiteData>,
    chrom_ref: Option<&Vec<u8>>,
    debug_info: &BamScannerDebugInfo,
    progress_reporter: &ProgressReporter,
) -> (ChromSiteMethInfo, ChromSiteMethTrainingData) {
    let (tx, rx) = channel();

    rayon::scope(move |scope| {
        for (start, end) in get_region_segments(chrom_size, mrp_options.segment_size) {
            let worker_thread_data = worker_thread_data.clone();
            let tx = tx.clone();

            scope.spawn(move |_| {
                let worker_id = rayon::current_thread_index().unwrap();

                let segment_data = scan_chromosome_segment(
                    &mut worker_thread_data[worker_id].lock().unwrap(),
                    chrom_index,
                    start as i64,
                    end as i64,
                    chrom_size,
                    mrp_options,
                    chrom_training_sites,
                    chrom_ref,
                    debug_info,
                    progress_reporter,
                );
                tx.send(segment_data).unwrap();
            });
        }
    });

    let mut chrom_meth_info = ChromSiteMethInfo::new();
    let mut chrom_site_meth_training_data = ChromSiteMethTrainingData::new();
    for (segment_meth_info, segment_site_meth_training_data) in rx {
        chrom_meth_info.merge(segment_meth_info);
        chrom_site_meth_training_data.merge(segment_site_meth_training_data);
    }
    (chrom_meth_info, chrom_site_meth_training_data)
}

/// Data that are persistent to each worker thread
struct WorkerThreadData<'a> {
    bam_reader: bam::IndexedReader,
    tflite: Option<TFLiteModelData<'a>>,
}

impl WorkerThreadData<'_> {
    fn new(
        bam_filename: &str,
        ref_filename: &Option<String>,
        model_selection: &ModelSelection,
    ) -> Self {
        let tflite = match model_selection {
            ModelSelection::None => None,
            _ => Some(TFLiteModelData::new(model_selection)),
        };
        let mut val = Self {
            bam_reader: bam::IndexedReader::from_path(bam_filename).unwrap(),
            tflite,
        };

        if let Some(ref_filename) = ref_filename {
            val.bam_reader.set_reference(ref_filename).unwrap();
        }

        val
    }
}

/// Site methylation data for the whole genome
///
pub struct GenomeSiteMethInfo {
    pub chrom_list: ChromList,
    pub chroms: Vec<ChromSiteMethInfo>,
}

/// Read single alignment file and translate into methylation pileup data structures
///
pub fn get_meth_pileup_from_bam(
    bam_filename: &str,
    ref_filename: &Option<String>,
    thread_count: usize,
    mrp_options: &MethReadsProcessorOptions,
    model_selection: &ModelSelection,
    training_sites: &Option<SiteData>,
    genome_ref: &Option<GenomeRef>,
) -> (
    Option<GenomeSiteMethInfo>,
    Option<GenomeSiteMethTrainingData>,
) {
    assert!(thread_count > 0);

    info!("Processing alignment file '{}'", bam_filename);

    // Setup per-worker-thread data structures:
    let mut worker_thread_data = Vec::new();
    for _ in 0..thread_count {
        worker_thread_data.push(Mutex::new(WorkerThreadData::new(
            bam_filename,
            ref_filename,
            model_selection,
        )));
    }
    let worker_thread_data = Arc::new(worker_thread_data);

    let chrom_list = {
        // Temporarily access the first bam_reader to get the chromosome list:
        let bam_reader = &worker_thread_data[0].lock().unwrap().bam_reader;
        ChromList::from_bam_header(bam_reader.header())
    };

    assert!(
        !chrom_list.data.is_empty(),
        "No chromosome names in alignment file header. \
        Input alignment file must be mapped."
    );

    let debug_info = BamScannerDebugInfo {
        bam_filename: bam_filename.to_string(),
        chrom_list: chrom_list.clone(),
    };
    let debug_info = &debug_info;

    let worker_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(thread_count)
        .build()
        .unwrap();

    let chrom_count = chrom_list.data.len();

    let genome_kb = chrom_list.data.iter().map(|x| x.length).sum::<u64>() / 1000;
    let progress_reporter = ProgressReporter::new(genome_kb, "Processed", "ref genome kb", false);

    let (tx, rx) = channel();

    let progress_reporter = &progress_reporter;
    let chrom_list_ref = &chrom_list;
    worker_pool.scope(move |scope| {
        for chrom_index in 0..chrom_count {
            let worker_thread_data = worker_thread_data.clone();
            let chrom_info = &chrom_list_ref.data[chrom_index];
            let chrom_size = chrom_info.length;

            let chrom_training_sites = match training_sites {
                Some(x) => x.chroms.get(&chrom_info.label),
                None => None,
            };

            let chrom_ref = genome_ref
                .as_ref()
                .map(|x| x.chroms.get(&chrom_info.label).unwrap());

            let tx = tx.clone();
            scope.spawn(move |_| {
                let chroms_site_meth_info = scan_chromosome_segments(
                    worker_thread_data,
                    chrom_index,
                    chrom_size,
                    mrp_options,
                    chrom_training_sites,
                    chrom_ref,
                    debug_info,
                    progress_reporter,
                );
                tx.send((chrom_index, chroms_site_meth_info)).unwrap();
            });
        }
    });

    let mut genome_site_meth_info_chroms = vec![ChromSiteMethInfo::new(); chrom_count];
    let mut genome_site_meth_training_data_chroms = BTreeMap::new();
    for (chrom_index, (chrom_site_meth_info, chrom_site_meth_training_data)) in rx {
        genome_site_meth_info_chroms[chrom_index] = chrom_site_meth_info;
        let chrom_info = &chrom_list_ref.data[chrom_index];
        let val = genome_site_meth_training_data_chroms
            .insert(chrom_info.label.clone(), chrom_site_meth_training_data);
        assert!(val.is_none(), "Chromosome label collision");
    }

    let genome_site_meth_info = if training_sites.is_some() {
        None
    } else {
        Some(GenomeSiteMethInfo {
            chrom_list,
            chroms: genome_site_meth_info_chroms,
        })
    };

    let genome_site_meth_training_data = if training_sites.is_none() {
        None
    } else {
        Some(GenomeSiteMethTrainingData {
            chroms: genome_site_meth_training_data_chroms,
        })
    };

    progress_reporter.clear();

    (genome_site_meth_info, genome_site_meth_training_data)
}
