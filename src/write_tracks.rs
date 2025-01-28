use log::info;
use rust_vc_utils::bigwig_utils;
use unwrap::unwrap;

use crate::bam_scanner::GenomeSiteMethInfo;
use crate::cli::{DerivedSettings, Settings};
use crate::meth_read_processor::{
    self, get_python_script_bed_labels, MethInfoTrackType, MethReadsProcessorPileupMode,
    SiteMethInfo,
};

/// Write a single bigwig track of site methylation probabilities
fn write_site_meth_prob_bigwig_track(
    output_prefix: &str,
    genome_site_meth_info: &GenomeSiteMethInfo,
    tt: MethInfoTrackType,
) {
    let label = meth_read_processor::get_meth_info_track_label(tt);
    let bw_filename = output_prefix.to_owned() + "." + label + ".bw";
    info!(
        "Writing {} site methylation to bigwig file: '{}'",
        label, bw_filename
    );

    let chrom_list = &genome_site_meth_info.chrom_list;
    let mut bigwig_writer = bigwig_utils::get_new_writer(&bw_filename, chrom_list);

    for (chrom_index, chrom_site_meth_info) in genome_site_meth_info.chroms.iter().enumerate() {
        let chrom_name = chrom_list.data[chrom_index].label.as_str();
        let mut start_pos: Vec<u32> = Vec::new();
        let mut site_meth_scores: Vec<f32> = Vec::new();
        for (&ref_pos, site_meth_info) in chrom_site_meth_info.get_track(tt).iter() {
            start_pos.push(ref_pos as u32);
            site_meth_scores.push(site_meth_info.pileup_output_score() as f32);
        }

        bigwig_writer
            .add_interval_spans(chrom_name, &mut start_pos, 1, &mut site_meth_scores)
            .unwrap();
    }
}

/// Logic translated from the original aligned_bam_to_cpg_scores.py script
pub fn discretize_model_score(mut score: f64, cov: u32) -> (u32, u32, f64) {
    /*
    """
        Apply a small correction to the model probability to make it
        compatible with the number of reads at that site. Allows the number
        of modified and unmodified reads to be estimated.

        :param score: Modification probability, from model. (float)
        :param coverage: Number of reads. (int)
        :return mod_reads: Estimated number of modified reads. (int)
        :return unmod_reads: Estimated number of unmodified reads. (int)
        :return adjusted_score: Adjusted probability score, based on percent modified reads. (float)
        """
         */
    //need to round up or round down modified read numbers based on score
    //which allows a push towards 0/50/100 for adjusted score

    let is_floor = if score > 50.0 {
        score < 65.0
    } else {
        score <= 35.0
    };
    score = (score / 100.0) * cov as f64;
    let mod_reads = if is_floor {
        score.floor() as u32
    } else {
        score.ceil() as u32
    };
    let unmod_reads = cov - mod_reads;

    let adj_score = if mod_reads == 0 {
        0.0
    } else {
        (1000.0 * mod_reads as f64 / cov as f64).round() / 10.0
    };

    (mod_reads, unmod_reads, adj_score)
}

fn write_site_meth_prob_bed_track_file(
    bed_filename: &str,
    genome_site_meth_info: &GenomeSiteMethInfo,
    tt: MethInfoTrackType,
    bed_header: &[String],
) {
    use rust_htslib::bgzf;
    use std::io::{BufWriter, Write};

    let f = unwrap!(
        bgzf::Writer::from_path(bed_filename),
        "Unable to create bed file: '{bed_filename}'",
    );
    let mut f = BufWriter::new(f);

    for bed_header_line in bed_header {
        writeln!(f, "{}", bed_header_line).unwrap();
    }

    let col_label = get_python_script_bed_labels(tt);

    for (chrom_index, chrom_site_meth_info) in genome_site_meth_info.chroms.iter().enumerate() {
        let chrom_name = genome_site_meth_info.chrom_list.data[chrom_index]
            .label
            .as_str();
        for (&ref_pos, site_meth_info) in chrom_site_meth_info.get_track(tt).iter() {
            let score = site_meth_info.pileup_output_score();
            write!(
                f,
                "{}\t{}\t{}\t{:.1}\t{}\t{}\t",
                chrom_name,
                ref_pos,
                ref_pos + 1,
                score,
                col_label,
                site_meth_info.cov()
            )
            .unwrap();

            match site_meth_info {
                SiteMethInfo::Model { cov, .. } => {
                    let (mod_reads, unmod_reads, adj_score) = discretize_model_score(score, *cov);
                    write!(f, "{}\t{}\t{:.1}", mod_reads, unmod_reads, adj_score).unwrap();
                }
                SiteMethInfo::Count {
                    meth_cov,
                    unmeth_cov,
                    mean_meth_prob,
                    mean_unmeth_prob,
                    ..
                } => {
                    write!(
                        f,
                        "{}\t{}\t{:.3}\t{:.3}",
                        meth_cov, unmeth_cov, mean_meth_prob, mean_unmeth_prob,
                    )
                    .unwrap();
                }
            };
            writeln!(f).unwrap();
        }
    }
}

fn make_bed_file_tabix_index(fname: &str) {
    use rust_htslib::htslib;

    let cfname = std::ffi::CString::new(fname.as_bytes()).unwrap();
    let min_shift = 0;
    let threads = 1;

    let ret = unsafe {
        htslib::tbx_index_build3(
            cfname.as_ptr(),
            std::ptr::null_mut(),
            min_shift,
            threads,
            &htslib::tbx_conf_bed,
        )
    };

    match ret {
        0 => (),
        -2 => panic!("[tabix] the compression of {fname} is not BGZF"),
        _ => panic!("tbx_index_build3 failed: {fname}"),
    }
}

/// Write a single bed track of site methylation probabilities
///
/// bed file will be bgzip compressed and tabix indexed
///
pub fn write_site_meth_prob_bed_track(
    output_prefix: &str,
    genome_site_meth_info: &GenomeSiteMethInfo,
    tt: MethInfoTrackType,
    bed_header: &[String],
) {
    let label = meth_read_processor::get_meth_info_track_label(tt);
    let bed_filename = output_prefix.to_owned() + "." + label + ".bed.gz";
    info!(
        "Writing {} site methylation to bed file: '{}'",
        label, bed_filename
    );

    write_site_meth_prob_bed_track_file(&bed_filename, genome_site_meth_info, tt, bed_header);

    make_bed_file_tabix_index(&bed_filename);
}

fn write_site_meth_prob_tracks(
    output_prefix: &str,
    genome_site_meth_info: &GenomeSiteMethInfo,
    tt: MethInfoTrackType,
    bed_header: &[String],
) {
    if meth_read_processor::skip_if_empty(tt) {
        let mut is_empty = true;
        for chrom_site_meth_info in genome_site_meth_info.chroms.iter() {
            if !chrom_site_meth_info.get_track(tt).is_empty() {
                is_empty = false;
                break;
            }
        }
        if is_empty {
            return;
        }
    }

    rayon::scope(|scope| {
        scope.spawn(move |_| {
            write_site_meth_prob_bigwig_track(output_prefix, genome_site_meth_info, tt);
        });
        scope.spawn(move |_| {
            write_site_meth_prob_bed_track(output_prefix, genome_site_meth_info, tt, bed_header);
        });
    });
}

/// Create the full header shared by all bed files
fn get_shared_bed_header(settings: &Settings, derived_settings: &DerivedSettings) -> Vec<String> {
    use crate::globals::PROGRAM_VERSION;

    let mut header = Vec::new();
    header.push(format!("##pb-cpg-tools-version={PROGRAM_VERSION}"));

    let cmdline = std::env::args().collect::<Vec<_>>().join(" ");
    header.push(format!("##cmdline={cmdline}"));

    header.push(format!("##pileup-mode={}", settings.pileup_mode));
    header.push(format!("##modsites-mode={}", settings.modsites_mode));
    header.push(format!("##min-coverage={}", settings.min_coverage));
    header.push(format!("##min-mapq={}", settings.min_mapq));

    for basemod_info in derived_settings.detected_basemod_info.iter() {
        header.push(format!(
            "##basemod-source={}\t{}",
            basemod_info.program, basemod_info.version
        ));
    }

    let extra_columns = if settings.pileup_mode == MethReadsProcessorPileupMode::model {
        "est_mod_count\test_unmod_count\tdiscretized_mod_score"
    } else {
        "mod_count\tunmod_count\tavg_mod_score\tavg_unmod_score"
    };
    header.push(format!(
        "#chrom\tbegin\tend\tmod_score\ttype\tcov\t{extra_columns}"
    ));

    header
}

/// Write all bed and bigwig outputs
///
/// The full settings object is passed through these routines to fill in run metadata in the bed file output
///
pub fn write_all_output_files(
    settings: &Settings,
    derived_settings: &DerivedSettings,
    genome_site_meth_info: &GenomeSiteMethInfo,
) {
    // Set io thread count independent of core count, because each of these threads
    // shouldn't take much cpu:
    let io_thread_count = 6;
    let worker_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(io_thread_count)
        .build()
        .unwrap();
    let bed_header = &get_shared_bed_header(settings, derived_settings);

    worker_pool.scope(|scope| {
        use MethInfoTrackType::*;
        for tt in [Combined, Hap1, Hap2] {
            scope.spawn(move |_| {
                write_site_meth_prob_tracks(
                    &settings.output_prefix,
                    genome_site_meth_info,
                    tt,
                    bed_header,
                );
            });
        }
    });
}
