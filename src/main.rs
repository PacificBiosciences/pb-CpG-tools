mod bam_scanner;
mod chrom_sites;
mod cli;
mod detect_basemod_info;
mod globals;
mod meth_read_processor;
mod site_meth_prob;
mod tflite_model;
mod write_tracks;

use std::{error, process};

use hhmmss::Hhmmss;
use log::info;
use meth_read_processor::MethReadsProcessorPileupMode;
use rust_vc_utils::genome_ref;
use tflite_model::ModelSelection;

use crate::globals::{PROGRAM_NAME, PROGRAM_VERSION};
use crate::meth_read_processor::{MethReadsProcessorOptions, ModSitesMode};

fn setup_logger(output_prefix: &str) -> Result<(), fern::InitError> {
    let log_filename = output_prefix.to_owned() + ".log";
    fern::Dispatch::new()
        .format(move |out, message, record| {
            out.finish(format_args!(
                "{}[{}][{}] {}",
                chrono::Local::now().format("[%Y-%m-%d][%H:%M:%S]"),
                PROGRAM_NAME,
                record.level(),
                message
            ))
        })
        .level(log::LevelFilter::Info)
        .chain(std::io::stderr())
        .chain(fern::log_file(log_filename)?)
        .apply()?;
    Ok(())
}

/// Boil down all the settings to select the methylation pileup model
pub fn get_model_selection_from_settings(settings: &cli::Settings) -> ModelSelection {
    if settings.pileup_mode == MethReadsProcessorPileupMode::count {
        ModelSelection::None
    } else if let Some(x) = &settings.override_model_filename {
        ModelSelection::OverrideFilename(x.clone())
    } else {
        ModelSelection::Builtin
    }
}

/// Write out model training features, plus some extra metadata
///
pub fn write_all_training_data(
    output_prefix: &str,
    genome_site_meth_training_data: &bam_scanner::GenomeSiteMethTrainingData,
) {
    let json_filename = output_prefix.to_owned() + ".features.json";
    info!("Writing model features to json file: '{}'", json_filename);

    // Serialize it to a JSON string.
    let file = std::fs::File::create(json_filename).unwrap();
    serde_json::to_writer(file, &genome_site_meth_training_data).unwrap();
}

fn run(
    settings: &cli::Settings,
    derived_settings: &cli::DerivedSettings,
) -> Result<(), Box<dyn error::Error>> {
    info!("Starting {PROGRAM_NAME} {PROGRAM_VERSION}");
    info!(
        "cmdline: {}",
        std::env::args().collect::<Vec<_>>().join(" ")
    );
    info!("Running on {} threads", derived_settings.thread_count);
    let start = std::time::Instant::now();

    cli::validate_settings_data(settings);

    // Read the fasta reference
    let genome_ref = if let ModSitesMode::reference = settings.modsites_mode {
        settings
            .ref_filename
            .as_ref()
            .map(|x| genome_ref::get_genome_ref_from_fasta(x))
    } else {
        None
    };

    // Read in the sites file
    let training_sites = settings
        .training_sites_filename
        .as_ref()
        .map(|x| chrom_sites::read_sites_file(x));

    // Scan the entire genome in parallel
    let scan_result = {
        // This defines how large of a chromosome segment should be processed by a single thread:
        let segment_size = 20_000_000;

        let mrp_options = MethReadsProcessorOptions {
            min_coverage: settings.min_coverage,
            min_mapq: settings.min_mapq,
            hap_tag: settings.hap_tag.as_bytes().try_into().unwrap(),
            pileup_mode: settings.pileup_mode,
            modsites_mode: settings.modsites_mode,
            segment_size,
        };

        let model_selection = get_model_selection_from_settings(settings);
        bam_scanner::get_meth_pileup_from_bam(
            &settings.bam_filename,
            &settings.ref_filename,
            derived_settings.thread_count,
            &mrp_options,
            &model_selection,
            &training_sites,
            &genome_ref,
        )
    };

    info!("Finished processing alignment files.");

    if let Some(genome_site_meth_training_data) = scan_result.1 {
        let safe_frac = |a, b| match b {
            0 => 0f64,
            _ => a as f64 / b as f64,
        };

        let in_site_count = training_sites.unwrap().site_count();
        let out_site_count = genome_site_meth_training_data.site_count();
        info!(
            "Features generated for {} out of {} requested sites ({:.5})",
            out_site_count,
            in_site_count,
            safe_frac(out_site_count, in_site_count)
        );

        // Write out all training_data
        write_all_training_data(&settings.output_prefix, &genome_site_meth_training_data);
    }

    if let Some(genome_site_meth_info) = scan_result.0 {
        // Write out all bed and bigwig files
        write_tracks::write_all_output_files(settings, derived_settings, &genome_site_meth_info);
    }

    info!(
        "{PROGRAM_NAME} completed. Total Runtime: {}",
        start.elapsed().hhmmssxxx()
    );
    Ok(())
}

fn main() {
    let settings = cli::parse_settings();
    setup_logger(&settings.output_prefix).unwrap();
    let derived_settings = cli::validate_and_fix_settings(&settings);

    if let Err(err) = run(&settings, &derived_settings) {
        eprintln!("{err}");
        process::exit(2);
    }
}
