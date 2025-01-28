use chrono::Datelike;
use clap::Parser;
use log::{error, warn};
use simple_error::{bail, SimpleResult};

use crate::detect_basemod_info::{detect_bam_basemod_info, BasemodProgramInfo};
use crate::globals::PROGRAM_VERSION;
use crate::meth_read_processor::{MethReadsProcessorPileupMode, ModSitesMode};

#[derive(Parser)]
#[command(
    author,
    version = PROGRAM_VERSION,
    about,
    after_help = format!("Copyright (C) 2004-{}     Pacific Biosciences of California, Inc.
This program comes with ABSOLUTELY NO WARRANTY; it is intended for
Research Use Only and not for use in diagnostic procedures.", chrono::Utc::now().year()),
    help_template = "\
{before-help}{name} {version}
{author-with-newline}{about-with-newline}
{usage-heading} {usage}

{all-args}{after-help}"
)]
#[clap(rename_all = "kebab_case")]
pub struct Settings {
    /// Alignment file for input sample in BAM or CRAM format. Alignment file must be mapped and
    /// indexed, with MM/ML tags specifying 5mC methylation. If a CRAM file is provided, then
    /// `--ref` must also be specified.
    #[arg(long = "bam", value_name = "FILE")]
    pub bam_filename: String,

    /// Method to estimate site methylation from the read pileup
    #[arg(long, value_enum, value_name = "MODE", default_value_t = Default::default())]
    pub pileup_mode: MethReadsProcessorPileupMode,

    /// Method to pick 5mC scoring sites
    #[arg(long, value_enum, value_name = "MODE", default_value_t = Default::default())]
    pub modsites_mode: ModSitesMode,

    /// Genome reference in FASTA format. This is required if either (1) 'reference' modsite
    /// mode is selected or (2) input alignments are in CRAM format
    #[arg(long = "ref", value_name = "FILE")]
    pub ref_filename: Option<String>,

    /// Prefix used for all file output. If the prefix includes a directory, the directory must
    /// already exist.
    #[arg(long, value_name = "PREFIX", default_value = env!("CARGO_PKG_NAME"))]
    pub output_prefix: String,

    /// Number of threads to use. Defaults to all logical cpus detected.
    #[arg(long = "threads", value_name = "THREAD_COUNT")]
    thread_count_option: Option<usize>,

    /// Minimum site coverage. The tensorflow prediction models have their own hard-coded minimum
    /// coverage of 4, so any value below this level should have no effect in 'model' pileup mode.
    #[arg(long, default_value_t = 4)]
    pub min_coverage: u32,

    /// Minimum read mapping quality
    #[arg(long, default_value_t = 1)]
    pub min_mapq: u32,

    /// The 2-letter SAM aux tag used for per-read haplotype ids in the input bam
    #[arg(long, default_value = "HP")]
    pub hap_tag: String,

    /// Developer option to specify Tensorflow-lite model file (*.tflite) used to generate site methylation probability
    ///
    /// This model does not need to be specified for any standard use-case, the recommended model is already built into
    /// pb-CpG-tools directly.
    ///
    #[arg(hide = true, long = "override-model-file", value_name = "FILE")]
    pub override_model_filename: Option<String>,

    /// A tsv file containing chrom, position, expected meth frequency for all training sites
    ///
    /// Either this option or the model option must be specified, depending on whether
    /// the task is model inference or outputting features for model training.
    #[arg(hide = true, long = "training-sites")]
    pub training_sites_filename: Option<String>,

    /// Enables additional output for debugging
    #[arg(hide = true, long)]
    pub debug: bool,
}

/// Values immediately computed from the user settings, but not part of direct user inputs
///
pub struct DerivedSettings {
    /// Global thread count for pb-CpG-tools to use
    pub thread_count: usize,

    pub detected_basemod_info: Vec<BasemodProgramInfo>,
}

/// Validate settings and use these to produce derived settings
///
fn validate_and_fix_settings_impl(settings: &Settings) -> SimpleResult<DerivedSettings> {
    use std::ffi::OsStr;
    use std::path::Path;

    fn check_required_filename(filename: &str, label: &str) -> SimpleResult<()> {
        if filename.is_empty() {
            bail!("Must specify {label} file");
        }
        if !Path::new(&filename).exists() {
            bail!("Can't find specified {label} file: '{filename}'");
        }
        Ok(())
    }

    fn check_optional_filename(filename_opt: &Option<String>, label: &str) -> SimpleResult<()> {
        if let Some(filename) = filename_opt {
            if !Path::new(&filename).exists() {
                bail!("Can't find specified {label} file: '{filename}'");
            }
        }
        Ok(())
    }

    check_required_filename(&settings.bam_filename, "alignment")?;

    check_optional_filename(&settings.ref_filename, "reference")?;

    // Check if alignment file is cram
    let alignment_file_ext = Path::new(&settings.bam_filename)
        .extension()
        .and_then(OsStr::to_str);
    let cram_input = alignment_file_ext == Some("cram");

    if settings.ref_filename.is_none() {
        if cram_input {
            bail!(
                "The `--ref` option must be provided when input alignment file is in CRAM format"
            );
        }

        if settings.modsites_mode == ModSitesMode::reference {
            bail!("The `--ref` option must be provided when modsite mode is 'reference'");
        }
    }

    if let Some(training_sites_filename) = &settings.training_sites_filename {
        check_required_filename(training_sites_filename, "training sites")?;

        match settings.pileup_mode {
            MethReadsProcessorPileupMode::model => {}
            _ => {
                bail!("Training feature output can only be selected when pileup mode is 'model'");
            }
        };
    }

    if let Some(override_model_filename) = &settings.override_model_filename {
        check_required_filename(override_model_filename, "tflite model")?;
        if settings.pileup_mode == MethReadsProcessorPileupMode::count {
            warn!("A tflite model file has been specified, but will not be used in the 'count' pileup mode");
        }
    }

    if settings.pileup_mode == MethReadsProcessorPileupMode::model && cfg!(not(feature = "tflite"))
    {
        bail!("The model pileup mode has been specified without enabling the corresponding tflite build feature.");
    }

    if settings.hap_tag.len() != 2 {
        bail!("SAM haplotype tag must have a length of 2");
    }

    let thread_count = match settings.thread_count_option {
        Some(count) => {
            if count == 0 {
                bail!("--threads argument must be greater than 0");
            }
            count
        }
        None => num_cpus::get(),
    };

    let detected_basemod_info = detect_bam_basemod_info(&settings.bam_filename);

    Ok(DerivedSettings {
        thread_count,
        detected_basemod_info,
    })
}

pub fn validate_and_fix_settings(settings: &Settings) -> DerivedSettings {
    match validate_and_fix_settings_impl(settings) {
        Ok(x) => x,
        Err(msg) => {
            error!("Invalid command-line setting: {}", msg);
            std::process::exit(exitcode::USAGE);
        }
    }
}

/// Extended input data/settings validation that's too complex/slow to put in the cmdline parser
///
pub fn validate_settings_data(settings: &Settings) {
    use rust_htslib::bam::{self, Read};
    use rust_vc_utils::ChromList;

    // Test that the BAM file has an index recognized by htslib.
    //
    if let Err(error) = bam::IndexedReader::from_path(&settings.bam_filename) {
        error!("Failed to open input bam file: {}", error);
        std::process::exit(exitcode::USAGE);
    }

    // Pull chromosome list from bam
    let chrom_list = {
        let bam_reader = bam::Reader::from_path(&settings.bam_filename).unwrap();
        ChromList::from_bam_header(bam_reader.header())
    };
    if chrom_list.data.is_empty() {
        error!(
            "No chromosome names found in header of alignment file `{}`.\n\
        Input alignment file must be mapped.",
            settings.bam_filename
        );
        std::process::exit(exitcode::DATAERR);
    }
}

pub fn parse_settings() -> Settings {
    Settings::parse()
}
