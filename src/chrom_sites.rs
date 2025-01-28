use std::collections::BTreeMap;

use log::info;

pub struct ChromSiteData {
    /// Map represents chromosome position (0-indexed) and expected methylation fraction
    pub sites: BTreeMap<i64, f64>,
}

impl ChromSiteData {
    pub fn new() -> Self {
        Self {
            sites: BTreeMap::new(),
        }
    }
}

/// Specifies sites for training feature generation:
pub struct SiteData {
    pub chroms: BTreeMap<String, ChromSiteData>,
}

impl SiteData {
    pub fn new() -> Self {
        Self {
            chroms: BTreeMap::new(),
        }
    }

    pub fn site_count(&self) -> usize {
        self.chroms.values().map(|x| x.sites.len()).sum()
    }
}

pub fn read_sites_file(sites_filename: &str) -> SiteData {
    info!("Processing sites file '{}'", sites_filename);

    let mut sites = SiteData::new();

    // Build the CSV reader and iterate over each record.
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .comment(Some(b'#'))
        .from_path(sites_filename)
        .unwrap();

    for result in reader.records() {
        let record = result.unwrap();

        let chrom = record[0].to_owned();

        // Input file uses 1-index positions. switch to 0-indexed for internal methods
        let pos = record[1].parse::<i64>().ok().unwrap() - 1;
        let frac = record[2].parse::<f64>().ok().unwrap();

        let chrom_sites = sites.chroms.entry(chrom).or_insert_with(ChromSiteData::new);
        chrom_sites.sites.insert(pos, frac);
    }
    sites
}
