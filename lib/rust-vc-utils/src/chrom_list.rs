use std::collections::HashMap;

use rust_htslib::bam;

#[derive(Clone, Debug)]
pub struct ChromInfo {
    pub label: String,
    pub length: u64,
}

impl PartialEq for ChromInfo {
    fn eq(&self, other: &Self) -> bool {
        (self.label == other.label) && (self.length == other.length)
    }
}

/// Facilitates translation between an ordered list of chromosomes and codes
///
#[derive(Clone, Debug, Default)]
pub struct ChromList {
    pub data: Vec<ChromInfo>,
    pub label_to_index: HashMap<String, usize>,
}

impl ChromList {
    pub fn from_bam_header(header: &bam::HeaderView) -> Self {
        let chrom_count = header.target_count();
        let mut chrom_list = ChromList::default();
        for tid in 0..chrom_count {
            chrom_list.add_chrom(
                std::str::from_utf8(header.tid2name(tid)).unwrap(),
                header.target_len(tid).unwrap(),
            );
        }
        chrom_list
    }

    pub fn from_bam_filename(bam_filename: &str) -> Self {
        use rust_htslib::bam::Read;

        let bam_reader = bam::Reader::from_path(bam_filename).unwrap();
        Self::from_bam_header(bam_reader.header())
    }

    pub fn add_chrom(&mut self, label: &str, length: u64) {
        assert!(!self.label_to_index.contains_key(label));

        self.label_to_index
            .insert(label.to_string(), self.data.len());
        self.data.push(ChromInfo {
            label: label.to_string(),
            length,
        });
    }
}

impl PartialEq for ChromList {
    fn eq(&self, other: &Self) -> bool {
        let data_size = self.data.len();
        if data_size != other.data.len() {
            return false;
        }

        for i in 0..data_size {
            if self.data[i] != other.data[i] {
                return false;
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chrom_list_from_bam() {
        //eprintln!("cwd for test: {:?}", std::env::current_dir());
        let test_bam_file = "./test_data/test_empty.bam";

        let chrom_list = ChromList::from_bam_filename(test_bam_file);

        assert_eq!(chrom_list.data.len(), 4);
        assert_eq!(
            chrom_list.data[3],
            ChromInfo {
                label: "chr4".to_string(),
                length: 14,
            }
        );
        assert_eq!(*chrom_list.label_to_index.get("chr4").unwrap(), 3);
    }
}
