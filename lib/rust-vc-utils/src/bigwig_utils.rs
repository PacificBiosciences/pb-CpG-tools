use rust_libbigwig::BigWigWriter;

use crate::ChromList;

/// Convenience function to create new bigwig writer from the ChromList object
///
pub fn get_new_writer(filename: &str, chrom_list: &ChromList) -> BigWigWriter {
    BigWigWriter::new(
        filename,
        &chrom_list
            .data
            .iter()
            .map(|s| s.label.clone())
            .collect::<Vec<_>>(),
        &chrom_list
            .data
            .iter()
            .map(|s| s.length as u32)
            .collect::<Vec<_>>(),
    )
    .unwrap()
}
