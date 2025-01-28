use std::collections::HashMap;
use std::fs::File;

use bio::io::fasta;
use log::info;
use unwrap::unwrap;

#[derive(Default)]
pub struct GenomeRef {
    /// A map from chrom name to chrom sequence
    pub chroms: HashMap<String, Vec<u8>>,
}

pub fn get_genome_ref_from_fasta_fp(file: File) -> GenomeRef {
    let reader = fasta::Reader::new(file);

    let mut genome_ref = GenomeRef::default();

    for result in reader.records() {
        let record = result.expect("Error during fasta record parsing");

        genome_ref
            .chroms
            .insert(record.id().to_string(), record.seq().to_ascii_uppercase());

        /*
        let chrom = record.id().to_string();
        let seq = genome_ref.chroms.get(chrom.as_str()).unwrap();
        let mini = seq.iter().take(10).copied().collect::<Vec<_>>();
        println!("Read chrom {} length {}", chrom, seq.len());
        println!("Read mini {}", std::str::from_utf8(&mini).unwrap());
         */
    }
    genome_ref
}

/// Read fasta file into GenomeRef data structure
///
pub fn get_genome_ref_from_fasta(filename: &str) -> GenomeRef {
    info!("Reading reference genome from file '{filename}'");

    let file = unwrap!(
        File::open(filename),
        "Unable to open reference fasta file: '{}'",
        filename,
    );

    get_genome_ref_from_fasta_fp(file)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{Seek, SeekFrom, Write};

    #[test]
    fn test_get_genome_ref_from_fasta_fp() {
        let mut file = tempfile::tempfile().unwrap();

        let cname = "foo";
        let seq = "ACGTACGT";
        writeln!(file, ">{cname}").unwrap();
        writeln!(file, "{seq}").unwrap();
        file.seek(SeekFrom::Start(0)).unwrap();
        let result = get_genome_ref_from_fasta_fp(file);

        assert_eq!(result.chroms.len(), 1);
        assert_eq!(result.chroms.keys().next().unwrap(), cname);
        assert_eq!(
            std::str::from_utf8(result.chroms.values().next().unwrap()).unwrap(),
            seq
        );
    }
}
