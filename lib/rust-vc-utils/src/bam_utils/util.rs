//! Miscelanious BAM record processing utilities that don't fit any other modules
//!

use rust_htslib::{bam, htslib};

use super::cigar::{
    update_ref_and_complete_read_pos, update_ref_and_hard_clipped_read_pos, update_ref_pos,
};

/// Recreation of htslib hts_reg2bin in rust, I've copied this from noodles (MIT License)
/// with minor type adaptions.
///
/// begin and end should follow bed zero-based half-closed format
///
fn hts_reg2bin(begin: usize, end: usize, min_shift: u8, depth: u8) -> usize {
    let end = end - 1;
    let mut l = depth;
    let mut s = min_shift;
    let mut t = ((1 << (depth * 3)) - 1) / 7;

    while l > 0 {
        if begin >> s == end >> s {
            return t + (begin >> s);
        }

        l -= 1;
        s += 3;
        t -= 1 << (l * 3);
    }

    0
}

/// Recreation of htslib bam_reg2bin in rust
///
/// begin and end should follow bed zero-based half-closed format
///
pub fn bam_reg2bin(begin: usize, end: usize) -> u16 {
    hts_reg2bin(begin, end, 14, 5) as u16
}

/// Check if the alignment record should be filtered from consideration in any part of the variant
/// calling pipeline
///
pub fn filter_out_alignment_record(record: &bam::Record) -> bool {
    static FLAG_FILTER: u32 =
        htslib::BAM_FUNMAP | htslib::BAM_FSECONDARY | htslib::BAM_FQCFAIL | htslib::BAM_FDUP;

    ((record.flags() as u32) & FLAG_FILTER) != 0
}

/// Report the end reference position of a bam record
///
/// The end position is the zero-indexed right-most mapped position + 1
///
pub fn get_alignment_end(record: &bam::Record) -> i64 {
    let mut ref_pos = record.pos();
    for c in record.cigar().iter() {
        update_ref_pos(c, &mut ref_pos);
    }
    ref_pos
}

/// Get gap-compressed sequence identity from the bam record cigar
///
/// Metric is discussed here:
/// <https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity>
///
/// This method requires the use of the match (=) and mismatch (X) cigar
/// operations. The alignment match (M) operation will trigger a runtime error.
///
pub fn get_gap_compressed_identity(record: &bam::Record) -> f64 {
    use rust_htslib::bam::record::Cigar::*;

    let mut mismatch_events = 0u32;
    let mut match_bases = 0u32;
    for c in record.cigar().iter() {
        match c {
            Ins(_) | Del(_) => {
                mismatch_events += 1;
            }
            Diff(len) => {
                mismatch_events += len;
            }
            Equal(len) => {
                match_bases += len;
            }
            Match(_) => {
                panic!("Method assumes alignment CIGAR strings use seq match/mismatch (=/X) instead of alignment match (M)");
            }
            RefSkip(_) | SoftClip(_) | HardClip(_) | Pad(_) => {}
        }
    }

    match_bases as f64 / (match_bases + mismatch_events) as f64
}

/// Create an array of read length, mapping from the read to the reference position
///
/// Read positions which are not mapped to a reference position are None. Read positions include
/// hard-clipped segments.
///
/// Ref positions are 0-indexed
///
pub fn get_complete_read_to_ref_pos_map(record: &bam::Record) -> Vec<Option<i64>> {
    use rust_htslib::bam::record::Cigar::*;

    assert!(!record.is_unmapped());

    let mut read_to_ref = vec![None; record.seq_len()];
    let mut read_pos = 0usize;
    let mut ref_pos = record.pos();

    for c in record.cigar().iter() {
        match c {
            Diff(len) | Equal(len) | Match(len) => {
                let len = *len as usize;
                for i in 0..len {
                    read_to_ref[read_pos + i] = Some(ref_pos + i as i64);
                }
            }
            _ => {}
        }
        update_ref_and_complete_read_pos(c, &mut ref_pos, &mut read_pos);
    }
    read_to_ref
}

/// Create an array of read length, mapping from the read to the reference position
///
/// Read positions which are not mapped to a reference position are None. Read positions do not
/// include hard-clipped segments.
///
/// Ref positions are 0-indexed
///
pub fn get_hard_clipped_read_to_ref_pos_map(record: &bam::Record) -> Vec<Option<i64>> {
    use rust_htslib::bam::record::Cigar::*;

    assert!(!record.is_unmapped());

    let mut read_to_ref = vec![None; record.seq_len()];
    let mut read_pos = 0usize;
    let mut ref_pos = record.pos();

    for c in record.cigar().iter() {
        match c {
            Diff(len) | Equal(len) | Match(len) => {
                let len = *len as usize;
                for i in 0..len {
                    read_to_ref[read_pos + i] = Some(ref_pos + i as i64);
                }
            }
            _ => {}
        }
        update_ref_and_hard_clipped_read_pos(c, &mut ref_pos, &mut read_pos);
    }
    read_to_ref
}

/// Extract sample name from bam header
///
/// This uses the sample name from the first read group found in the header, and does not
/// check for additional read groups. If no read group is found, `default_name` will
/// be used.
///
pub fn get_sample_name(header: &bam::HeaderView, default_name: &str) -> String {
    for line in std::str::from_utf8(header.as_bytes()).unwrap().split('\n') {
        for (i, word) in line.split('\t').enumerate() {
            if i == 0 {
                if word != "@RG" {
                    break;
                }
            } else if let Some(sample_name) = word.strip_prefix("SM:") {
                return sample_name.to_string();
            }
        }
    }
    default_name.to_string()
}

/// Translate the (zero-indexed) read position into the corresponding read position in the reverse orientation
///
pub fn get_reverse_read_position(record: &bam::Record, read_pos: usize) -> usize {
    let read_len = record.seq_len();
    if read_pos >= read_len {
        let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
        panic!(
            "Invalid read position {read_pos}, exceeds the read_length {read_len}, in read {qname}"
        );
    }
    read_len - (read_pos + 1)
}

/// Translate the (zero-indexed) read position in fwd-aligned read order to the read position for the
/// read in sequencer order
///
pub fn get_seq_order_read_position(record: &bam::Record, read_pos: usize) -> usize {
    if record.is_reverse() {
        get_reverse_read_position(record, read_pos)
    } else {
        read_pos
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::{header, Header, HeaderView};

    fn get_test_header() -> HeaderView {
        let mut _header = Header::new();
        _header.push_record(
            header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr1")
                .push_tag(b"LN", 10000000),
        );
        HeaderView::from_header(&_header)
    }

    #[test]
    fn test_filter_out_alignment_record() {
        let header = get_test_header();

        // Unmapped read:
        let sam_line =
            b"qname\t4\t*\t0\t255\t*\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        assert!(filter_out_alignment_record(&rec));

        // Mapped read:
        let sam_line =
            b"qname\t0\tchr1\t10\t60\t20M\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        assert!(!filter_out_alignment_record(&rec));
    }

    #[test]
    fn test_get_alignment_end() {
        let header = get_test_header();

        // Mapped read:
        let sam_line = b"qname\t0\tchr1\t10\t60\t5S5M10D5I5M\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        assert_eq!(get_alignment_end(&rec), 29);
    }

    #[test]
    fn test_get_hard_clipped_read_to_ref_pos_map() {
        let header = get_test_header();
        let sam_line = b"qname\t0\tchr1\t10\t60\t2H2M1I1M\t*\t0\t0\tACTG\tDDDD";
        let record = bam::Record::from_sam(&header, sam_line).unwrap();

        let rval = get_hard_clipped_read_to_ref_pos_map(&record);
        assert_eq!(rval, vec![Some(9), Some(10), None, Some(11)]);
    }

    #[test]
    fn test_get_seq_order_read_position() {
        let header = get_test_header();
        let sam_line =
            b"qname\t16\tchr1\t10\t60\t20M\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        assert_eq!(get_seq_order_read_position(&rec, 1), 18);
    }
}
