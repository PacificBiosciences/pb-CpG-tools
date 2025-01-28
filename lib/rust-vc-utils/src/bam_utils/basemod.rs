use std::collections::BTreeMap;

use rust_htslib::bam::{self, record::Aux};
use unwrap::unwrap;

fn get_mm_tag(rec: &bam::Record) -> Option<Aux> {
    if let Ok(value) = rec.aux(b"MM") {
        Some(value)
    } else if let Ok(value) = rec.aux(b"Mm") {
        Some(value)
    } else {
        None
    }
}

fn get_ml_tag(rec: &bam::Record) -> Option<Aux> {
    if let Ok(value) = rec.aux(b"ML") {
        Some(value)
    } else if let Ok(value) = rec.aux(b"Ml") {
        Some(value)
    } else {
        None
    }
}

fn base_comp(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'N',
    }
}

pub fn decode_ml(ml: u8) -> f32 {
    static BIN_SIZE: f64 = 1.0 / 256.0;
    static HALF_BIN_SIZE: f64 = 1.0 / 512.0;
    (((ml as f64) * BIN_SIZE) + HALF_BIN_SIZE) as f32
}

/// Captures the interpretation of bases that are skipped in the MM tag base position sequence.
///
/// 'Unknown' corresponds to the '?' character, indicating that unspecified bases have unknown
/// methylation probability.
///
/// 'LowProb' corresponds to the '.' character, indicating that unspecified bases have low
/// methylation probability.
///
/// Although the default state is listed as 'LowProb' in the spec, it is preserved as a separate
/// state here because many older tools were not setting this default state with intent.
///
/// '''
/// Following the base modification codes is a recommended but optional ‘.’ or ‘?’ describing how skipped
/// seq bases of the stated base type should be interpreted by downstream tools. When this flag is ‘?’
/// there is no information about the modification status of the skipped bases provided. When this flag is
/// not present, or it is ‘.’, these bases should be assumed to have low probability of modification.
/// '''
///
#[derive(Default, PartialEq)]
pub enum CpGMethSkippedBaseMode {
    #[default]
    Default,
    Unknown,
    LowProb,
}

#[derive(Default)]
pub struct CpgMethInfo {
    /// Key is read coordinates of the C-base from the CpG on the fwd strand, such that these can
    /// be matched up to the cigar string. Value is methylation prob.
    ///
    pub pos_prob: BTreeMap<usize, f32>,
    pub skip_mode: CpGMethSkippedBaseMode,
}

impl CpgMethInfo {
    pub fn new(skip_mode: CpGMethSkippedBaseMode) -> Self {
        Self {
            skip_mode,
            ..Default::default()
        }
    }
}

/// Return information on the pattern of 5mC methylation at CpG motifs for the input bam record
///
/// Note that because this is restricted to CpG contexts, only positive strand methylation (C+m) is
/// reported, and negative strand is ignored if present (G-m). This is done under the assumption
/// that each CpG has a uniform methylation status on both strands.
///
/// Also note that for any read mapped to the reverse strand, the position of the methylation
/// observation will be adjusted so that it occurs at the forward-strand "C" base of the CpG.
///
/// This function should tolerate MM/ML tags that have methylation entries other than C+m, but it
/// cannot handle the 'multi-modification' format, such as 'C+hm'.
///
/// Return Err with an error code indicating one of the following issues:
/// (1) either of the MM (Mm) or ML (Ml) tags are missing from the bam record
/// (2) There are no entries in the MM tag
/// (3) there are no entries corresponding to 5mC methylation
/// (4) 5mC entry is present, but blank
/// (5) No 5mC is present in the CpG context
///
pub fn decode_cpg_meth_info(record: &bam::Record) -> Result<CpgMethInfo, i32> {
    // get MM amd ML tags from record
    let mm_tag = get_mm_tag(record);
    let ml_tag = get_ml_tag(record);

    // Some reads will have missing MM/ML tags, just skip these cases:
    if mm_tag.is_none() || ml_tag.is_none() {
        return Err(1);
    }

    // qname is only used to report errors:
    let qname = std::str::from_utf8(record.qname()).unwrap().to_string();

    // Extract underlying mm_tag string:
    let mm_tag = mm_tag.unwrap();
    let mm_tag = match mm_tag {
        Aux::String(tag) => tag,
        _ => panic!("Unexpected MM tag format in read {}: {:?}", qname, mm_tag),
    };

    // A blank entry after the MM tag is part of normal ccs output:
    if mm_tag.is_empty() {
        return Err(2);
    }

    // Get offset array for 5mC basemods:
    let mut ml_offset = 0;
    let mut offsets = None;
    let mut skip_mode = CpGMethSkippedBaseMode::default();
    for mm_segment in mm_tag.split(';') {
        let mut mm_iter = mm_segment.split(',');
        if let Some(word) = mm_iter.next() {
            if word.starts_with("C+m") {
                offsets = Some(
                    mm_iter
                        .map(|n| n.parse::<usize>().unwrap())
                        .collect::<Vec<_>>(),
                );
                if word.len() > 3 {
                    skip_mode = match word.chars().nth(3).unwrap() {
                        '?' => CpGMethSkippedBaseMode::Unknown,
                        '.' => CpGMethSkippedBaseMode::LowProb,
                        _ => panic!("Unexpected MM tag format in read {}: {:?}", qname, mm_tag),
                    };
                }
                break;
            } else {
                ml_offset += mm_iter.count()
            }
        } else {
            panic!("Unexpected MM tag format in read {}: {:?}", qname, mm_tag);
        }
    }

    let offsets = match offsets {
        None => {
            // Expected if no 5mC basemods are present in the record:
            return Err(3);
        }
        Some(x) => x,
    };

    // Expected if the 5mC modification is present in the record, but has no entries:
    if offsets.is_empty() {
        return Err(4);
    }

    // Extract corresponding segment of the ml array:
    let ml_tag = ml_tag.unwrap();
    let ml_vals = match ml_tag {
        Aux::ArrayU8(tag) => tag
            .iter()
            .skip(ml_offset)
            .take(offsets.len())
            .collect::<Vec<_>>(),
        _ => panic!("Unexpected MM tag format in read {}: {:?}", qname, mm_tag),
    };

    // Used for error messages:
    let mm_tag_offset_count = offsets.len();
    let ml_tag_value_count = ml_vals.len();
    assert_eq!(
        mm_tag_offset_count, ml_tag_value_count,
        "Error: bam record C+m MM and ML counts disagree ({} vs {}) in bam record: {}",
        mm_tag_offset_count, ml_tag_value_count, qname
    );

    // Convert offsets into read positions
    let mut read = record.seq().as_bytes();
    if record.is_reverse() {
        read = read.into_iter().rev().map(base_comp).collect::<Vec<_>>();
    }
    //eprintln!("read: {:?}", read.iter().map(|&x| x as char ).collect::<String>());
    let mut read_iter = read.iter().enumerate();
    let mut basemod_read_indexes = Vec::new();
    for &offset in offsets.iter() {
        let mut c_count = 0;
        loop {
            let (read_index, &base) = unwrap!(
                read_iter.next(),
                "Read sequence is too short for MM tag offsets in bam record: {qname}"
            );
            if base == b'C' {
                if c_count == offset {
                    basemod_read_indexes.push(read_index);
                    break;
                } else {
                    c_count += 1;
                }
            }
            assert!(c_count <= offset);
        }
    }
    assert_eq!(basemod_read_indexes.len(), offsets.len());

    // Convert into final data structure
    //
    // Filter out any non-CpG 'C''s during conversion
    let mut cpg_meth_info = CpgMethInfo::new(skip_mode);

    let read_len = read.len();
    for (&read_index, &ml_val) in basemod_read_indexes.iter().zip(ml_vals.iter()) {
        assert!(read_index < read_len);
        assert_eq!(read[read_index], b'C');

        // Test for CpG
        if read_index + 1 >= read_len || read[read_index + 1] != b'G' {
            continue;
        }

        let mut ref_strand_read_index = read_index;
        if record.is_reverse() {
            // Note this isn't a standard index reversal, there is an additional -1 so that we flip
            // from the "G" to the "C" index of the "CpG" context on the forward strand.
            ref_strand_read_index = read_len - (read_index + 1) - 1;
        }

        cpg_meth_info
            .pos_prob
            .insert(ref_strand_read_index, decode_ml(ml_val));
    }

    if cpg_meth_info.pos_prob.is_empty() {
        return Err(5);
    }

    Ok(cpg_meth_info)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::{header, Header, HeaderView};

    #[test]
    fn test_mm_ml_parse() {
        let mut _header = Header::new();
        _header.push_record(
            header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr1")
                .push_tag(b"LN", 10000000),
        );
        let header = HeaderView::from_header(&_header);

        //
        // test an unmapped read with no MM/ML tags:
        //
        let sam_line =
            b"qname\t4\t*\t0\t255\t*\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let info = decode_cpg_meth_info(&rec);
        assert_eq!(1, info.err().unwrap());

        //
        // test an unmapped read with blank MM/ML tags:
        //
        let sam_line = b"qname\t4\t*\t0\t255\t*\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE\tMM:Z:\tMl:B:C";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let info = decode_cpg_meth_info(&rec);
        assert_eq!(2, info.err().unwrap());

        //
        // test an unmapped read with empty basemod section in MM tags:
        //
        let sam_line = b"qname\t4\t*\t0\t255\t*\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE\tMM:Z:C+m;\tMl:B:C";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let info = decode_cpg_meth_info(&rec);
        assert_eq!(4, info.err().unwrap());

        //
        // test an unmapped read with non-5mC MM/ML tags:
        //
        let sam_line = b"qname\t4\t*\t0\t255\t*\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE\tMM:Z:A+m,1,0;\tMl:B:C,100,150";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let info = decode_cpg_meth_info(&rec);
        assert_eq!(3, info.err().unwrap());

        //
        // test an unmapped read with MM/ML tags:
        //
        let sam_line = b"qname\t4\t*\t0\t255\t*\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE\tMM:Z:A+m,1,0;C+m,0,1,1;\tMl:B:C,100,150,200,220,240";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let info = decode_cpg_meth_info(&rec);

        assert!(info.is_ok());
        let info = info.unwrap();
        assert!(info.skip_mode == CpGMethSkippedBaseMode::Default);
        let mut pp_iter = info.pos_prob.iter();

        let val = pp_iter.next();
        assert!(val.is_some());
        let (&read_index, &prob) = val.unwrap();
        assert_eq!(read_index, 1);
        approx::assert_ulps_eq!(prob, decode_ml(200), max_ulps = 4);

        let val = pp_iter.next();
        assert!(val.is_some());
        let (&read_index, &prob) = val.unwrap();
        assert_eq!(read_index, 4);
        approx::assert_ulps_eq!(prob, decode_ml(220), max_ulps = 4);

        let val = pp_iter.next();
        assert!(val.is_none());

        //
        // test an unmapped read with MM/ML tags using the optional skipped-base character '?':
        //
        let sam_line = b"qname\t4\t*\t0\t255\t*\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE\tMM:Z:A+m,1,0;C+m?,0,1,1;\tMl:B:C,100,150,200,220,240";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let info = decode_cpg_meth_info(&rec);

        assert!(info.is_ok());
        let info = info.unwrap();
        assert!(info.skip_mode == CpGMethSkippedBaseMode::Unknown);
        let mut pp_iter = info.pos_prob.iter();

        let val = pp_iter.next();
        assert!(val.is_some());
        let (&read_index, &prob) = val.unwrap();
        assert_eq!(read_index, 1);
        approx::assert_ulps_eq!(prob, decode_ml(200), max_ulps = 4);

        let val = pp_iter.next();
        assert!(val.is_some());
        let (&read_index, &prob) = val.unwrap();
        assert_eq!(read_index, 4);
        approx::assert_ulps_eq!(prob, decode_ml(220), max_ulps = 4);

        let val = pp_iter.next();
        assert!(val.is_none());

        //
        // test an unmapped read with MM/ML tags using the optional skipped-base character '.':
        //
        let sam_line = b"qname\t4\t*\t0\t255\t*\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE\tMM:Z:A+m,1,0;C+m.,0,1,1;\tMl:B:C,100,150,200,220,240";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let info = decode_cpg_meth_info(&rec);

        assert!(info.is_ok());
        let info = info.unwrap();
        assert!(info.skip_mode == CpGMethSkippedBaseMode::LowProb);
        let mut pp_iter = info.pos_prob.iter();

        let val = pp_iter.next();
        assert!(val.is_some());
        let (&read_index, &prob) = val.unwrap();
        assert_eq!(read_index, 1);
        approx::assert_ulps_eq!(prob, decode_ml(200), max_ulps = 4);

        let val = pp_iter.next();
        assert!(val.is_some());
        let (&read_index, &prob) = val.unwrap();
        assert_eq!(read_index, 4);
        approx::assert_ulps_eq!(prob, decode_ml(220), max_ulps = 4);

        let val = pp_iter.next();
        assert!(val.is_none());

        //
        // test a fwd-strand mapped read with MM/ML tags:
        //
        let sam_line = b"qname\t0\tchr1\t10\t60\t20M\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE\tMM:Z:A+m,1,0;C+m,0,1,1;\tMl:B:C,100,150,200,220,240";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let info = decode_cpg_meth_info(&rec);

        assert!(info.is_ok());
        let info = info.unwrap();
        let mut pp_iter = info.pos_prob.iter();

        let val = pp_iter.next();
        assert!(val.is_some());
        let (&read_index, &prob) = val.unwrap();
        assert_eq!(read_index, 1);
        approx::assert_ulps_eq!(prob, decode_ml(200), max_ulps = 4);

        let val = pp_iter.next();
        assert!(val.is_some());
        let (&read_index, &prob) = val.unwrap();
        assert_eq!(read_index, 4);
        approx::assert_ulps_eq!(prob, decode_ml(220), max_ulps = 4);

        let val = pp_iter.next();
        assert!(val.is_none());

        //
        // test a rev-strand mapped read with MM/ML tags:
        //
        let sam_line = b"qname\t16\tchr1\t10\t60\t20M\t*\t0\t0\tTCCTCGAGACGATACGGCGT\tEEEEEDDDDDEEEEEDDDDD\tMM:Z:A+m,1,0;C+m,0,1,1;\tMl:B:C,100,150,200,220,240";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let info = decode_cpg_meth_info(&rec);

        assert!(info.is_ok());
        let info = info.unwrap();
        let mut pp_iter = info.pos_prob.iter();

        let val = pp_iter.next();
        assert!(val.is_some());
        let (&read_index, &prob) = val.unwrap();
        assert_eq!(read_index, 14);
        approx::assert_ulps_eq!(prob, decode_ml(220), max_ulps = 4);

        let val = pp_iter.next();
        assert!(val.is_some());
        let (&read_index, &prob) = val.unwrap();
        assert_eq!(read_index, 17);
        approx::assert_ulps_eq!(prob, decode_ml(200), max_ulps = 4);

        let val = pp_iter.next();
        assert!(val.is_none());

        //
        // test an unmapped read with MM/ML tags and final C
        //
        let sam_line =
            b"qname\t4\t*\t0\t255\t*\t*\t0\t0\tACGC\tDDDD\tMM:Z:C+m,0,0;\tMl:B:C,100,100";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let info = decode_cpg_meth_info(&rec);

        assert!(info.is_ok());
        let info = info.unwrap();
        assert!(info.skip_mode == CpGMethSkippedBaseMode::Default);
        let mut pp_iter = info.pos_prob.iter();

        let val = pp_iter.next();
        assert!(val.is_some());
        let (&read_index, &prob) = val.unwrap();
        assert_eq!(read_index, 1);
        approx::assert_ulps_eq!(prob, decode_ml(100), max_ulps = 4);

        let val = pp_iter.next();
        assert!(val.is_none());
    }
}
