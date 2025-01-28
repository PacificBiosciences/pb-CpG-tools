//! BAM record cigar-processing utilities
//!

use rust_htslib::bam::record::{self, Cigar};

pub fn get_cigarseg_complete_read_offset(c: &Cigar) -> usize {
    use Cigar::*;
    match c {
        HardClip(len) | Ins(len) | SoftClip(len) | Diff(len) | Equal(len) | Match(len) => {
            *len as usize
        }
        _ => 0,
    }
}

pub fn get_cigarseg_hard_clipped_read_offset(c: &Cigar) -> usize {
    use Cigar::*;
    match c {
        Ins(len) | SoftClip(len) | Diff(len) | Equal(len) | Match(len) => *len as usize,
        _ => 0,
    }
}

pub fn get_cigarseg_ref_offset(c: &Cigar) -> i64 {
    use Cigar::*;
    match c {
        Del(len) | RefSkip(len) | Diff(len) | Equal(len) | Match(len) => *len as i64,
        _ => 0,
    }
}

/// A utility method to track read positions while iterating through a cigar string
///
/// Read position here is the position in the original read, before hard-clipping
pub fn update_complete_read_pos(c: &Cigar, read_pos: &mut usize) {
    *read_pos += get_cigarseg_complete_read_offset(c);
}

/// A utility method to track read positions while iterating through a cigar string
///
/// Read position here is the position in the truncated read after hard-clipping
pub fn update_hard_clipped_read_pos(c: &Cigar, read_pos: &mut usize) {
    *read_pos += get_cigarseg_hard_clipped_read_offset(c);
}

/// A utility method to track ref positions while iterating through a cigar string
pub fn update_ref_pos(c: &Cigar, ref_pos: &mut i64) {
    *ref_pos += get_cigarseg_ref_offset(c);
}

/// A utility method to track ref and read positions while iterating through a cigar string
///
/// Read position here is the position in the original read, before hard-clipping
///
/// # Example
/// ```ignore
/// let mut ref_pos = 100;
/// let mut read_pos = 100;
/// for (index, c) in record.cigar().iter().enumerate() {
///     update_ref_and_complete_read_pos(c, &mut ref_pos, &mut read_pos);
/// }
/// ```
pub fn update_ref_and_complete_read_pos(c: &Cigar, ref_pos: &mut i64, read_pos: &mut usize) {
    update_complete_read_pos(c, read_pos);
    update_ref_pos(c, ref_pos);
}

/// A utility method to track ref and read positions while iterating through a cigar string
///
/// Read position here is the position in the original read, before hard-clipping
///
/// # Example
/// ```ignore
/// let mut ref_pos = 100;
/// let mut read_pos = 100;
/// for (index, c) in record.cigar().iter().enumerate() {
///     update_ref_and_complete_read_pos(c, &mut ref_pos, &mut read_pos);
/// }
/// ```
pub fn update_ref_and_hard_clipped_read_pos(c: &Cigar, ref_pos: &mut i64, read_pos: &mut usize) {
    update_hard_clipped_read_pos(c, read_pos);
    update_ref_pos(c, ref_pos);
}

/// Report the following positions in read coordinates, all cases ignoring any hard-clipping:
/// 1. The first position after all left-side soft clipping
/// 2. The first position of all right-side soft clipping
/// 3. The read length
///
pub fn get_hard_clipped_read_clip_positions(cigar: &[Cigar]) -> (usize, usize, usize) {
    use Cigar::*;

    let mut ref_pos = 0;
    let mut read_pos = 0;
    let mut left_clip_size = 0;
    let mut right_clip_size = 0;
    let mut left_clip = true;
    for c in cigar.iter() {
        match c {
            SoftClip(len) => {
                if left_clip {
                    left_clip_size += *len as usize;
                } else {
                    right_clip_size += *len as usize;
                }
            }
            HardClip(_) => {}
            _ => {
                left_clip = false;
            }
        };
        update_ref_and_hard_clipped_read_pos(c, &mut ref_pos, &mut read_pos);
    }
    (left_clip_size, read_pos - right_clip_size, read_pos)
}

/// Return a 3-tuple of values related to read clipping
///
/// The first two values are in zero-index read coordinates:
/// 1. The first position after all left-side hard and soft clipping
/// 2. The first position after the last non-clipped base
///
/// ...and the 3rd value:
/// 3. The read length
///
pub fn get_complete_read_clip_positions(cigar: &[Cigar]) -> (usize, usize, usize) {
    use Cigar::*;

    let mut ref_pos = 0;
    let mut read_pos = 0;
    let mut left_clip_size = 0;
    let mut right_clip_size = 0;
    let mut left_clip = true;
    for c in cigar.iter() {
        match c {
            HardClip(len) | SoftClip(len) => {
                if left_clip {
                    left_clip_size += *len as usize;
                } else {
                    right_clip_size += *len as usize;
                }
            }
            _ => {
                left_clip = false;
            }
        };
        update_ref_and_complete_read_pos(c, &mut ref_pos, &mut read_pos);
    }
    (left_clip_size, read_pos - right_clip_size, read_pos)
}

/// Report the following positions in read coordinates:
/// 1. The first position after all left-side hard clipping
/// 2. The first position of all right-side hard clipping
/// 3. The read length
///
pub fn get_complete_read_hard_clip_positions(cigar: &[Cigar]) -> (usize, usize, usize) {
    use Cigar::*;

    let mut ref_pos = 0;
    let mut read_pos = 0;
    let mut left_clip_size = 0;
    let mut right_clip_size = 0;
    let mut left_clip = true;
    for c in cigar.iter() {
        match c {
            HardClip(len) => {
                if left_clip {
                    left_clip_size += *len as usize;
                } else {
                    right_clip_size += *len as usize;
                }
            }
            _ => {
                left_clip = false;
            }
        };
        update_ref_and_complete_read_pos(c, &mut ref_pos, &mut read_pos);
    }
    (left_clip_size, read_pos - right_clip_size, read_pos)
}

/// Report the reference and read offset of the cigar alignment
///
/// Read position here is the position in the original read, before hard-clipping
pub fn get_cigar_ref_and_complete_read_offset(cigar: &[Cigar]) -> (i64, usize) {
    let mut read_pos = 0;
    let mut ref_pos = 0;
    for c in cigar.iter() {
        update_ref_and_complete_read_pos(c, &mut ref_pos, &mut read_pos);
    }
    (ref_pos, read_pos)
}

/// Report the complete read offset of the cigar alignment
///
pub fn get_cigar_complete_read_offset(cigar: &[Cigar]) -> usize {
    let mut read_pos = 0;
    for c in cigar.iter() {
        update_complete_read_pos(c, &mut read_pos);
    }
    read_pos
}

/// Report the read offset of the cigar alignment after hard-clipping is removed
///
pub fn get_cigar_hard_clipped_read_offset(cigar: &[Cigar]) -> usize {
    let mut read_pos = 0;
    for c in cigar.iter() {
        update_hard_clipped_read_pos(c, &mut read_pos);
    }
    read_pos
}

/// Report the reference offset of the cigar alignment
///
pub fn get_cigar_ref_offset(cigar: &[Cigar]) -> i64 {
    let mut ref_pos = 0;
    for c in cigar.iter() {
        update_ref_pos(c, &mut ref_pos);
    }
    ref_pos
}

/// Return true if any part of the alignment is hard-clipped
///
pub fn is_hard_clipped(cigar: &[Cigar]) -> bool {
    cigar
        .iter()
        .any(|x| matches!(x, record::Cigar::HardClip { .. }))
}

/// Convert CIGAR in string format into the format used in the this library
///
/// This is a convenience function so that client code doesn't need to add their own (potentially
/// conflicting) rust-htslib dependency
///
pub fn get_cigar_from_string(cigar_str: &str) -> Vec<Cigar> {
    record::CigarString::try_from(cigar_str.as_bytes())
        .unwrap()
        .into()
}

/// Convert any matching adjacent cigar elements into a single element
///
pub fn compress_cigar(cigar_in: &[Cigar]) -> Vec<Cigar> {
    let mut cigar_out = Vec::new();
    let mut last_elem = Cigar::Match(0);
    for new_elem in cigar_in {
        if std::mem::discriminant(new_elem) == std::mem::discriminant(&last_elem) {
            use Cigar::*;
            if let Match(ref mut n) | Equal(ref mut n) | Diff(ref mut n) | Del(ref mut n)
            | Ins(ref mut n) | HardClip(ref mut n) | SoftClip(ref mut n)
            | RefSkip(ref mut n) = last_elem
            {
                *n += new_elem.len();
            }
        } else {
            if last_elem != Cigar::Match(0) {
                cigar_out.push(last_elem)
            }
            last_elem = *new_elem;
        }
    }
    if last_elem != Cigar::Match(0) {
        cigar_out.push(last_elem)
    }

    cigar_out
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::{self, header, Header, HeaderView};

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
    fn test_update_ref_and_complete_read_pos() {
        let header = get_test_header();

        // Mapped read:
        let sam_line = b"qname\t0\tchr1\t10\t60\t5H5S5M5D5I5=5N5X5S\t*\t0\t0\tACGCCGTATCGTCTCGAGGACTCTAGAGCT\tDDDDDEEEEEDDDDDEEEEEDDDDDEEEEE";
        let record = bam::Record::from_sam(&header, sam_line).unwrap();

        let ref_pos_expected = [100, 100, 105, 110, 110, 115, 120, 125, 125];
        let read_pos_expected = [5, 10, 15, 15, 20, 25, 25, 30, 35];

        let mut ref_pos = 100;
        let mut read_pos = 0;
        for (index, c) in record.cigar().iter().enumerate() {
            update_ref_and_complete_read_pos(c, &mut ref_pos, &mut read_pos);
            assert_eq!(ref_pos, ref_pos_expected[index]);
            assert_eq!(read_pos, read_pos_expected[index]);
        }
        assert_eq!(record.cigar().iter().count(), ref_pos_expected.len());
    }

    #[test]
    fn test_update_hard_clipped_read_pos() {
        let header = get_test_header();

        // Mapped read:
        let sam_line = b"qname\t0\tchr1\t10\t60\t5H5S5M5D5I5=5N5X5S\t*\t0\t0\tACGCCGTATCGTCTCGAGGACTCTAGAGCT\tDDDDDEEEEEDDDDDEEEEEDDDDDEEEEE";
        let record = bam::Record::from_sam(&header, sam_line).unwrap();

        let read_pos_expected = [0, 5, 10, 10, 15, 20, 20, 25, 30];

        let mut read_pos = 0;
        for (index, c) in record.cigar().iter().enumerate() {
            update_hard_clipped_read_pos(c, &mut read_pos);
            assert_eq!(read_pos, read_pos_expected[index]);
        }
    }

    #[test]
    fn test_get_hard_clipped_read_clip_positions() {
        let cigar = record::CigarString::try_from("10H10S10M10S10H".as_bytes()).unwrap();
        let result = get_hard_clipped_read_clip_positions(&cigar);
        assert_eq!(result, (10, 20, 30));
    }

    #[test]
    fn test_get_complete_read_clip_positions() {
        let cigar = record::CigarString::try_from("10H10S10M10S10H".as_bytes()).unwrap();
        let result = get_complete_read_clip_positions(&cigar);
        assert_eq!(result, (20, 30, 50));
    }

    #[test]
    fn test_get_complete_read_hard_clip_positions() {
        let cigar = record::CigarString::try_from("10H10S10M10S10H".as_bytes()).unwrap();
        let result = get_complete_read_hard_clip_positions(&cigar);
        assert_eq!(result, (10, 40, 50));
    }

    #[test]
    fn test_is_hard_clipped() {
        let cigar = record::CigarString::try_from("10H10S10M10S10H".as_bytes()).unwrap();
        assert!(is_hard_clipped(&cigar));

        let cigar = record::CigarString::try_from("10S10M10S".as_bytes()).unwrap();
        assert!(!is_hard_clipped(&cigar));
    }

    #[test]
    fn test_compress_cigar() {
        use Cigar::*;
        let cigar = vec![
            HardClip(1),
            HardClip(1),
            SoftClip(1),
            SoftClip(1),
            Match(1),
            Match(1),
            Diff(1),
            Diff(1),
            Equal(1),
            Equal(1),
            Ins(1),
            Ins(1),
            Del(1),
            Del(1),
            Match(1),
            Match(1),
        ];
        assert_eq!(
            compress_cigar(&cigar),
            record::CigarString::try_from("2H2S2M2X2=2I2D2M")
                .unwrap()
                .to_vec()
        );
    }
}
