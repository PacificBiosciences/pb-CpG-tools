use rust_htslib::bam;

fn unexpected_aux_val_err(
    record: &bam::Record,
    aux_tag: &[u8],
    aux_val: bam::record::Aux<'_>,
) -> ! {
    let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
    let str_aux_tag = std::str::from_utf8(aux_tag).unwrap();
    panic!(
        "Unexpected {str_aux_tag} tag format in read {qname}: {:?}",
        aux_val,
    );
}

fn missing_aux_tag_err(record: &bam::Record, aux_tag: &[u8]) -> ! {
    let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
    let str_aux_tag = std::str::from_utf8(aux_tag).unwrap();
    panic!("Missing {str_aux_tag} tag in read {qname}",);
}

/// Retrieve an aux tag int value from bam record, if the tag exists.
///
/// In this version the tag itself is optional, but the function will still panic if the tag is present but has a non-int value
///
pub fn get_optional_int_aux_tag(record: &bam::Record, aux_tag: &[u8]) -> Option<i32> {
    match record.aux(aux_tag) {
        Ok(aux_val) => Some(match aux_val {
            bam::record::Aux::U8(val) => val as i32,
            bam::record::Aux::U16(val) => val as i32,
            bam::record::Aux::U32(val) => val as i32,
            bam::record::Aux::I8(val) => val as i32,
            bam::record::Aux::I16(val) => val as i32,
            bam::record::Aux::I32(val) => val,
            _ => unexpected_aux_val_err(record, aux_tag, aux_val),
        }),
        _ => None,
    }
}

/// Retrieve an aux tag int value from bam record
///
/// Function will panic if the tag is missing or has a non-int value
///
pub fn get_int_aux_tag(record: &bam::Record, aux_tag: &[u8]) -> i32 {
    get_optional_int_aux_tag(record, aux_tag)
        .unwrap_or_else(|| missing_aux_tag_err(record, aux_tag))
}

/// Retrieve an aux tag float value from bam record, if the tag exists.
///
/// In this version the tag itself is optional, but the function will still panic if the tag is present but has a non-float value
///
pub fn get_optional_float_aux_tag(record: &bam::Record, aux_tag: &[u8]) -> Option<f32> {
    match record.aux(aux_tag) {
        Ok(aux_val) => Some(match aux_val {
            bam::record::Aux::Float(val) => val,
            _ => unexpected_aux_val_err(record, aux_tag, aux_val),
        }),
        _ => None,
    }
}

/// Retrieve an aux tag float value from bam record
///
/// Function will panic if the tag is missing or has a non-float value
///
pub fn get_float_aux_tag(record: &bam::Record, aux_tag: &[u8]) -> f32 {
    get_optional_float_aux_tag(record, aux_tag)
        .unwrap_or_else(|| missing_aux_tag_err(record, aux_tag))
}

/// Retrieve an aux tag string value from bam record, if the tag exists.
///
/// In this version the tag itself is optional, but the function will still panic if the tag is present but has a non-string value
///
pub fn get_optional_string_aux_tag(record: &bam::Record, aux_tag: &[u8]) -> Option<String> {
    match record.aux(aux_tag) {
        Ok(aux_val) => Some(match aux_val {
            bam::record::Aux::String(val) => val.to_string(),
            _ => unexpected_aux_val_err(record, aux_tag, aux_val),
        }),
        _ => None,
    }
}

/// Retrieve an aux tag string value from bam record
///
/// Function will panic if the tag is missing or has a non-string value
///
pub fn get_string_aux_tag(record: &bam::Record, aux_tag: &[u8]) -> String {
    get_optional_string_aux_tag(record, aux_tag)
        .unwrap_or_else(|| missing_aux_tag_err(record, aux_tag))
}

pub fn is_aux_tag_found(record: &bam::Record, aux_tag: &[u8]) -> bool {
    record.aux(aux_tag).is_ok()
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
    fn test_is_aux_tag_found() {
        let header = get_test_header();

        let sam_line =
            b"qname\t4\t*\t0\t255\t*\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE\tHP:i:1";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();

        assert!(is_aux_tag_found(&rec, b"HP"));
        assert!(!is_aux_tag_found(&rec, b"HX"));
    }

    #[test]
    fn test_get_int_aux_tag() {
        let header = get_test_header();

        // Expected
        let sam_line =
            b"qname\t4\t*\t0\t255\t*\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE\tHP:i:1";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let val = get_int_aux_tag(&rec, b"HP");
        assert_eq!(val, 1);

        // Wrong format
        let sam_line =
                b"qname\t4\t*\t0\t255\t*\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE\tHP:f:0.001";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let result = std::panic::catch_unwind(|| get_int_aux_tag(&rec, b"HP"));
        assert!(result.is_err());

        // Missing tag:
        let sam_line =
            b"qname\t4\t*\t0\t255\t*\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let result = std::panic::catch_unwind(|| get_int_aux_tag(&rec, b"HP"));
        assert!(result.is_err());
    }

    #[test]
    fn test_get_optional_int_aux_tag() {
        let header = get_test_header();

        // Expected
        let sam_line =
            b"qname\t4\t*\t0\t255\t*\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE\tHP:i:1";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let val = get_optional_int_aux_tag(&rec, b"HP");
        assert_eq!(val, Some(1));

        // Wrong format
        let sam_line =
            b"qname\t4\t*\t0\t255\t*\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE\tHP:f:0.001";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let result = std::panic::catch_unwind(|| get_optional_int_aux_tag(&rec, b"HP"));
        assert!(result.is_err());

        // Missing tag:
        let sam_line =
            b"qname\t4\t*\t0\t255\t*\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let result = get_optional_int_aux_tag(&rec, b"HP");
        assert!(result.is_none());
    }

    #[test]
    fn test_get_float_aux_tag() {
        let header = get_test_header();

        // Expected
        let sam_line =
            b"qname\t4\t*\t0\t255\t*\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE\trq:f:0.999";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let val = get_float_aux_tag(&rec, b"rq");
        approx::assert_ulps_eq!(val, 0.999f32, max_ulps = 4);

        // Wrong format
        let sam_line =
            b"qname\t4\t*\t0\t255\t*\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE\trq:i:1";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let result = std::panic::catch_unwind(|| get_float_aux_tag(&rec, b"rq"));
        assert!(result.is_err());

        // Missing tag:
        let sam_line =
            b"qname\t4\t*\t0\t255\t*\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let result = std::panic::catch_unwind(|| get_float_aux_tag(&rec, b"rq"));
        assert!(result.is_err());
    }

    #[test]
    fn test_get_optional_float_aux_tag() {
        let header = get_test_header();

        // Expected
        let sam_line =
            b"qname\t4\t*\t0\t255\t*\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE\trq:f:0.999";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let val = get_optional_float_aux_tag(&rec, b"rq");
        assert!(val.is_some());
        approx::assert_ulps_eq!(val.unwrap(), 0.999f32, max_ulps = 4);

        // Wrong format
        let sam_line =
            b"qname\t4\t*\t0\t255\t*\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE\trq:i:1";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let result = std::panic::catch_unwind(|| get_optional_float_aux_tag(&rec, b"rq"));
        assert!(result.is_err());

        // Missing tag:
        let sam_line =
            b"qname\t4\t*\t0\t255\t*\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let result = get_optional_float_aux_tag(&rec, b"rq");
        assert!(result.is_none());
    }

    #[test]
    fn test_get_string_aux_tag() {
        let header = get_test_header();

        // Expected
        let sam_line =
            b"qname\t4\t*\t0\t255\t*\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE\tSA:Z:FOO";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let val = get_string_aux_tag(&rec, b"SA");
        assert_eq!(val, "FOO");

        // Wrong format
        let sam_line =
            b"qname\t4\t*\t0\t255\t*\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE\tSA:i:1";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let result = std::panic::catch_unwind(|| get_string_aux_tag(&rec, b"SA"));
        assert!(result.is_err());

        // Missing tag:
        let sam_line =
            b"qname\t4\t*\t0\t255\t*\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let result = std::panic::catch_unwind(|| get_int_aux_tag(&rec, b"SA"));
        assert!(result.is_err());
    }

    #[test]
    fn test_get_optional_string_aux_tag() {
        let header = get_test_header();

        // Expected
        let sam_line =
        b"qname\t4\t*\t0\t255\t*\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE\tSA:Z:FOO";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let val = get_optional_string_aux_tag(&rec, b"SA");
        assert_eq!(val, Some("FOO".to_string()));

        // Wrong format
        let sam_line =
            b"qname\t4\t*\t0\t255\t*\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE\tSA:i:1";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let result = std::panic::catch_unwind(|| get_optional_string_aux_tag(&rec, b"SA"));
        assert!(result.is_err());

        // Missing tag:
        let sam_line =
            b"qname\t4\t*\t0\t255\t*\t*\t0\t0\tACGCCGTATCGTCTCGAGGA\tDDDDDEEEEEDDDDDEEEEE";
        let rec = bam::Record::from_sam(&header, sam_line).unwrap();
        let result = get_optional_string_aux_tag(&rec, b"SA");
        assert!(result.is_none());
    }
}
