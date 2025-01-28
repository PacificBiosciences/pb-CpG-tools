use rust_htslib::bam::{self, Header, HeaderView, Read};

pub struct BasemodProgramInfo {
    pub program: String,
    pub version: String,
}

/// Parse details on which program(s) generated basemods in the bam input
///
pub fn detect_bam_header_basemod_info(bam_header: &HeaderView) -> Vec<BasemodProgramInfo> {
    let header_hashmap = Header::from_template(bam_header).to_hashmap();
    let header_programs = header_hashmap.get("PG").unwrap();

    let mut basemod_info = Vec::new();
    for header_program in header_programs {
        let pn = match header_program.get("PN") {
            Some(x) => x,
            None => {
                // In general we check PN for the program name, but there is an exception here for primrose,
                // for which some versions enter the primrose label in ID but not PN, with the ID entry
                // potentially being a prefix:
                match header_program.get("ID") {
                    Some(x) if x.starts_with("primrose") => x,
                    _ => continue,
                }
            }
        };

        let program = if pn.starts_with("jasmine") {
            Some("jasmine".to_string())
        } else if pn.starts_with("primrose") {
            Some("primrose".to_string())
        } else {
            None
        };

        if let Some(program) = program {
            let version = match header_program.get("VN") {
                Some(x) => x,
                None => "",
            }
            .to_string();

            basemod_info.push(BasemodProgramInfo { program, version });
        }
    }

    basemod_info
}

/// Parse details on which program generated basemods in the bam input
///
pub fn detect_bam_basemod_info(bam_filename: &str) -> Vec<BasemodProgramInfo> {
    let bam_reader = bam::Reader::from_path(bam_filename).unwrap();
    detect_bam_header_basemod_info(bam_reader.header())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_bam_header_model_version() {
        // jasmine
        {
            let mut header = Header::new();
            header.push_record(
                bam::header::HeaderRecord::new(b"PG")
                    .push_tag(b"ID", "pbmm2")
                    .push_tag(b"PN", "pbmm2")
                    .push_tag(b"VN", "1.13.1"),
            );
            header.push_record(
                bam::header::HeaderRecord::new(b"PG")
                    .push_tag(b"ID", "jasmine")
                    .push_tag(b"PN", "jasmine")
                    .push_tag(b"VN", "2.2.2"),
            );
            let header_view = HeaderView::from_header(&header);

            let basemod_program_info = detect_bam_header_basemod_info(&header_view);
            assert_eq!(basemod_program_info.len(), 1);
            assert_eq!(&basemod_program_info[0].program, "jasmine");
            assert_eq!(&basemod_program_info[0].version, "2.2.2");
        }

        // primrose
        {
            let mut header = Header::new();
            header.push_record(
                bam::header::HeaderRecord::new(b"PG")
                    .push_tag(b"ID", "primrose")
                    .push_tag(b"PN", "primrose")
                    .push_tag(b"VN", "1.3.99"),
            );
            let header_view = HeaderView::from_header(&header);

            let basemod_program_info = detect_bam_header_basemod_info(&header_view);
            assert_eq!(basemod_program_info.len(), 1);
            assert_eq!(&basemod_program_info[0].program, "primrose");
            assert_eq!(&basemod_program_info[0].version, "1.3.99");
        }

        // no meth caller version
        {
            let mut header = Header::new();
            header.push_record(
                bam::header::HeaderRecord::new(b"PG")
                    .push_tag(b"ID", "pbmm2")
                    .push_tag(b"PN", "pbmm2")
                    .push_tag(b"VN", "1.13.1"),
            );
            let header_view = HeaderView::from_header(&header);

            let basemod_program_info = detect_bam_header_basemod_info(&header_view);
            assert!(basemod_program_info.is_empty());
        }

        // multiple meth caller versions
        {
            let mut header = Header::new();
            header.push_record(
                bam::header::HeaderRecord::new(b"PG")
                    .push_tag(b"ID", "jasmine")
                    .push_tag(b"PN", "jasmine")
                    .push_tag(b"VN", "2.2.0"),
            );
            header.push_record(
                bam::header::HeaderRecord::new(b"PG")
                    .push_tag(b"ID", "jasmine")
                    .push_tag(b"PN", "jasmine")
                    .push_tag(b"VN", "2.3.0"),
            );
            let header_view = HeaderView::from_header(&header);

            let basemod_program_info = detect_bam_header_basemod_info(&header_view);
            assert_eq!(basemod_program_info.len(), 2);
            assert_eq!(&basemod_program_info[0].program, "jasmine");
            assert_eq!(&basemod_program_info[0].version, "2.2.0");
            assert_eq!(&basemod_program_info[1].program, "jasmine");
            assert_eq!(&basemod_program_info[1].version, "2.3.0");
        }

        // Another case in the wild:
        // @PG     ID:primrose-6B67824C    VN:1.3.0 (commit v1.3.0)        CL:primrose --min-passes 3 -j 8 --keep-kinetics temp/PS00075_1/ccs.736-of-736.bam temp/PS00075_1/primrose.736-of-736.bam
        {
            let mut header = Header::new();
            header.push_record(
                bam::header::HeaderRecord::new(b"PG")
                    .push_tag(b"ID", "primrose-6B67824C")
                    .push_tag(b"VN", "1.3.0 (commit v1.3.0)"),
            );
            let header_view = HeaderView::from_header(&header);

            let basemod_program_info = detect_bam_header_basemod_info(&header_view);
            assert_eq!(&basemod_program_info[0].program, "primrose");
            assert_eq!(&basemod_program_info[0].version, "1.3.0 (commit v1.3.0)");
        }
    }
}
