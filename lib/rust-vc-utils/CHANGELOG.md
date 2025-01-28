## v0.29.0 - 2024-10-31

### Added
- Add simple mean tracker
- Add new force-periodic option to `ProgressReporter`, breaks api
- Add new bam record utils

## v0.28.0 - 2024-09-04

### Changed
- Add changes to bam aux tag parsing
    - CR-380 Reduce the chance of a 'panic-in-panic' leading to a SIGILL
    - Add new parser for float aux tags

## v0.27.0 - 2024-08-29

### Changed
- Changed to github version of rust-libbigwig

## v0.26.0 - 2024-03-21

### Fixed
- Fix fasta reader error in formatted message

## v0.25.0 - 2024-03-15

### Added
- Add progress reporter

## v0.24.1 - 2023-10-05

### Changed
- Updated build to remove unused rust-htslib dependencies

## v0.24.0 - 2023-08-21

### Added
- Add new sparse window sum container

### Fixed
- Fix basemod parsing error messages

## v0.23.0 - 2023-07-26

### Added
- Add new cigar methods

### Changed
- Changed Cargo.toml to more flexibly allow dependency changes from client code, principally so that this library no
  longer imposes such an exact rust-htslib requirement.
- Multiple API breaking changes:
  - Refactored parameterless new methods into default
  - Refactored bam_util submodule names

## v0.22.0 - 2023-06-21

### Added
- Add new cigar convenience functions

### Changed
- Updated to rust-htslib 0.44.1

## v0.21.0 - 2023-06-06

### Added
- Added deterministic vector downsample
- Added new cigar processing methods.
    - Hard-clip test
    - Individual cigar segment offset functions

### Changed
- Updated cigar interfaces
  - Used `hard_clip` more consistently in all functions with an underscore. This update **intentionally breaks api**.
  - Accept `&[Cigar]` instead of `CigarString` to generalize input cigar requirements to more data types.

## v0.20.0 - 2023-05-25

### Added
- Added several new cigar processing methods.
  - This update **intentionally breaks api** for all cigar processing functions where hard-clip handling was previously
    unspecified, and any downstream bam processing functions in the library which used the previously unspecified
    behavior. In most of these cases they have been replaced with hard-clip and non-hard-clip versions. Updating code will
    have to explicitly consider which version to update to.

## v0.19.0 - 2023-05-25

### Added
- Added new cigar clipping method

### Fixed
- Reorganized internal bam utility structure

## v0.18.0 - 2023-05-25

### Added
- Added bam aux string method
- Added bam reg2bin

## v0.17.0 - 2023-05-18

### Added
- New ChromList ctor and test

## v0.16.0 - 2023-05-10

### Added
- Added new aux parse options 

## v0.15.0 - 2023-05-01

### Added
- Added bam record read to ref map

## v0.14.0 - 2023-04-13

### Added
- Added new bam aux tag utils
- Updated to 2021 rest edition

## v0.13.0 - 2023-04-07

### Added
- rev comp

## v0.12.0 - 2023-04-07

### Added
- New bam cigar tracking util

## v0.11.0 - 2023-04-04

### Fixed
- Fixed bug in `decode_cpg_meth_info` that would produce an invalid read_pos for a final C in the read with methylation when the read was reverse-mapped.

