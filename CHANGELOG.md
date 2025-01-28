## v3.0.0 - 2025-01-27

### Added
- CR-412 Add run metadata header to all bed output files

### Changed
- CR-413 Compress and index all bed file outputs
- CR-411 Build default tflite model directly into binary
  - Users no longer need to provide a model file for the default model pileup mode
- CR-407 Release source code, update code to build in standard linux context

### Fixed
- Improve error message when BAM index is missing

## v2.3.2 - 2023-08-21

### Fixed
- Improve error message when read length is too short for MM tag

## v2.3.1 - 2023-05-10

### Fixed
- Fix HP tag parse
  - The HP tag parser was not handling all possible binary int representations in BAM format
  - Fixes #47

## v2.3.0 - 2023-05-02

### Added
- Support alignment input in CRAM format

## v2.2.0 - 2023-04-14

### Added
- Added progress bar to bam processing step

### Fixed
- Clarified error message when an unmapped alignment file is input

## v2.1.1 - 2023-04-04

### Fixed
- Fix parse error in MM/ML methylation tags
  - An infrequent error occurred when a reverse-mapped read ends (in sequenced orientation) with a C base annotated by the MM tag
  - Fixes #41

## v2.1.0 - 2023-03-27

Initial update to the new binary version of `aligned_bam_to_cpg_scores`, replacing the similarly named python script
with a faster solution that should be simpler to install and run on linux systems. The new version has
near feature parity with the previous python script. See updated [README](README.md) contents for details.

## v1.2.0 - 2023-03-23

Final released state of the python implementation prior to new binary v2 release. Includes various minor bug fixes and
runtime improvements compared to v1.1.0.

## v1.1.0 - 2022-06-21

Includes several new features, such as denovo modsites mode, and addresses previously reported bugs

## v1.0.0 - 2022-03-25

Initial release of the python implementation.
