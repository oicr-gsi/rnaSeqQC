# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.4.0] - 2024-06-25
### Added
- [GRD-797](https://jira.oicr.on.ca/browse/GRD-797) - Add vidarr labels to outputs (changes to medata only)

## [1.3.0] - 2023-06-28
### Changed
- Moving assembly-specific settings inside the wdl

## [1.2.3] - 2022-09-01
### Changed
- Reverting back to genomic-aligned bam for bamQC.
- Running a separate bamQC on the transcriptome bam.
- Modified to ensure for bam input that both genomic and transcriptomic bam files are provided.
- Collecting bamQC output from transcriptome bam into collated json with some jq code, under bamqc_transcriptomic.
- Task for picard InsertSizeMetrics included, but not currently in use. this is an alternate way to calculate insert size metrics.


## [1.2.2] - 2022-01-28
### Changed
- Switching to transcriptome-aligned bam for bamQC (insert size calculation).

## [1.2.1] - 2021-06-01
### Changed
- Migration to Vidarr

## [1.2.0] - 2020-05-05
### Changed
- GRD-50 Making sure rnaSeqQC is in good shape
- Removed duplicated code
- Updating README, pom and test.json


## [1.1.1] - 2020-03-30
### Changed
- GP-2327 Use RNA kit strandedness to calculate strandedness values
- GP-2328 Rewrite Picard strand metrics when strandSpecificity=NONE
- Update bam-qc-metrics dependency

## [1.1.0] - 2020-02-13
### Changed
- GP-2132 Run BWA streamed, without writing fastq to disk

### Fixed
- Bugfix for `bwa mem -p`; see comment on GP-2132

## [1.0.2] - 2020-02-06
### Fixed
- GP-2270 Reduce disk usage for bamToFastq task

## [1.0.1] - 2020-02-05
### Fixed
- Fixup; pass `modules` string as a parameter to the Picard task.

## [1.0.0] - 2020-02-04
### Added
- GR-1036 Initial release of rebuilt RNASeqQC workflow.

## [0.5.0] - 2020-02-04
### Fixed
- Take reference paths instead of module names as input.

## [0.4.4] - 2020-02-03
### Fixed
- Fixup; correctly set REFFLAT_ROOT.

## [0.4.3] - 2020-02-03
### Fixed
- Fixup; correct BASH syntax again.

## [0.4.2] - 2020-02-03
### Fixed
- Fixup; correct BASH syntax; separate params for each reference module.

## [0.4.1] - 2020-02-03
### Fixed
- Fixup; modify picard command to set correct reference variables.

## [0.4.0] - 2020-02-03
### Added
- GR-1007 Add a picardModules parameter; allows reference to be set by olive.

## [0.3.1] - 2020-01-31
### Fixed
- Fixup; increment version in POM.

## [0.3.0] - 2020-01-31
### Fixed
Fixes issues with BWA alignment:
- GR-1031 Set pipefail for the bwa command; run with -p for smart pairing.
- GR-1035 Include BWA index files in the reference data module.

## [0.2.0] - 2020-01-29
### Changed
- Rename workflow from RNASeqQC to rnaSeqQC for Shesmu compatibility.


## [0.1.0] - 2020-01-29
### Added
- Development version of workflow for testing on cluster.
