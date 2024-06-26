CHANGELOG
=========
## [1.4.0] - 2024-06-25
- [GRD-797](https://jira.oicr.on.ca/browse/GRD-797) - add vidarr labels to outputs (changes to medata only)

## [1.3.0] - 2023-06-28
------------------
- moving assembly-specific settings inside the wdl

## [1.2.3] - 2022-09-01
------------------ 
- reverting back to genomic-aligned bam for bamQC
- running a separate bamQC on the transcriptome bam
- modified to ensure for bam input that both genomic and transcriptomic bam files are provided
- collecting bamQC output from transcriptome bam into collated json with some jq code, under bamqc_transcriptomic
- task for picard InsertSizeMetrics included, but not currently in use. this is an alternate way to calculate insert size metrics


## [1.2.2] - 2022-01-28
------------------
- Switching to transcriptome-aligned bam for bamQC (insert size calculation)

## [1.2.1] - 2021-06-01
------------------
- Migration to Vidarr

## [1.2.0] - 2020-05-05
------------------

- GRD-50 Making sure rnaSeqQC is in good shape
- Removed duplicated code
- updating README, pom and test.json


## [1.1.1] - 2020-03-30
------------------

- GP-2327 Use RNA kit strandedness to calculate strandedness values
- GP-2328 Rewrite Picard strand metrics when strandSpecificity=NONE
- Update bam-qc-metrics dependency

## [1.1.0] - 2020-02-13
------------------

- GP-2132 Run BWA streamed, without writing fastq to disk
- Bugfix for `bwa mem -p`; see comment on GP-2132

## [1.0.2] - 2020-02-06
------------------

- GP-2270 Reduce disk usage for bamToFastq task

## [1.0.1] - 2020-02-05
------------------

- Fixup; pass `modules` string as a parameter to the Picard task

## [1.0.0] - 2020-02-04
------------------

- GR-1036 Initial release of rebuilt RNASeqQC workflow

## [0.5.0] - 2020-02-04
------------------

- Take reference paths instead of module names as input

## [0.4.4] - 2020-02-03
------------------

- Fixup; correctly set REFFLAT_ROOT

## [0.4.3] - 2020-02-03
------------------

- Fixup; correct BASH syntax again

## [0.4.2] - 2020-02-03
------------------

- Fixup; correct BASH syntax; separate params for each reference module

## [0.4.1] - 2020-02-03
------------------

- Fixup; modify picard command to set correct reference variables

## [0.4.0] - 2020-02-03
------------------

- GR-1007 Add a picardModules parameter; allows reference to be set by olive

## [0.3.1] - 2020-01-31
------------------

- Fixup; increment version in POM

## [0.3.0] - 2020-01-31
------------------

Fixes issues with BWA alignment:
- GR-1031 Set pipefail for the bwa command; run with -p for smart pairing
- GR-1035 Include BWA index files in the reference data module

## [0.2.0] - 2020-01-29
------------------

- Rename workflow from RNASeqQC to rnaSeqQC for Shesmu compatibility


## [0.1.0] - 2020-01-29
------------------

- Development version of workflow for testing on cluster
