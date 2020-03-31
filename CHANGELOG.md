CHANGELOG
=========

v1.1.1: 2020-03-30
------------------

- GP-2327 Use RNA kit strandedness to calculate strandedness values
- GP-2328 Rewrite Picard strand metrics when strandSpecificity=NONE
- Update bam-qc-metrics dependency

v1.1.0: 2020-02-13
------------------

- GP-2132 Run BWA streamed, without writing fastq to disk
- Bugfix for `bwa mem -p`; see comment on GP-2132

v1.0.2: 2020-02-06
------------------

- GP-2270 Reduce disk usage for bamToFastq task

v1.0.1: 2020-02-05
------------------

- Fixup; pass `modules` string as a parameter to the Picard task

v1.0.0: 2020-02-04
------------------

- GR-1036 Initial release of rebuilt RNASeqQC workflow

v0.5.0: 2020-02-04
------------------

- Take reference paths instead of module names as input

v0.4.4: 2020-02-03
------------------

- Fixup; correctly set REFFLAT_ROOT

v0.4.3: 2020-02-03
------------------

- Fixup; correct BASH syntax again

v0.4.2: 2020-02-03
------------------

- Fixup; correct BASH syntax; separate params for each reference module

v0.4.1: 2020-02-03
------------------

- Fixup; modify picard command to set correct reference variables

v0.4.0: 2020-02-03
------------------

- GR-1007 Add a picardModules parameter; allows reference to be set by olive

v0.3.1: 2020-01-31
------------------

- Fixup; increment version in POM

v0.3.0: 2020-01-31
------------------

Fixes issues with BWA alignment:
- GR-1031 Set pipefail for the bwa command; run with -p for smart pairing
- GR-1035 Include BWA index files in the reference data module

v0.2.0: 2020-01-29
------------------

- Rename workflow from RNASeqQC to rnaSeqQC for Shesmu compatibility


v0.1.0: 2020-01-29
------------------

- Development version of workflow for testing on cluster
