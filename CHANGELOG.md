CHANGELOG
=========

Unreleased
----------

- GR-980 Initial release of rebuilt RNASeqQC workflow

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