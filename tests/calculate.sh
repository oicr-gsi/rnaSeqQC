#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

module load jq
# remove the Picard header because it includes temporary paths
#for f in $(find . -xtype f);do jq 'del(.picard | .header)' $f --sort-keys | md5sum;done

for f in $(find . -xtype f);do jq '.bamqc' $f --sort-keys | md5sum;done
for f in $(find . -xtype f);do jq '.bamqc_transcriptome."insert size average"' $f | md5sum;done
#for f in $(find . -xtype f);do jq '.picard' $f --sort-keys | md5sum;done
for f in $(find . -xtype f);do jq '.rrna_contamination' $f --sort-keys | md5sum;done
for f in $(find . -xtype f);do jq '.unique_reads' $f --sort-keys | md5sum;done




