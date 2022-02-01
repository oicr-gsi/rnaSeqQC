#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

module load jq
# remove the Picard header because it includes temporary paths
for f in $(find . -xtype f);do jq 'del(.picard | .header)' $f --sort-keys | grep -v '"[[:digit:]]*": [[:digit:]]' | sed 's/\.[0-9]*//' | md5sum;done
