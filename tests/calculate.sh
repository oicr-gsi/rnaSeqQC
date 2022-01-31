#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

module load jq
# remove the Picard header because it includes temporary paths
for f in $(find . -xtype f | grep -v star);do jq 'del(.picard | .header)' $f | python3 -mjson.tool --sort-keys | md5sum;done
for f in $(find . -xtype f | grep star);do jq 'del(.picard | .header)' $f | python3 -mjson.tool --sort-keys | grep -v "\"[[:digit:]]*\": [[:digit:]]" | md5sum;done
