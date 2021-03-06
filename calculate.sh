#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

module load jq
# remove the Picard header because it includes temporary paths
find . -xtype f -exec jq 'del(.picard | .header)' {} \; | python3 -mjson.tool --sort-keys
