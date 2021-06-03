#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

module load jq
# remove the Picard header because it includes temporary paths
find . -type f -exec jq 'del(.picard | .header)' {} \; | python3 -mjson.tool --sort-keys
