#! /bin/bash

set -eou pipefail

~cwpears/repos/matrix-market/build/tools/mtx-stats /vscratch1/cwpears/*.mtx | tee /vscratch1/cwpears/stats.csv
