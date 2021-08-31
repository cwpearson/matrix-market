#! /bin/bash

set -eou pipefail

for m in /vscratch1/cwpears/*.mtx; do
  p="${m%.*}";
  p="$p.ppm"
  #echo $p;
  ~cwpears/repos/matrix-market/build/tools/mtx-to-ppm $m $p 2048
done
