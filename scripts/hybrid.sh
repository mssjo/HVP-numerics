#!/bin/bash

set -e

file=../plots/prec$1/hybrid_err_t${2:--2.0}

for var in SJ1 SJ2 Hreg E1h2
do
    python3 symlog.py ${file}_${var} -18 -1
done
