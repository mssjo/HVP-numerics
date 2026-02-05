#!/bin/bash

set -e

file=$1

# Split input by target integral
awk -v file="$file" -f hybrid.awk $file.dat

# Take relative difference with value at reference chi

# This one uses chi=10
awk -v ref=-0.89742805932224248607174450633168257048  -f reldiff.awk ${file}_SJ1.dat > ${file}_SJ1_reldiff.dat
# These use chi=200
awk -v ref=0.0078068847596524203538465700442538148928 -f reldiff.awk ${file}_SJ2.dat > ${file}_SJ2_reldiff.dat
awk -v ref=-6.1146175849768580502608259590204633402   -f reldiff.awk ${file}_Hreg.dat > ${file}_Hreg_reldiff.dat
awk -v ref=2.8271779496429843696330956338416423182    -f reldiff.awk ${file}_E1h2.dat > ${file}_E1h2_reldiff.dat

# Convert those to symlog
for var in SJ1 SJ2 Hreg E1h2
do
    python3 symlog.py ${file}_${var}_reldiff -15 -1
done
