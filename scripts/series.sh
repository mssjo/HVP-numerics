#!/bin/bash

set -e

order=${1:-5}

cat > series.m \
<< EOF
order   = $order;
verbose = $(if [[ $2 ]] then echo True; else echo False; fi);
debug   = $(if [[ $3 ]] then echo True; else echo False; fi);
EOF

math -script series.wl | while read -r line; do echo -e "$line"; done

maxlen=$(ls series/*.txt | wc -L)

for file in series/*.txt
do
    target=${file%.txt}.py
    type=$(printf %.1s ${file#series/})

    printf "%-*s -> %-*s(%s)\n" $maxlen $file $maxlen $target $type

    printf "$( case $type in t) printf "[";; b) printf "{";; esac ) \n" > $target
    cat $file \
        | case $type in
            t) tee;;
            b) nl -s': ' -v0 -i-1;;
        esac \
        | sed -E \
        -e '/^ *-?[0-9]+: 0$/d' \
        -e 's/Zeta\[3\]/zeta3/g'  \
        -e 's/Pi/pi/g' \
        -e 's/\^/**/g' \
        -e 's/((^|[- (])[0-9]+)\/([0-9]+)/\1\/mpf(\3)/g' \
        -e 's/$/,/; s/^ */        /' \
        >> $target
    printf "$( case $type in t) printf "\n        ]";; b) printf "        }";; esac ),\n" >> $target
done


target=../HVPpy/series_expansion.py

cat > $target \
<< EOF
from mpmath import zeta, pi, mpf
zeta3 = zeta(3)

# Below are tabulated series expansions,
#  computed by $(realpath $0 --relative-to=${target%/*})
#  on $(date)
#  to order $order in t (twice that in 1/beta).

# DO NOT MODIFY THIS FILE - CHANGES WILL BE OVERWRITTEN

#===============================================================================
# Expansions around t=0

Ebar_series = {
$(for n in {1..5}; do printf "    $n: "; cat series/tseriesE${n}bar.py; done)
   }

E_2d_series = {
$(for n in {1..3}; do printf "    $n: "; cat series/tseriesE${n}.py; done)
   }

#===============================================================================
# Expansions around beta=infinity needed in the E5bar computation

SJ_gn_series = {
$(for n in {1..2}; do printf "    $n: "; cat series/bseriesSJg${n}.py; done)
}

E1h2_series = $(cat series/bseriesE1h2.py)

Hreg_series = $(cat series/bseriesHreg.py)


#===============================================================================
# Expansions around t=infinity (externally computed, incomplete)

E_2d_series_inf = {
    1: {
        -1: [0, 0, 0, 4],
        -2: [48, 0, -72, 16],
        -3: [522, 432, -684, 112],
        -4: [43072/mpf(9), 6480, -7008, 1024]
        },
    2: {
        -1: [0, 0, -3],
        -2: [12, 36, -30, 4],
        -3: [153, 558, -426, 56],
        -4: [5908/mpf(3), 8364, -6024, 768]
        },
    3: {
        -1: [0, 0, 3/mpf(2)],
        -2: [-6, -6, 3],
        -3: [-105/mpf(2), 57, -27, 4],
        -4: [-1847/mpf(3), 2202, -1074, 120],
        }
    }

Ebar_series_inf = {}
EOF

echo -e "\033[90mWrote to $target\033[0m"
