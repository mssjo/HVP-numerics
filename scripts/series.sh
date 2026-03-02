#!/bin/bash

set -e

order=${1:-5}

cat > series.m \
<< EOF
order   = $order;
verbose = $(if [[ $2 ]] then echo True; else echo False; fi);
debug   = $(if [[ $3 ]] then echo True; else echo False; fi);
EOF

math -script series.wl | tee series.log | while read -r line; do echo -e "$line"; done

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
        -e 's/Zeta\[3\]/zeta(3)/g'  \
        -e 's/Pi/pi/g' \
        -e 's/\^/**/g' \
        -e 's/((^|[- (])[0-9]+)\/([0-9]+)/\1\/mpf(\3)/g' \
        -e 's/$/,/; s/^ */        /' \
        >> $target
    printf "$( case $type in t) printf "\n        ]";; b) printf "        }";; esac )\n" >> $target
done


target=../HVPpy/series_expansion.py

cat > $target \
<< EOF
from mpmath import zeta, pi, mpf

# Below are tabulated series expansions,
#  computed by $(realpath $0 --relative-to=${target%/*})
#  on $(date)
#  to order $order in t (twice that in 1/beta).
# They are implemented as functions rather than dicts to enable dynamic precision.

# DO NOT MODIFY THIS FILE - CHANGES WILL BE OVERWRITTEN

#===============================================================================
# Expansions around t=0

def Ebar_series(n):
    match n:
$(for n in {1..6}; do printf "        case $n: return "; cat series/tseriesE${n}bar.py; done)
        case _: raise IndexError(str(n))

def E_2d_series(n):
    match n:
$(for n in {1..3}; do printf "        case $n: return "; cat series/tseriesE${n}.py; done)
        case _: raise IndexError(str(n))

#===============================================================================
# Expansions around beta=infinity needed in the E5bar computation

def SJ_gn_series(n):
    match n:
$(for n in {1..2}; do printf "        case $n: return "; cat series/bseriesSJg${n}.py; done)
        case _: raise IndexError(str(n))

def E1h2_series():
    return $(cat series/bseriesE1h2.py)

def Hreg_series():
    return $(cat series/bseriesHreg.py)


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

printf "\e[7m\e[32mPASSED ($(grep -c 'PASS' series.log)):\e[0m $(grep --color=never -F 'PASS' series.log | tr ':\n' ', ' | sed -E 's/PASS //g;s/, *$//')\n"
printf "\e[7m\e[31mFAILED ($(grep -c 'FAIL' series.log)):\e[0m $(grep --color=never -F 'FAIL' series.log | tr ':\n' ', ' | sed -E 's/FAIL //g;s/, *$//')\n"
