#!/bin/bash

# Arguments: order margin verbose debug

set -e

# Series are computed to t^order but output to t^(order-margin)
# to account for inaccuracies in highest terms
order=${1:-5}
margin=${2:-0}

# Setup the parameters for the Mathematica program
cat > series.m \
<< EOF
order   = $order;
verbose = $(if [[ $3 ]] then echo True; else echo False; fi);
debug   = $(if [[ $4 ]] then echo True; else echo False; fi);
EOF

# Run Pierre's large-t series if needed
# inffile=series/iseriesE1_$order.m
# if [[ ! -f $inffile ]]
# then
#     echo -e "\033[90mGenerating large-t expansion with SageMath...\033[34m"
#     sage infseries.sage $order | tee >(tail -n 1 | sed -E 's/log\(-1\/t\)/Log[-1\/t]/g; s/zeta\(3\)/Zeta[3]/g' > $inffile)
#     echo -e "\033[90mWrote to $inffile\033[0m"
# fi

# Run the Mathematica program to generate the series
IFS=''
math -script series.wl | tee series.log | while read -r line; do echo -e "$line"; done

# Get longest length of any filename
maxlen=$(ls series/*.txt | wc -L)

# Convert series.wl output to Mathematica
for file in series/*.txt
do
    # Skip shift files
    if [[ $file == *shift* ]]
    then
        continue
    fi

    target=${file%.txt}.py
    type=$(printf %.1s ${file#series/})

    # Read shift if present
    shiftfile=${file%.txt}_shift.txt
    if [[ -f $shiftfile ]]
    then
        shift=$(cat $shiftfile)
    else
        shift=0
    fi

    printf "%-*s -> %-*s(%s)\n" $maxlen $file $maxlen $target $type

    # Convert file, dropping the last margin terms
    # Use nl to label the series powers if nontrivial
    printf "$( case $type in t) printf "[";; b | i) printf "{";; esac ) \n" > $target
    head -n -$margin $file \
        | case $type in
            t) tee;;
            b | i) nl -s': ' -v$shift -i-1;;
        esac \
        | sed -E \
        -e '/^ *-?[0-9]+: 0$/d' \
        -e 's/Zeta\[3\]/zeta(3)/g'  \
        -e 's/Pi/pi/g' \
        -e 's/\^/**/g' \
        -e 's/((^|[- ({])[0-9]+)\/([0-9]+)/\1\/mpf(\3)/g' \
        -e 's/\{/[/g; s/\}/]/g' \
        -e 's/$/,/; s/^ */        /' \
        >> $target
    printf "$( case $type in t) printf "\n        ]";; b | i) printf "        }";; esac )\n" >> $target
done

# Paste all into python file
target=../HVPpy/series_expansion_$((order-margin)).py
cat > $target \
<< EOF
from mpmath import zeta, pi, mpf

# Below are tabulated series expansions,
#  computed by $(realpath $0 --relative-to=${target%/*})
#  on $(date)
#  to order $order in t and 1/t (twice that in 1/beta).
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
# Expansions around t=infinity

def E_2d_series_inf(n):
    match n:
$(false && for n in {1..3}; do printf "        case $n: return "; cat series/iseriesE${n}.py; done)
        case _: raise IndexError(str(n))

def Ebar_series_inf(n):
    match n:
$(false && for n in {1..4}; do printf "        case $n: return "; cat series/iseriesE${n}bar.py; done)
        case _: raise IndexError(str(n))

EOF




echo -e "\033[90mWrote to $target\033[0m"

# Point the main series expansion file to it
ln -sf $target ../HVPpy/series_expansion.py

printf "\e[7m\e[32mPASSED ($(grep -c 'PASS' series.log)):\e[0m $(grep --color=never -F 'PASS' series.log | tr ':\n' ', ' | sed -E 's/PASS //g;s/, *$//')\n"
printf "\e[7m\e[31mFAILED ($(grep -c 'FAIL' series.log)):\e[0m $(grep --color=never -F 'FAIL' series.log | tr ':\n' ', ' | sed -E 's/FAIL //g;s/, *$//')\n"
