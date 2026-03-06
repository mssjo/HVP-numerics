#!/bin/bash

set -e

var=$1
shift

case $var in
    0)  cd ..
        ./HVP bits=64 plot=masters,-8,24,2 methods=0 mso=5 $@
#         ./HVP bits=64 plot=masters,-8,24,2 methods=0 mso=8 $@
#         ./HVP bits=64 plot=masters,-8,24,2 methods=0 mso=16 $@
        ./HVP bits=64 plot=masters,-8,24,2 methods=0 mso=32 $@
        ;;

    1)  cd ..
        ./HVP bits=64 ut timeout=60 plot=masters,-8,24,2 indices=1 methods=e,s,0,f,b,l,p samples=1000_000 lambda=.1 $@
        for file in plots/prec64/E1?.dat
        do
            python scripts/symlog.py $file -15 -1
            python scripts/symlog.py $file -24 -6
        done
        ;;

    234) cd ..
        ./HVP bits=64 ut timeout=60 plot=masters,-8,24,2 indices=2,3,4 methods=e,s,0,f samples=1000_000 lambda=.1 $@
        for file in plots/prec64/E2?.dat plots/prec64/E3?.dat plots/prec64/E4?.dat
        do
            python scripts/symlog.py $file -15 -1
        done
        ;;

    56) cd ..
        ./HVP bits=64 ut timeout=60 plot=masters,-8,24,2 indices=5,6 methods=e,s,0 awc clp=.5 $@
        for file in plots/prec64/E4?.dat plots/prec64/E5?.dat plots/prec64/E6?.dat
        do
            python scripts/symlog.py $file -15 -1
        done
        ;;

    4)  E=../plots/prec64/E
        paste \
            <(./print_column "t"  ${E}1e.dat) \
            <(./print_column "Re" ${E}1e.dat "Re1") \
            <(./print_column "Im" ${E}1e.dat "Im1") \
            <(./print_column "Re" ${E}2e.dat "Re2") \
            <(./print_column "Im" ${E}2e.dat "Im2") \
            <(./print_column "Re" ${E}3e.dat "Re3") \
            <(./print_column "Im" ${E}3e.dat "Im3") \
            > ${E}1234e.dat
        ;;

    Pi) cd ..
        ./HVP bits=64 ut timeout=60 plot=Pi,-8,24,4 method=e awc clp=.5 $@
        ./HVP bits=64 ut timeout=60 plot=Pi,-8,24,4 method=0 mso=5 $@
        ./HVP bits=64 ut timeout=60 plot=Pi,-8,24,4 method=0 mso=32 $@
        ;;

    amflow) cd ../amflow
        wolframscript generate_AMFlow_data.wl
        cd ..
        python scripts/convert_AMFlow.py
        for file in plots/prec64/E?a.dat
        do
            python scripts/symlog.py $file -15 -1
        done
        ;;

    all)
        for var in 0 1 23 456 Pi amflow
        do
            ./generate_data.sh $var $@
        done
        ;;
esac


