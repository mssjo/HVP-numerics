# Computes the large-t expansion of E1 using Pierre Vanhove's expansion
#  https://github.com/pierrevanhove/AllLoopSunset
# whose path should be entered as the variable below
all_loop_sunset = "../all_loop_sunset"

load(f"{all_loop_sunset}/routines.sage")
load(f"{all_loop_sunset}/number_theory_routines.sage")
load(f"{all_loop_sunset}/system_solver.sage")
load(f"{all_loop_sunset}/multiloop_sunset_periods.sage")
load(f"{all_loop_sunset}/multiloop_sunset_mellin.sage")
load(f"{all_loop_sunset}/multiloop_sunset_equal_masses.sage")

from ore_algebra import *
from sys import argv

loop = 3
order = int(argv[1])

OA, t, Dt = DifferentialOperators(QQ, 't')
E1_intL = expand(int_sunset_L(loop, t, [1]*(loop+1), order, simplify=True))

print(E1_intL)
