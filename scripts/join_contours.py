#!/bin/python3

# Isolines are generated in arbitrarily separated segments by e.g. GNUplot.
# Those lines are plotted correctly as continuous lines, but the segment gaps
# lead to artifacts when plotted with more sophisticated linestyles or filled.
# This identifies those gaps and joins the segments into continous lines.

# The arguments should be a list of x y .dat files (additional columns are dropped).
# Each is replaced by a new x y .dat file with segments joined.
# A .log file is generated with the coordinates of all removed gaps,
# for visualization/debugging.


import sys

def is_close(p1, p2, tol=.01, log=None):
    x1,y1 = p1
    x2,y2 = p2

    if abs(x1-x2) < tol and abs(y1-y2) < tol:
        if log is not None:
            print(f'{x1} {y1}', file=log)
        return True
    return False

for filename in sys.argv[1:]:
    seg = []
    with open(filename, 'r') as r:
        s = []
        for line in r:
            if not line.strip() or line[0] == '#':
                if s:
                    seg.append(s)
                    s = []
                continue

            coords = line.strip().split()

            try:
                s.append( (float(coords[0]), float(coords[1])) )
            except ValueError:
                raise ValueError(f"Expected x y [...], got '{line.strip()}'")

    with open(f'{filename}.log', 'w') as j:
        n1 = 0
        n2 = 1
        while n2 < len(seg):
            if is_close( seg[n1][-1], seg[n2][0], log=j):
                seg[n1] += seg[n2]
                del seg[n2]
            elif is_close( seg[n1][-1], seg[n2][-1], log=j):
                seg[n1] += list(reversed(seg[n2]))
                del seg[n2]
            elif is_close( seg[n1][0], seg[n2][-1], log=j):
                seg[n1] = seg[n2] + seg[n1]
                del seg[n2]
            elif is_close( seg[n1][0], seg[n2][0], log=j):
                seg[n1] =  list(reversed(seg[n1])) + seg[n2]
                del seg[n2]
            else:
                n2 += 1

            if n2 >= len(seg):
                n1 += 1
                n2 = n1 + 1

    with open(filename, 'w') as w:
        for s in seg:
            for x,y in s:
                print(f'{x} {y}', file=w)
            print('', file=w)



