#/bin/python3
from math import log, exp

# Converts plot data to symmmetric log for pgfplots
#
# Call as python symlog.py FILE MIN MAX
#
# Reads FILE.dat, writes converted data to FILE_symlog.dat
# Writes TeX commands for yticks to FILE_symlog.tex
#
# Conversion:
#   10^MIN < y <  10^MAX on log scale
#  -10^MIN < y <  10^MIN on linear scale
#  -10^MAX < y < -10^MIN on inverted log scale


def symlog(x, basis):
    if x < 0:
        return -symlog(-x, basis)
    return x/basis if x < basis else log(+x/basis) + 1

def symexp(y, basis):
    if y < 0:
        return -symexp(-y, basis)
    return y*basis if y < 1 else basis*exp(+y-1)

def tikz_yline(y):
    return r"(\pgfkeysvalueof{/pgfplots/xmax}, " + str(y) + r") -- (\pgfkeysvalueof{/pgfplots/xmin}, " + str(y) + r")"

def convert_symlog(filename, logmin, logmax):
    basis = 10**logmin

    with open(f'{filename}.dat', 'r') as infile, open(f'{filename}_symlog.dat', 'w') as outfile:
        for n, line in enumerate(infile):
            if n == 0 or not line.strip():
                print(line, file=outfile)
                continue

            values = [float(val) for val in line.split('\t')]

            print(f'{values[0]}' + ''.join('\t'+str(symlog(val, basis)) for val in values[1:]), file=outfile)
        print(f"Read linear data from {infile.name}")
        print(f"Wrote symlog data to {outfile.name}")

    with open(f'{filename}_symlog.tex', 'w') as texfile:
        print(r"\def\symlogymax{" + str(symlog(+10**logmax, basis)) + "}", file=texfile)
        print(r"\def\symlogymin{" + str(symlog(-10**logmax, basis)) + "}", file=texfile)


        for name, step in [("",1), ("sparse",3)]:
            print(r"\def\symlog" + name + "ytick{" + ','.join(
                [str(symlog(-10**log, basis)) for log in reversed(range(logmin,logmax,step))] +
                ['0'] +
                [str(symlog(+10**log, basis)) for log in range(logmin,logmax,step)]
                ) + r"}", file=texfile)
            print(r"\def\symlog" + name + "yticklabels{" + ','.join(
                [f'$-10^{{{log}}}$' for log in reversed(range(logmin,logmax,step))] +
                ['$0$'] +
                [f'$+10^{{{log}}}$' for log in range(logmin,logmax,step)]
                ) + r"}", file=texfile)

        print(r"\def\symlogminorytick{" + ','.join(
            [str(symlog(mult*10**log, basis)) for log in reversed(range(logmin,logmax)) for mult in range(-9,-1)] +
            [str(mult/10) for mult in range(-9,-1)] +
            [str(mult/10) for mult in range(2,10)] +
            [str(symlog(+mult*10**log, basis)) for log in range(logmin,logmax) for mult in range(2,10)]
            ) + r"}", file=texfile)

        #
        #     print(r"\def\drawsymlog" + name + "majorgrids{%", file=texfile)
        #     print(r"    \draw[/pgfplots/every major grid] " + tikz_yline(0) + r"node[left] {$0$};", file=texfile)
        #     for log in range(logmin,logmax,step):
        #         for sign in (-1,1):
        #             print(
        #                 r"    \draw[/pgfplots/every major grid] "
        #                 + tikz_yline(symlog(sign*10**log, basis))
        #                 + (r" node[left] {$" + ('-' if sign<1 else '') + r"10^{" + str(log) + r"}$}" if not (log == logmin and name=="sparse") else '')
        #                 + ";", file=texfile)
        #     print(r"    }", file=texfile)
        #
        # print(r"\def\drawsymlogminorgrids{%", file=texfile)
        # for log in range(logmin,logmax):
        #     for mult in range(-9,10):
        #         if mult in (-1,0,1):
        #             continue
        #         print(r"    \draw[/pgfplots/every minor grid] " + tikz_yline(symlog(mult*10**log, basis)) + ";", file=texfile)
        # for mult in range(-9,10):
        #     if mult in (-1,0,1):
        #         continue
        #     print(r"    \draw[/pgfplots/every minor grid] " + tikz_yline(mult/10) + ";", file=texfile)
        # print(r"    }", file=texfile)
        print(f"Wrote symlog tex commands to {texfile.name}")

if __name__ == '__main__':
    import sys

    convert_symlog(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))
