# HVP-numerics

[![Author: Mattias Sjö](https://img.shields.io/badge/author-Mattias%20Sj%C3%B6-blue)](https://inspirehep.net/authors/2747078)
[![Author: Pierre Vanhove](https://img.shields.io/badge/author-Pierre_Vanhove-blue)](https://pierrevanhove.github.io)
[![Author: Alessandro Lupo](https://img.shields.io/badge/author-Alessandro_Lupo-blue)](https://inspirehep.net/authors/1982690)

This is the implementation behind `PAPER TBA`.
It is mainly written in Python using mpmath, with some auxiliary scripts for preparing its output for inclusion in the paper.

## Dependencies

Python 3.10 or later (not tested on anything older than 3.12.3)

[pySecDec](https://github.com/secdec-research/secdec), installed as a Python package

[AMFlow](https://gitlab.com/multiloop-pku/amflow), installed as a Mathematica package.
AMFlow can use a variety of different backends, of which the [FIRE](https://gitlab.com/feynmanintegrals/fire)+[LiteRed 2](https://github.com/rnlg/LiteRed2) combination is used here, but this can be changed in `amflow/AMFlow_defs.m`. The script `makeShortcut.m` included in `LiteRed2` can also be used to install AMFlow.

[SageMath](https://www.sagemath.org/), but only for some plotting functionality.
(This was originally a Sage program, but was converted to bare Python to avoid the overhead. If redone from scratch, something even faster would be used.)

## Usage

The main functionality is provided by the Python module `HVPpy`, which is installed with `pip --install .` or similar in the top-level directory.
Its interface is provided by `main.py`; running `python main.py help` prints out instructions for it.
Instructions for the generation of plots is found in `scprits/generate_plots.sh`.
