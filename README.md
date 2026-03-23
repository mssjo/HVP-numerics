# HVP-numerics

[![arXiv:2603.15252](https://img.shields.io/badge/arXiv-2603.15252-b31b1b.svg)](https://arxiv.org/abs/2603.15252?logo=arxiv)

[![Author: Mattias Sjö](https://img.shields.io/badge/author-Mattias%20Sj%C3%B6-blue)](https://inspirehep.net/authors/2747078)
[![Author: Pierre Vanhove](https://img.shields.io/badge/author-Pierre_Vanhove-blue)](https://pierrevanhove.github.io)
[![Author: Alessandro Lupo](https://img.shields.io/badge/author-Alessandro_Lupo-blue)](https://inspirehep.net/authors/1982690)

![Language: python](https://img.shields.io/badge/Language-python-yellow?logo=python)

This is the implementation behind [arXiv:2603.15252](https://arxiv.org/abs/2603.15252)
It is mainly written in Python using mpmath, with some auxiliary scripts for preparing its output for inclusion in the paper.

## Dependencies

Python 3.10 or later (not tested on anything older than 3.12.3)

[pySecDec](https://github.com/secdec-research/secdec), installed as a Python package.

[AMFlow](https://gitlab.com/multiloop-pku/amflow), installed as a Mathematica package.
AMFlow can use a variety of different backends, of which the [FIRE](https://gitlab.com/feynmanintegrals/fire)+[LiteRed 2](https://github.com/rnlg/LiteRed2) combination is used here, but this can be changed in `amflow/AMFlow_defs.m`. The script `makeShortcut.m` included in `LiteRed2` can also be used to install AMFlow.

[FeynTrop](https://github.com/michibo/feyntrop), with its Python interface installed as a package.

[SageMath](https://www.sagemath.org/), but only for some plotting functionality.
(This was originally a Sage program, but was converted to bare Python to avoid the overhead. If redone from scratch, something even faster would be used.)

## Installation

Run `scripts/series.sh N` with `N` the desired order of the series expansions; this will regenerate `HVPpy/series_expansion.py` accordingly.
Install the Python module `HVPpy`, e.g. by running `pip --install .` in the main directory.

## Usage

The main functionality is provided by the Python module `HVPpy`.
Its interface is provided by `main.py`; running `python main.py help` prints out instructions for it.
The executable `HVP` is an alias for `python main.py`.
Instructions for the generation of plots is found in `scripts/generate_plots.sh` and `scripts/generate_data.sh`.
