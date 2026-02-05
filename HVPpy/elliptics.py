"""
File: elliptics.py
Author: Mattias Sjö

This file implements the elliptic polylogs, Eisenstein series, etc.
"""

from math import lcm
from mpmath import mpf, mpc
from mpmath import pi, inf
from mpmath import exp, eta, zeta

from .method import Method
from .integration import QuadError

def eisen(a,b, q, char=None, d_logq=0, error=True):
    """
    Eisenstein series with a character.

    Arguments:
    a - the power of 1/n in the nth term of the series
    b - the power of 1/(1 - q^n) in nth term of the series
    q - the nome
    char - the character, an array of integers (default: Eisenstein series without a character)
    d_logq - the number of derivatives to take with respect to log(q) (default: zero)
    error - if true, give the result as a QuadError (default: true)

    The series is evaluated using mpmath.nsum, and the error estimate is taken from its tolerance.
    """

    if not isinstance(q, mpc):
        q = mpc(q)

    # Implement derivatives recursively
    if d_logq < 0 or not isinstance(d_logq, int):
        raise ValueError(f"Number of log(q)q-derivatives must be a nonnegative integer, got {d_logq}")
    if d_logq > 0:
        if b == 1:
            return eisen(a-1,2, q, char, d_logq-1)
        else:
            return b*eisen(a-1,b+1, q, char, d_logq-1) + (1-b)*eisen(a-1,b, q, char, d_logq-1)

    def character(n):
        return 1 if char is None else char[int(n-1) % len(char)]
    def summand(n):
        return character(n)/n**a * q**n/(1-q**n)**b

    result = QuadError.from_nsum(summand, [1,inf])
    return result if error else result.value()

def get_character(a, coeffs):
    """
    Convert the coefficients of a linear combination of Eisenstein series.

    Arguments:
    a - the power of 1/n in the nth term of the series
    coeffs - a mapping from k to the coefficient of eisen(a,b, q**k) in the linear combination

    This allows
    sum( coeff * eisen(a,b, q**k) for k,coeff in coeffs.items() )
    to be more efficiently evaluated as
    eisen(a,b, q, get_character(a, coeffs))
    """
    period = lcm(*coeffs.keys())
    return [sum(c*k**a for k,c in coeffs.items() if (n+1)%k == 0) for n in range(period)]

# Characters
psi_t = get_character(1, {1:1, 2:-1, 3:1, 6:-1}) # Factor of 6 extracted
psi_dt = get_character(-1, {1:3, 2:-3, 3:3, 6:-3}) # Factor of 2 extracted
psi_ddt = get_character(-2, {1:18, 2:-18, 3:18, 6:-18}) # Factor of 1/3 extracted
psi_varpi = get_character(1, {1:1, 2:-2, 3:1, 6:-2}) # Factor of -2 extracted
psi_dvarpi = get_character(-1, {1:2, 2:-4, 3:2, 6:-4}) # Factor of -1 extraced
psi_ddvarpi = get_character(-2, {1:2, 2:-4, 3:2, 6:-4}) # Factor of -1 extracted
psi_bball = [1, -15, -8, -15, 1, 120]

def tau_to_t(tau, method=None, error=True):
    if not isinstance(tau, mpc):
        tau = mpc(tau)

    if method == Method.EISENSTEIN:
        logq = 2j*pi*tau
        q = exp(logq)

        logt = -logq - 6*eisen(1,1, q, psi_t)
        return -QuadError.exp(logt if error else QuadError.decay(logt))
    else:
        try:
            return -( eta(tau) * eta(3*tau) / (eta(2*tau) * eta(6*tau)) )**6
        except ValueError as err:
            raise ValueError(f"{err} ({tau=})")


# The eta implementation goes unstable close to tau=0 (q=1),
# in particular with the mpmath eta (sage's is better).
# For abs(tau) below this cutoff, use the asymptotic expansion instead
# In sage, the difference between eta and asymptotic is less than the
# numerical noise for abs(tau) < 0.03 or so
varpi_1_cutoff = 0.01
def varpi_1(tau, method=None):

    if method == Method.EISENSTEIN:
        logq = 2j*pi*tau
        q = exp(logq)

        logw = logq + 2*eisen(1,1, q, psi_varpi)
        return QuadError.exp(logw)

    elif abs(tau) > varpi_1_cutoff:
        return complex((eta(2*tau)*eta(6*tau))**4 / (eta(tau)*eta(3*tau))**2)
    else:
        return -1/(48*tau**2)
