"""
File: t_tau_beta.py
Author: Mattias Sjö

This file handles conversions between t, tau and beta,
including the highly nontrivial t -> tau
"""

from mpmath import pi, gamma, exp, log, sqrt, hyp2f1, mpf, inf, eta

from .utilities import *
from .elliptics import tau_to_t # for completeness
from .method import Method

# beta <-> t
def beta_to_t(beta):
    """ Convert beta to t. """
    return 4/(1-beta**2)
def t_to_beta(t):
    """ Convert t to beta, choosing the upper half-plane. """
    if t == 0:
        return inf
    elif t.real < 0 or t.real >= 4:
        return sqrt(1 - 4/t)
    elif t.imag >= 0:
        return +1j*sqrt(4/t - 1)
    else:
        return -1j*sqrt(4/t - 1)

# beta <-> tau
def beta_to_tau(beta):
    """ Convert beta to tau via t. """
    return t_to_tau(beta_to_t(beta))
def tau_to_beta(tau, method=None, error=True):
    """ Convert tau to beta via t (see tau_to_t for details). """
    from .integration import QuadError
    return t_to_beta(QuadError.decay(tau_to_t(tau, method, error)))

# misc. derivatives
def dbeta_dt(beta, n=1):
    """ The nth derivative of t_to_beta. """
    match n:
        case 0:
            return beta
        case 1:
            return (beta**2-1)**2/(8*beta)
        case 2:
            return ((beta**2-1)**3 * (1 + 3*beta**2))/(64*beta**3)
        case _:
            raise NotImplementedError(f"Not implemented: {ordinal(n)} t-derivative of beta")

def dbeta_dlogt(beta, n=1):
    """ The logarithmic derivative of t_to_beta. """
    """ The nth derivative of t_to_beta. """
    match n:
        case 0:
            return beta
        case 1:
            return (1 - beta**2)/(2*beta)
        case 2:
            return (beta**4-1)/(4*beta**3)
        case _:
            raise NotImplementedError(f"Not implemented: {ordinal(n)} log(t)-derivative of beta")
    return
def dt_dtau(t, tau, method=None, error=False):
    """ The derivative of tau_to_t (which see for details). """
    from .elliptics import eisen, psi_t
    if method == Method.EISENSTEIN:
        q = exp(2j*pi*tau)
        return -2j*pi*t* (1 + 6*eisen(0,2, q, psi_t, error=error))
    else:
        h1 = eta(1*tau)
        h2 = eta(2*tau)
        h3 = eta(3*tau)
        h6 = eta(6*tau)
        return 2j*pi * ((h1*h3)/(h2*h6))**6 * ((h1*h2)**4 + 9*(h3*h6)**4)/(h1*h2*h3*h6)
def dbeta_dtau(beta, tau, method=None, error=True):
    """ The derivative of tau_to_beta. """
    return dbeta_dt(beta) * dt_dtau(beta_to_t(beta), tau, method, error)

# t -> tau for the remainder of this file

# @verify_inverse(tau_to_t, epsabs=10*epsilon, epsrel=100*epsilon)
def t_to_tau(t, error=False):
    """
    Convert t to tau. This is expensive.

    Parameters:
    t - the t value to be converted
    error - if True, the uncertainty of the conversion is estimated based on
        how closely the inverse conversion recovers t. The result is returned
        as a QuadError.

    This method normally follows sec. 3.4 in the paper.
    When t is real, the conversion instead uses theta.py
    """

    # # Apply tiny imaginary part if needed
    # if not is_complex(t) or t.imag == 0:
    #     tsun = t_sun(t + epsilon*1j)

    # Better: use the theta method for real values
    if is_real(t, epsabs=tolerance):
        from .theta import theta_to_tau, t_to_theta
        return theta_to_tau(t_to_theta(t, tol=tolerance, error=error))

    tsun = t_sun(t)
    tau = 1j * varpi_c(tsun) / varpi_r(tsun)

    # # Check that tau ends up withing the principal inverse image of the t-plane
    # if (   tau.imag < 0
    #     or abs(tau.real) > .5
    #     or min(abs(tau - 1/6), abs(tau + 1/6)) < 1/6
    #     or min(abs(tau - 1/2), abs(tau + 1/2)) < sqrt(3)/6
    #     ):
    #     raise ValueError(f"Invalid tau value ({tau}) obtained for t={t}")

    if error:
        return QuadError(tau, abs(tau_to_t(tau) - t))
    else:
        return tau

def t_sun(t, alternative=False):
    """ Tthe two-loop sunset parameter [eq. (3.28)]. """
    return (5*t - 32 + 4*rho2_corrected(t, rho2_region(t), alternative)) / t

def j_sun(z):
    """ The function j_sunset [eq. (3.33)]. """
    return ((z-3)**3 * (z**3 - 9*z**2 + 3*z - 3)**3)/((z-9) * (z-1)**3 * z**2)
def varpi_c(z):
    """ The complex period [eq. (3.29)]. """
    return (2*pi / ((z-3)**(1/4) * rho4_corrected(z, rho4_region(z)))
                * rhoF_corrected(1728/j_sun(z), rhoF_region(z)))
def r_map(z):
    """ The variable mapping between varpi_c and varpi_r. """
    return (9*(z-1)/(z-9))
def varpi_r(z):
    """ The real period [eq. (3.29)]. """
    return 12*sqrt(3)/(z-9) * varpi_c(r_map(z))

# The below functions handle the monodromy corrections for the rho functions.
# They do this by subdividing the complex plane into regions, in each one of
# which a different sheet is used.

def rho2_region(z):
    """ Discriminate the regions for correcting rho_2 [eq. (3.31)]. """

    if z.imag < 0:
        return -rho2_region(z.conjugate())

    # Regions numbered from left to right
    if z.real < 10:
        return 1
    else:
        return 2

def rho4_region(z):
    """ Discriminate the regions for correcting rho_4 [eq. (3.32)]. """

    if z.imag < 0:
        return -rho4_region(z.conjugate())

    # "Discriminant": the region changes when this changes sign
    # Derived from the imaginary part of the polynomial under the root
    discr = -z.imag**2 + 3*z.real**2 - 18*z.real + 3

    # Regions numbered from left to right
    if discr > 0:
        if z.real < 3:
            return 1
        else:
            return 3
    else:
        return 2

def rhoF_region(z):

    if z.imag < 0:
        return -rhoF_region(z.conjugate())

    # "Discriminants": functions whose sign determines the region
    # discr_0 has the same sign as Im(j_sunset(z))
    # and is equal to Im(z)*discr_1*discr_2*[horrible third factor]
    num = 1728*(z-9)*z**2*(z-1)**3
    den = (z-3)**3*(z**3-9*z**2+3*z-3)**3
    x = z.real
    y = z.imag
    discr_0 = num.imag*den.real - num.real*den.imag
    discr_1 = (9*x - 6*x**2 + x**3 - 4*y**2 + x*y**2)
    discr_2 = (-9 + 30*x**2 - 24*x**3 + 3*x**4 + 30*y**2 - 24*x*y**2 + 2*x**2*y**2 - y**4)

    # Ensure these are done at mpf precision
    x13 = 4**(1/mpf(3));
    x23 = 4**(2/mpf(3));

    # Where lines discr(z) = 0 cross the real axis (unused ones commented out)
    # cross_0 = 3 - 2*sqrt(3)
    # cross_1 = 3 + sqrt(3)*(1 - sqrt(3)*sqrt(3 + 2*sqrt(3)))
    # cross_2 = 3 + 2*sqrt(3)
    cross_3 = 3 + 2*x13 + x23
    # cross_4 = 3 + sqrt(3)*(1 + sqrt(3)*sqrt(3 + 2*sqrt(3)))

    # Where lines discr(z) = 0 intersect (unused ones commented out)
    int_0_re = 3 - x13 - 1/2*x23
    # int_0_im = sqrt(3)*(x13 - 1/2*x23)
    int_1_re = 3 - sqrt(3)
    # int_1_im = (lambda x: -9 + 72*x - 75*x**2 + 22*x**3 - 2*x**4)(int_1_re)

    # Six unbounded regions numbered from left to right,
    # then six bounded ones numbered from left to right
    if discr_0 > 0:
        if discr_2 > 0:
            if z.real < int_0_re:
                return 1
            elif z.real < int_1_re:
                return 10
            elif z.real < cross_3:
                return 5
            elif z.imag < .25: # Exact value unknown but there is some margin
                return 12
            else:
                return 5
        elif discr_1 > 0:
            if z.real < 1:
                return 8
            else:
                assert False, "This case should not be reached"
        else:
            return 3
    else:
        if discr_2 > 0:
            if z.real < 3:
                return 2
            else:
                return 6
        elif discr_1 > 0:
            if z.real < 1:
                return 9
            elif z.real < 3:
                return 11
            else:
                return 4
        else:
            return 7

def rho2_corrected(z, region, alternative=False):
    # The opposite sign convention is also valid if root4 and F1 are adjusted
    # accordingly. However, that adjustment requires them to be multivalued
    # functions, since the resulting image of the real line wraps around their
    # branch points.

    if alternative:
        return -rho2_corrected(z, region)

    root2 = sqrt(z**2 - 20*z + 64)
    match region:
        case 1 | -1:
            return -root2
        case 2 | -2:
            return +root2

def rho4_corrected(z, region):
    root4 = (z**3 - 9*z**2 + 3*z - 3)**(1/4)

    match region:
        case 1 | 2:
            return -root4
        case -1 | -2:
            return +root4
        case 3 | -3:
            return 1j*root4

def F1(z):
    return hyp2f1(1/mpf(12),5/mpf(12), 1, z)
def F2(z):
    return hyp2f1(1/mpf(12),5/mpf(12), 1/2, 1-z)
# def F3(z):
    # return hyp2f1(7/12,11/12, 3/2, 1-z)*sqrt(z-1)

F1_at_1 = gamma(1/2)/(gamma(7/mpf(12))*gamma(11/mpf(12)))

def rhoF_corrected(z, region):
    f1 = F1(z)
    f2 = F2(z)
    ll = 2*F1_at_1

    # C1: f1 -> ll*f2 - f1, f2 unchanged
    # C0: f2 -> f2 - 1j*f1/ll, f1 unchanged

    match region:
        case 10 | -10 | 11 | -11 | 12 | -12: # C1 monodromy
            return ll*f2 - f1
        case 9:     # ... plus C0 monodromy (or C0**-1 if Im(z) < 0)
            return ll*f2 - (1+1j)*f1
        case -9:
            return ll*f2 - (1-1j)*f1
        case 8:     # ... plus another C1
            return -1j*ll*f2 + (1+1j)*f1
        case -8:
            return +1j*ll*f2 + (1-1j)*f1
        case 7:     # ... plus another C0
            return -1j*ll*f2 + 1j*f1
        case -7:
            return +1j*ll*f2 - 1j*f1
        case _:
            return f1 # Principal branch

