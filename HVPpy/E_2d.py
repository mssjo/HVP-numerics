
from mpmath import inf, pi
from mpmath import exp, log, polylog, sin, cos

from .elliptics import *
from .integration import *
from .series_expansion import *
from .utilities import *

from .t_tau_beta import t_to_tau, t_to_beta, tau_to_t, tau_to_beta, beta_to_t, beta_to_tau
from .theta import *

from .bessel import bessel_integrand
from .method import Method


# The various d_logq's defined in the section on Eisenstein series
# The d_logq parameter indicates additional logq derivatives
def Dq_t(q, *, d_logq=0):
    from .elliptics import psi_t
    return (0 if d_logq else -1) - 6*eisen(1,1, q, psi_t, d_logq=d_logq+1)
def Dq_varpi_1(q, *, d_logq=0):
    from .elliptics import psi_varpi
    return (0 if d_logq else 1) + 2*eisen(1,1, q, psi_varpi, d_logq=d_logq+1)

# Elliptic polylog as defined in Bloch, Kerr & Vanhove
def Li_elliptic(r, q,logq, z):
    Li = polylog
    match r:
        case 3:
            return (Li(3,z)
                    + QuadError.from_nsum(lambda n: Li(3,q**n*z) + Li(3,q**n/z), (1,inf))
                    + log(z)**3/12 - logq*log(z)**2/24 + logq**3/720
                    )
        case 2:
            return (
                    + QuadError.from_nsum(lambda n: n*(Li(2,q**n*z) + Li(2,q**n/z)), (1,inf))
                    - log(z)**2/24 + logq**2/240
                    )
        case 1:
            return (
                    + QuadError.from_nsum(lambda n: n**2*(Li(1,q**n*z) + Li(1,q**n/z)), (1,inf))
                    + logq/120
                    )
        case 0:
            return (
                    + QuadError.from_nsum(lambda n: n**3*(Li(0,q**n*z) + Li(0,q**n/z)), (1,inf))
                    + 1/120
                    )
        case _:
            raise NotImplementedError(f"Elliptic polylogarithm of weight {r} not implemented")

# H_bball based on that in Bloch, Kerr & Vanhove
def H_bball(r, q,logq, ctx):
    match ctx.method:
        case Method.ELLIPTIC:
            zeta6 = exp(1j*pi/3)
            return (
                + 24*Li_elliptic(r, q,logq, zeta6   )
                + 21*Li_elliptic(r, q,logq, zeta6**2)
                +  8*Li_elliptic(r, q,logq, zeta6**3)
                +  7*Li_elliptic(r, q,logq, 1       )
                )

        case Method.EISENSTEIN:
            match r:
                case 3:
                    return logq**3/12 + 5*pi**2*logq/6 - zeta(3)/3 + eisen(3,1, q, psi_bball)
                case 2:
                    # return logq**2/4 + 5*pi**2/6 + eisen(2,2, q, psi_bball)
                    return logq**2/4 + 5*pi**2/6 + eisen(3,1, q, psi_bball, d_logq=1)
                case 1:
                    # return logq/2 + eisen(1,3, q, psi_bball) - eisen(1,2, q, psi_bball)
                    return logq/2 + eisen(3,1, q, psi_bball, d_logq=2)
                case 0:
                    return 1/2 + eisen(3,1, q, psi_bball, d_logq=3)
                case _:
                    raise NotImplementedError(f"H_bball of weight {r} not implemented")

        case _:
            raise NotImplementedError(f"H_bball not implemented for method {ctx.method}")

# The "H-block" is the parenthesized quantity containing H_bball
# which is common to the elliptic, Eisenstein and Poisson formulations.
def H_block(r, q,logq, ctx):

    if ctx.method == Method.POISSON:
        tau = ctx.tau

        def poisson_sum(r):
            def summand(n):
                x = exp(-1j*pi*int(n)/(3*tau))
                return int(n)**(3-r) * polylog(r, x) * sum(psi_bball[i] * cos(pi*n*(i+1)/3) for i in range(6))
            return QuadError.from_nsum(summand, [1,inf], method='alternating')

        match r:
            case 3:
                return tau**2 * (
                    - 8*poisson_sum(3)
                    - 336*zeta(3)
                    )
            case 2:
                return tau * (
                    2 * (
                        - 8*poisson_sum(3)
                        - 336*zeta(3)
                        )
                    - 1j*pi/3 * (
                        - 8*poisson_sum(2)
                        )
                    ) / (2j*pi)
            case 1:
                return (
                        (
                            2 * (
                                - 8*poisson_sum(3)
                                - 336*zeta(3)
                                )
                            # - 1j*pi/3 * (
                            #     - 8*poisson_sum(2)
                            #     )
                            )
                        - 1j*pi/3 * (
                            2*2 * ( # <- times two by absorbing omitted term above
                                - 8*poisson_sum(2)
                                )
                            - 1j*pi/3 * (
                                - 8*poisson_sum(1)
                                )
                            )
                    ) / (2j*pi)**2

            case 0:
                return -1j*pi/(3*tau) * (
                        2 * (
                            -8*poisson_sum(2)
                            )
                        - 1j*pi/3 * (
                            2*2 * (
                                - 8*poisson_sum(1)
                                )
                            - 1j*pi/3 * (
                                - 8*poisson_sum(0)
                                )
                            )
                    ) / (2j*pi)**3

            case _:
                raise ValueError(f"H-block of weight {r} not implemented")
    else:

        match r:
            case 3:
                return (40*pi**2*logq - 48*H_bball(3, q,logq, ctx))
            case 2:
                return (40*pi**2 - 48*H_bball(2, q,logq, ctx))
            case 1:
                return -48*H_bball(1, q,logq, ctx)
            case 0:
                return -48*H_bball(1, q,logq, ctx)
            case _:
                raise NotImplementedError(f"H-block of weight {r} not implemented")

# @timeout(100, float('nan'))
@log_errors
@add_QuadError
# @print_return
def E_2d(n, context=None, *, t=None, tau=None, beta=None, d_logt=0, method=Method.ELLIPTIC, **options):

    ctx = context if context else IntegrationContext(t,tau,beta, method=method, **options)

    if (n, d_logt) in ctx.cache:
        return ctx.cache[(n, d_logt)];

    # Some derivatives computed via IBP, selected with IBP_derivs
    # This is the only option available for external programs that don't
    #  support derivatives.
    if ctx.IBP_derivs:
        t = ctx.t
        match (n, d_logt):
            case (1,1):
                return -E_2d(1, ctx) - 4*E_2d(2, ctx)
            case (1,2):
                return t*E_2d(2, ctx) + 2*t*E_2d(3, ctx)
            case (1,3):
                return t*(
                        + 24
                        + 2*(t - 7)*E_2d(1, ctx)
                        - 3*(24 - 14*t + t**2)*E_2d(2, ctx)
                        + 6*(10*t - t**2)*E_2d(3, ctx)
                    )/((t - 16)*(t - 4))
            # Fallthrough otherwise

    match ctx.method:
        # double bessel only differs from bessel for E5,E6 purposes
        case Method.BESSEL | Method.DOUBLE_BESSEL:
            if d_logt:
                raise NotImplementedError(f"Not implemented: t-derivatives with Bessel")
            if not is_real(ctx.t) or Re(ctx.t) >= 16:
                return float('nan')
            return QuadError.from_quad(lambda x: bessel_integrand(n, ctx.t,x, d_logt), [(0,inf)])

        # These differ only in how the H-block (H_bball with associated factors) is computed
        case Method.ELLIPTIC | Method.EISENSTEIN | Method.POISSON:
            logq = (pi*2j*ctx.tau)
            q = exp(logq)

            w1 = varpi_1(ctx.tau, method=ctx.method)

            match (n, d_logt):
                case (1,0):
                    return -w1*H_block(3, q,logq, ctx)

                case (1,1):

                    Dqw = Dq_varpi_1(q)
                    Dqt = Dq_t(q)

                    return (
                        + Dqw * E_2d(1, ctx)
                        - w1 * H_block(2, q,logq, ctx)
                        ) / Dqt

                case (1,2):

                    Dqw = Dq_varpi_1(q)
                    Dqt = Dq_t(q)
                    Dq2w = Dq_varpi_1(q, d_logq=1)
                    Dq2t = Dq_t(q, d_logq=1)
                    E1 = E_2d(1, ctx)
                    dE1 = E_2d(1, ctx, d_logq=1)

                    return (
                        + (2*Dqw/Dqt - Dq2t/Dqt**2) * dE1
                        + (
                            + (Dq2w - Dqw**2)*E1
                            - w1*H_block(1, q,logq, ctx)
                            ) / Dqt**2
                        )

                case (1,3):

                    Dqw = Dq_varpi_1(q)
                    Dqt = Dq_t(q)
                    Dq2w = Dq_varpi_1(q, d_logq=1)
                    Dq2t = Dq_t(q, d_logq=1)
                    Dq3w = Dq_varpi_1(q, d_logq=2)
                    Dq3t = Dq_t(q, d_logq=2)
                    E1 = E_2d(1, ctx)
                    dE1 = E_2d(1, ctx, d_logq=1)
                    d2E1 = E_2d(1, ctx, d_logq=2)

                    return (
                        + (3*Dqw/Dqt - 3*Dq2w/Dqt**2) * d2E1
                        + (3*(Dq2w-Dqw**2)/Dqt**2 + (3*Dq2t*Dqw - Dq3t)/Dqt**3) * dE1
                        + (
                            (Dqw**3 - 3*Dqw*Dq2w + Dq3w) * E1
                            - w1*H_block(0, q,logq, ctx)
                            ) / Dqt**3
                        )

                case (1,_):
                    raise NotImplementedError(f"{ordinal(d_logt)} log(t)-derivative not implemented")

                case (_,d) if d != 0:
                    raise NotImplementedError(f"log(t)-derivatives not implemented for method {ctx.method}")


                case (2,0):

                    Dqw = Dq_varpi_1(q)
                    Dqt = Dq_t(q)

                    return (
                        - (Dqw/Dqt + 1)/4 * E_2d(1, ctx)
                        + w1/(4*Dqt) * H_block(2, q,logq, ctx)
                        )

                case (3,0):

                    t = ctx.t

                    Dqw = Dq_varpi_1(q)
                    Dqt = Dq_t(q)
                    Dq2w = Dq_varpi_1(q, d_logq=1)
                    Dq2t = Dq_t(q, d_logq=1)

                    Dqw_2 = Dqw * Dqw
                    Dqw_3 = Dqw * Dqw_2
                    Dqt_2 = Dqt * Dqt
                    Dqt_3 = Dqt_2 * Dqt

                    return (
                        + (
                            + 4*Dq2w*Dqt
                            - 4*Dq2t*Dqw
                            + 4*Dqt*Dqw_2
                            + Dqt_3*t
                            + Dqt_2*Dqw*t
                        ) / (8*Dqt_3 * t) * E_2d(1, ctx)
                        + (
                            + 4*Dq2t
                            - Dqt*(8*Dqw + Dqt*t)
                        ) / (8*Dqt_3 * t) * w1 * H_block(2, q,logq, ctx)
                        - w1 * H_block(1, q,logq, ctx) / (2*t * Dqt_2)
                        )

                case _:
                    raise NotImplementedError(f"Not implemented: E{n} in 2d")

        case Method.SECDEC:
            # Only import if needed
            from .psd_wrapper import pySecDec, nu_master

            if d_logt:
                if not ctx.IBP_derivs: # Fallback
                    return E_2d(n, ctx.but(IBP_derivs=True))
                raise NotImplementedError(f"Not implemented: t-derivatives with pySecDec")

            # NOTE: overall sign due to convention difference
            with pySecDec(nu_master[n], dim="2-2*eps", **ctx.options) as integral:
                return -QuadError(*integral(t=ctx.t)[0])

        case Method.AMFLOW:
            # Only import if needed
            from .amflow_wrapper import AMFlow, nu_master

            if d_logt:
                if not ctx.IBP_derivs: # Fallback
                    return E_2d(n, ctx.but(IBP_derivs=True))
                raise NotImplementedError(f"Not implemented: t-derivatives with AMFlow")

            with AMFlow(nu_master[n], dim="2", **ctx.options) as integral:
                return integral(t=ctx.t)[0]

        case Method.EXPANSION_0:
            try:
                return evaluate_series(E_2d_series[n], var=ctx.t, log_deriv=d_logt, error=+1)
            except IndexError:
                raise NotImplementedError(f"Not implemented: E{n} in 2d expanded around t=0")

        case Method.EXPANSION_INF:
            try:
                return evaluate_series(E_2d_series_inf[n], var=ctx.t, log_var=log(-ctx.t), log_deriv=d_logt, error=-1)
            except IndexError:
                raise NotImplementedError(f"Not implemented: E{n} in 2d expanded around t=∞")

        case _:
            raise NotImplementedError(f"Not implemented: E{n} in 2d using {ctx.method}")

