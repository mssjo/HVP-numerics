from enum import Enum

from mpmath import pi, inf, nan, mp, workprec
from mpmath import log, atan, zeta, arg, exp, sign

from .clogging import clogger
from .method import Method

from .integration import *
from .psd_wrapper import *
from .series_expansion import *
from .utilities import *
from .points import threshold

from .t_tau_beta import t_to_tau, t_to_beta, tau_to_t, tau_to_beta, beta_to_t, beta_to_tau
from .t_tau_beta import dbeta_dt, dbeta_dtau, dbeta_dlogt
from .theta import theta_integral, theta_to_beta, theta_to_tau, t_to_theta

from .bessel import bessel_integrand
from .E_2d import E_2d
from .Jbub import Jbub

# Shift purely real or imaginary beta into a quadrant where Im(t)>0
def beta_eps(beta, eps=None):
    if beta.real == 0 or beta.imag == 0:
        return beta + (eps if eps else mp_eps())*(1+1j if (beta.real > 0 or beta.imag > 0) else -1-1j)
    else:
        return beta

def log1(beta):
    return log((beta+1)/(beta-1))
def log2(beta):
    return log((beta**2 - 1)/beta**2)

def g1(beta):
    return beta**2 - 1
def g2(beta):
    return g1(beta)/4 * log1(beta) - beta/2
def gn(n, beta):
    match n:
        case 1: return g1(beta)
        case 2: return g2(beta)
        case _: raise ValueError(f"g{n} is not a function")

def s0(beta):
    # return (54*beta**6 - 57*beta**4 + 11*beta**2 - 6)/(3*beta**2*g1(beta)) <- mysterious erroneous version
    return 2 * (beta**2*(2*beta**2 - 1)) / (3*g1(beta)**2)
def s1(beta):
    return (2*beta/3) * (3*beta**4 - 7*beta**2 + 3)/g1(beta)
def s2(beta):
    return beta**2 * (4*beta**2 - 3)/2

def h1(xi):
    return (54*xi**6 - 57*xi**4 + 11*xi**2 - 6) / (3*xi**2 * g1(xi))
def h2(xi):
    # return (
    #         (54*xi**8 - 111*xi**6 + 68*xi**4 - 17*xi**2 + 6)*log1(xi)/2
    #         - 54*xi**7 - 93*xi**5 + 31*xi**3 + 6*xi
    #         ) / (6*xi**2 * g1(xi)**2)
    return (
        g1(xi) * (54*xi**6 - 57*xi**4 + 11*xi**2 - 6)*log1(xi)
        - 2*xi * (54*xi**6 - 93*xi**4 + 31*xi**2 + 6)
        ) / (12*xi**2 * g1(xi)**2)
def hn(n, xi):
    match n:
        case 1: return h1(xi)
        case 2: return h2(xi)
        case _: raise ValueError(f"h{n} is not a function")

def Srat(beta):
    return (
        + (68*beta**4 - 44*beta**2 + 24 - beta**2*(beta**2+5)*pi*pi)/(6*g1(beta)**2)
        + (2*beta**2*zeta(3))/(3*g1(beta))
        )

def SJ(beta):
    return (
        beta**2/(3*g1(beta)) * Jbub(3,beta)
        + 2/g1(beta) * Jbub(1,beta) * Jbub(2,beta)
        - (beta**4+5)/g1(beta)**2 * Jbub(2,beta)
        - 12/g1(beta)**2 * Jbub(1,beta)**2
        + ((pi*pi-4)*beta**4 + (8-pi*pi)*beta**2 - 52)/(2*g1(beta)**2) * Jbub(1,beta)
        )


# The function int_SX_gn is the integral of S_X(xi) g_n(xi) / xi^2
# as described in the text
@print_return
def int_Srat_gn(n, ctx):
    b = ctx.beta

    if n == 1:
        return (
            (pi*pi-8)/2 * log1(b)
            + (24 + (68 + 4*zeta(3) - pi*pi)*b**2)/(6*b)
            )
    else: # n == 2
        return (
            ((pi*pi-8)/16 * log1(b)
                + (24 + (68 + 4*zeta(3) - pi*pi)*b**2)/(24*b)) * log1(b)
            + 2*log2(b)
            + (8 - pi*pi)/(4*(b**2 - 1))
            - (68 - pi*pi + 4*zeta(3))/12 # Extra term to make c1=0
            )

DEFAULT_IJ_XOVER = 100
DEFAULT_IE_XOVER = 10

@print_return
def int_SJ_gn(n, ctx):
    IJ_xover = ctx.get("IJ_xover", DEFAULT_IJ_XOVER, mpf)
    return hybrid_integral(
        lambda xi: gn(n,xi)/xi**2 * SJ(xi), SJ_gn_series(n),
        ctx.beta, xover = IJ_xover,
        int_real = (ctx.t < 0 or n == 1),
        var_real = (ctx.t < 0),
        source = f'SJ{n}',
        **ctx.options
        )

# The divergent part of E1(beta=xi) as xi -> infty
def E1_div(xi):
    return ((7*zeta(3) - 6)/(4*xi**2) - 7*zeta(3))

# The divergent and regulated parts of int h1 E1
@print_return
def int_Hdiv(ctx):
    b = ctx.beta
    return (
        -42*b**3*zeta(3)
        + b*(77*zeta(3) - 54)/2
        + (7*zeta(3) + 2)*log1(b)/4
        + (203*zeta(3) - 30)/(b*12)
        + (6 - 7*zeta(3))/(6*b**3)) # ERROR CORRECTED 2025-11-26

# This is roughly the point where (E1 - div) becomes zero within numerical precision
# Continuing the integral above this point just has h1 amplify the residual noise, spoiling the convergence
# TODO: increase (decrease in tau) if using higher-precision numerics
# DEFAULT_HREG_XOVER = Re(t_to_beta(tau_to_t(0.05j))) # 8828
# Similarly, this is about where the tail begins to spazz out
# DEFAULT_E1h2_XOVER = 1e2
# DEFAULT_E1h2_XOVER = 1e5

@print_return
def int_Hreg(ctx):

    IE_xover = ctx.get("IE_xover", DEFAULT_IE_XOVER, mpf)

    if ctx.method == Method.DOUBLE_BESSEL:
        # FIXME this is broken by bad convergence
        # @tabulate_return
        def Hreg_integrand(xi,x):
            return h1(xi) * (bessel_integrand(1, beta_to_t(xi),x) - E1_div(xi))

        # TODO: figure out how to apply hybrid_integral here
        if t < 0:
            return -QuadError.from_quad(lambda xi,x: Hreg_integrand(xi,x),
                                        [(Re(beta), Hrex_xover), (0,inf)])
        else: # 0 < t < 4
            return -1j*QuadError.from_quad(lambda xi,x: Hreg_integrand(1j*xi,x),
                                           [(Im(beta), IE_xover), (0,inf)])
    else:

        if ctx.use_theta:
            # @tabulate_return
            # @print_return
            def Hreg_integrand(xi, tau):
                try:
                    return h1(xi) * (E_2d(1, ctx.but(beta=False, tau=tau)) - E1_div(xi))
                except ValueError as err:
                    raise ValueError(f"{err} @  tau={tau}, xi={xi}, t={beta_to_t(xi)}")

            return theta_integral(
                Hreg_integrand, Hreg_series(),
                ctx.beta, IE_xover*(1 if ctx.t < 0 else 1j),
                source = 'Hreg',
                **ctx.options
                )
        else:
            # @tabulate_return
            def Hreg_integrand(xi):
                try:
                    return h1(xi) * (E_2d(1, ctx.but(beta=xi)) - E1_div(xi))
                except ValueError as err:
                    raise ValueError(f"{err} @ xi={xi}, t={beta_to_t(xi)}, tau={beta_to_tau(xi)}")

            return hybrid_integral(
                Hreg_integrand, Hreg_series(),
                ctx.beta, xover = IE_xover, midpoints=[1e1],
                int_real = True,
                var_real = (ctx.t < 0),
                source = 'Hreg',
                **ctx.options
                )

# int_SE1_gn, but after integration by parts and without boundary terms
@print_return
def int_E1_hn(n, ctx):

    IE_xover = ctx.get("IE_xover", DEFAULT_IE_XOVER, mpf)

    if n == 1:
        return int_Hdiv(ctx) + int_Hreg(ctx)
    if n == 2:

        # Optimization of quad(E1) in the case where E1 is obtained through Bessel integration
        # by doing everything in one go as a double integrand
        if ctx.method == Method.DOUBLE_BESSEL:

            return hybrid_integral(
                lambda xi,x: h2(xi)*bessel_integrand(1, beta_to_t(xi), x), E1h2_series(),
                ctx.beta, xover=IE_xover, other_limits=[0,inf],
                int_real = (ctx.t < 0),
                var_real = (ctx.t < 0),
                source = 'E1h2',
                **ctx.options
                )

        else:
            if ctx.use_theta:
                # @tabulate_return
                def E1h2_integrand(xi, tau):
                    try:
                        return h2(xi) * E_2d(1, ctx.but(beta=False, tau=tau))
                    except ValueError as err:
                        raise ValueError(f"{err} @ xi={xi}, tau={tau}, t={beta_to_t(xi)}")

                return theta_integral(
                    E1h2_integrand, E1h2_series(),
                    ctx.beta, IE_xover*(1 if ctx.t < 0 else 1j),
                    source = 'E1h2',
                    **ctx.options
                    )

            else:
                # @tabulate_return
                def E1h2_integrand(xi):
                    try:
                        E1 = E_2d(1, ctx.but(beta=xi))
                        return h2(xi) * E1
                    except ValueError as err:
                        raise ValueError(f"{err} @, xi={xi}, t={beta_to_t(xi)}")

                return hybrid_integral(
                    E1h2_integrand, E1h2_series(),
                    ctx.beta, xover = IE_xover, midpoints=[1e1],
                    int_real = (ctx.t < 0),
                    var_real = (ctx.t < 0),
                    source = 'E1h2',
                    **ctx.options
                    )

# This is the function \mathcal I[beta; G1,G2] in the text.
# Works for beta real or purely imaginary (within tolerance)
@print_return
def E5_subthr_integral(n, t,beta, context):
    ctx = IntegrationContext(t,None,beta, context)

    Erat = int_Srat_gn(n, ctx)
    EJ   = int_SJ_gn  (n, ctx)
    EE1  = int_E1_hn  (n, ctx)

    clogger.debug(f"> Erat({n}) = {Erat}")
    clogger.debug(f"> EJ({n})   = {EJ}")
    clogger.debug(f"> EE1({n})  = {EE1}")

    return Erat + EJ + EE1

T_STAR = -2.0
BETA_STAR = sqrt(3)
TAU_STAR = 1j/(2*sqrt(3))

# T_STAR = -8.0
# BETA_STAR = sqrt(3/2)
# TAU_STAR = 1j/sqrt(6)

# T_STAR = -32.0
# BETA_STAR = 3/sqrt(8)
# TAU_STAR = 1j/sqrt(3)

assert approx_equal(beta_to_t(BETA_STAR), T_STAR)
assert approx_equal(tau_to_t(TAU_STAR), T_STAR)

# The function values are extremely expensive to compute at the star point,
#  more so than the contour integrals, and may also require higher working precision.
# Therefore, this method maintains a persistent cache
#  (a global varaible within a session, a JSON file between them)
#  storing those values.
# It always loads the highest-precision values possible, to benefit from
#  higher working precision during the one-time computation of the cache.
STAR_POINT_INTEGRALS = None
STAR_POINT_FILE = "star_point_integrals.json"
def get_star_point_integrals(ctx):
    import json

    global STAR_POINT_INTEGRALS

    def compute_star_point(ctx):
        ctx_g = ctx.but(t=T_STAR, tau=TAU_STAR, beta=BETA_STAR)
        result = {
            'J': {1: int_SJ_gn(1, ctx_g), 2: int_SJ_gn(2, ctx_g)},
            'E': {1: int_Hreg (   ctx_g), 2: int_E1_hn(2, ctx_g)}
            }
        log = clogger.debug if ctx.no_cache else clogger.info
        for x in "JE":
            for n in (1,2):
                log(f"I{x}({n}),{ctx.method.value} at star point: {result[x][n]}")

        return result

    # Bypass cache
    if ctx.no_cache:
        return compute_star_point(ctx)

    # Obtain the persistent cache if that has not already been done
    if STAR_POINT_INTEGRALS is None:
        try:
            with open(STAR_POINT_FILE, 'r') as cache:
                cached = json.load(cache)
            STAR_POINT_INTEGRALS = map_nested(cached, json_to_QE)
        except IOError:
            clogger.warning(f"JSON cache not found, creating {STAR_POINT_FILE} anew...")
            STAR_POINT_INTEGRALS = {}
        except json.decoder.JSONDecodeError:
            clogger.warning(f"JSON cache corrupted, creating {STAR_POINT_FILE} anew...")
            STAR_POINT_INTEGRALS = {}

    # Find a value in the cache with the correct method and at least the current precision
    if ctx.method.value not in STAR_POINT_INTEGRALS:
        STAR_POINT_INTEGRALS[ctx.method.value] = {}
    max_prec = max((int(p) for p in STAR_POINT_INTEGRALS[ctx.method.value].keys())) if STAR_POINT_INTEGRALS[ctx.method.value] else 0
    if max_prec < mp.prec or ctx.overwrite_cache:
        max_prec = mp.prec

        # Compute the value if needed, then update the cache
        STAR_POINT_INTEGRALS[ctx.method.value][str(max_prec)] = compute_star_point(ctx)
        with open(STAR_POINT_FILE, 'w') as cache:
            json.dump(map_nested(STAR_POINT_INTEGRALS, QE_to_json), cache)

    clogger.debug(f"Using {max_prec}-bit precision cache")

    # Return the desired value
    return {key: {int(k): v for k,v in val.items()} for key,val in STAR_POINT_INTEGRALS[ctx.method.value][str(max_prec)].items()}

# STAR_POINT_INTEGRALS = {
#     Method.EISENSTEIN: {
#             'J': (0, QuadError(-0.897428059321723, 5.00874246940289e-8), QuadError(0.00780688475965246, 3.48650160216168e-21)),
#             'E': (0, QuadError(-6.11461758292908, 7.13528816344632e-7), QuadError(2.82717794964298, 1.69269914886564e-9))
#         }
#     }

# The contour integral from the given kinematics to the star point
@print_return
def E5_contour_integral(n, ctx):

    # @tabulate_return
    def integrand_J(xi):
        return gn(n, xi)/xi**2 * SJ(xi)

    sgn_real = -1 if ctx.beta.real < 0 else +1
    sgn_imag = -1 if ctx.beta.imag < 0 else +1

    if abs(ctx.beta.real) > 1 or abs(ctx.beta.imag) + abs(ctx.beta.real)/BETA_STAR > 1:
        xi_points = (ctx.beta,
                       BETA_STAR*sgn_real)
    else:
        xi_points = (ctx.beta,
                       sgn_real+1j*sgn_imag*(1-1/BETA_STAR),
                       BETA_STAR*sgn_real)

    # Minus sign since beta is supposed to be the upper integration limit
    integral_J = -QuadError.from_line_contour(
        integrand_J, xi_points, **ctx.options)

    if ctx.no_tau_contour:
        if n == 1:
            def integrand_E(xi):
                return hn(1, xi) * (E_2d(1, ctx.but(t=None, tau=None, beta=xi)) - E1_div(xi))
        else:
            def integrand_E(xi):
                return hn(2, xi) *  E_2d(1, ctx.but(t=None, tau=None, beta=xi))

        integral_E = -QuadError.from_line_contour(
            integrand_E, xi_points, **ctx.options)
    else:
        if n == 1:
            def integrand_E(tau):
                xi = tau_to_beta(tau)
                return hn(1, xi) * (E_2d(1, ctx.but(t=None, tau=tau, beta=None)) - E1_div(xi)) * dbeta_dtau(xi, tau)
        else:
            def integrand_E(tau):
                xi = tau_to_beta(tau)
                return hn(2, xi) *  E_2d(1, ctx.but(t=None, tau=tau, beta=None)) * dbeta_dtau(xi, tau)

        integral_E = -QuadError.from_line_contour(
            integrand_E, (ctx.tau, TAU_STAR), **ctx.options)

    return integral_J, integral_E


# @timeout(100, float('nan'))
@log_errors
@add_QuadError
@print_return
# @real_in_real_out
def Ebar(n, context=None, *,
         tau=None,t=None,beta=None,
         method=Method.EISENSTEIN, div=0, d_logt=0,
         **options):

    global STAR_POINT_INTEGRALS

    ctx = context if context else IntegrationContext(t,tau,beta, method=method, **options)

    # Implement imaginary_error
    if ctx.imaginary_error and (is_real(ctx.t) and Re(ctx.t) < threshold[n]):
        result = Ebar(n, ctx.but(imaginary_error=False), div=div, d_logt=d_logt)
        result.add_error(result.imag)
        return result.real
    # Implement timeout - note that underlying calls to E_2d don't get timers
    if ctx.timeout:
        return timeout(int(ctx.timeout), nan)(Ebar)(n, ctx.but(timeout=False), div=div, d_logt=d_logt)

    # Peek values indicative of zero
    if ctx._t == 0 or ctx._tau == 0 or ctx._beta == inf:
        return Ebar_series(n)[0]

    if ctx.subtract_0 and not d_logt: # zero value vanishes under d_logt
        return Ebar(n, ctx.but(subtract_0=False), div=div, d_logt=d_logt) - Ebar_series(n)[0]

    if (n, d_logt) in ctx.cache and not ctx.no_cached:
        return ctx.cache[(n, d_logt)]

    if ctx.IBP_derivs:
        pass
        # match (n, d_logt):
            # case (5,1):
                # TODO
            # case (5,2):
                # TODO
            # Fallthrough otherwise

    if d_logt and method.needs_IBP_derivs():
        if not ctx.IBP_derivs: # Try to use the hardcoded IBP derivatives
            return Ebar(n, ctx.but(IBP_derivs=True), d_logt=d_logt, div=div)
        raise NotImplementedError(f"Not implemented: t-derivatives with {ctx.method}")

    if div < 0:
        raise ValueError(f"Positive epsilon powers not supported, got 1/eps^{div}")
    if div and ctx.method != Method.SECDEC:
        from .divergences import Ebar_divergence

        return Ebar_divergence(n, div, ctx.t, ctx.beta)

    # Some method-specific preparations/results
    match ctx.concrete_method():
        case Method.DUMMY:
            return mpf("nan")

        case Method.GRID:
            from .integration import interpolate_grid
            return interpolate_grid(
                lambda c: Ebar(n, c, div=div, d_logt=d_logt),
                f"{f'e{div}' if div else ''}{f'd{d_logt}' if d_logt else ''}E{n}bar",
                ctx.but(method=Method.AUTO))

        case Method.SECDEC:
            # Only import if needed
            from .psd_wrapper import pySecDec

            with pySecDec(n, dim="4-2*eps", **ctx.options) as integral:
                return -integral(t=ctx.t)[0]

        case Method.AMFLOW:
            # Only import if needed
            from .amflow_wrapper import AMFlow

            with AMFlow(n, dim="4",**ctx.options) as integral:
                return integral(t=ctx.t)[0]

        case Method.FEYNTROP if n in {5,6}:
            # Only import if needed
            from .feyntrop_wrapper import FeynTrop

            if not n in {1,2,3,4}:
                raise NotImplementedError(f"Finite part of divergent integrals with FeynTrop")
            # Fall-through otherwise

            # with FeynTrop(n, dim="4",**ctx.options) as integral:
                # return integral(t=ctx.t)[0]

        case Method.EXPANSION_0:
            try:
                return evaluate_series(Ebar_series(n), var=ctx.t, log_deriv=d_logt, error=True, max_order=ctx.max_series_order)
            except IndexError:
                raise NotImplementedError(f"Not implemented: Ebar{n} expanded around t=0")
        case Method.EXPANSION_INF:
            try:
                return evaluate_series(Ebar_series_inf(n), var=ctx.t,  log_deriv=d_logt, error=True, max_order=ctx.max_series_order)
            except IndexError:
                raise NotImplementedError(f"Not implemented: Ebar{n} expanded around t=inf")

        case Method.BESSEL | Method.DOUBLE_BESSEL:
            if not is_real(ctx.t) or Re(ctx.t) >= 16:
                return float('nan')

    t = ctx.t
    beta = ctx.beta

    match (n, d_logt):

        # Most of the cases are simply given in terms of the 2d integrals

        case (1, 0):
            return (
                + E_2d(1, ctx) * (t**3 + 24*t**2 - 600*t + 896)/288
                - E_2d(2, ctx) * (t**4 + 12*t**3 - 792*t**2 + 3968*t + 1536)/288
                - E_2d(3, ctx) * (t-16)*(t-4)*(t**2 + 40*t + 64)/144
                + (71*t**2 + 18*t + 6102)/216
                + (23 - t)*pi*pi/12
                - 2*zeta(3)
                )

        case (2, 0):
            return (
                + E_2d(1, ctx) * (5*t**2 - 80*t + 96)/96
                - E_2d(2, ctx) * (5*t**3 - 116*t**2 + 376*t + 896)/96
                - E_2d(3, ctx) * (t-16)*(t-4)*(5*t + 16)/48
                + (95*t + 112 + (28 - t)*pi*pi)/48
                - zeta(3)
                )

        case (3, 0):
            return (
                + E_2d(1, ctx) * (t**2 - 12*t - 8)/96
                - E_2d(2, ctx) * (t**3 - 24*t**2 + 88*t + 160)/96
                - E_2d(3, ctx) * (t-16)*(t-4)*(t+4)/48
                + (39*t - 232 - 10*pi*pi)/48
                )

        case (4, 0):
            return (
                + E_2d(1, ctx) * (t**2 - 20*t + 160)/96
                - E_2d(2, ctx) * (t**2 - 28*t + 120)*(t-16)/96
                - E_2d(3, ctx) * (t-4)*(t-16)**2/48
                + (0 if ctx.elliptic_e4 else
                    - Jbub(1, beta) * (59 + 3*pi*pi)/8
                    + Jbub(2, beta) * 17/8
                    - Jbub(3, beta) / 4
                    + (25*t - 136 + 14*pi*pi)/24
                    - zeta(3)
                        )
                )

        # E5 and E6 are done using the main approach implemented in this file

        case (5, 0):
            boundary = (
                + s2(beta)/beta**2 * E_2d(1, ctx)
                    if beta else -3/2 * E_2d(1, ctx)
                )
            G1 = g1(beta)
            G2 = g2(beta)

        case (5, 1):
            boundary = (
                + (t-16)/(2*t)                   * E_2d(1, ctx, d_logt=1)
                + (t**2 - 28*t + 48)/(3*t*(t-4)) * E_2d(1, ctx, d_logt=0)
                )

            dbdt = dbeta_dlogt(beta)

            G1 = dbdt * 2*beta
            G2 = dbdt * Jbub(1,beta)/2

        case (5, 2): #FIXME
            boundary = (
                + (t-16)/(2*t)                   * E_2d(1, ctx, d_logt=2)
                + 16/t                           * E_2d(1, ctx, d_logt=1)
                - (t**3 - 42*t**2 + 56*t + 192)/(6*t*(t-4)**2)
                + 16/t                           * E_2d(1, ctx, d_logt=0)
                - 4/(t-4)**2 * (SJ(beta) + Srat(beta))
                )

            dbdt = dbeta_dlogt(beta, n=1)
            d2bdt2 = dbeta_dlogt(beta, n=2)
            jb = Jbub(1,beta)

            G1 = (
                + dbdt**2 * 2
                + d2bdt2  * 2*beta
                )
            G2 = (
                + dbdt**2 * (jb/2 + 1/(1-beta**2))/beta
                + d2bdt2  * jb/2
                )

        # For E6, the "boundary" term also includes all the extra terms from its relation to E5
        case (6, 0):

            # Cruder method, just eq. (4.9) directly
            # return (
            #     + E_2d(1, ctx, d_logt=2) * (t-16)*(t-4)/(12*t**2)
            #     + E_2d(1, ctx, d_logt=1) * ((t-10)/(12*t))
            #     + E_2d(1, ctx, d_logt=0) * (1/24)
            #     + (
            #         - Jbub(3,beta) * 2
            #         + Jbub(2,beta) * 6
            #         + Jbub(1,beta) * 3*(4-pi*pi)
            #         ) / (24*t)
            #     + (pi*pi - 4*zeta(3) - 20) / (12*t)
            #     - Ebar(5, ctx, d_logt=0) / t
            #     - Ebar(5, ctx, d_logt=1) / t
            #     )

            # The two invocations of E5(') merged into a single integration
            # G1 cancels, leading to E6 cheaper than E5
            boundary = (
                + E_2d(1, ctx, d_logt=2) * (t-16)*(t-4)/(12*t**2)
                + E_2d(1, ctx, d_logt=1) * ((t-10)/(12*t) - (t-16)/(2*t**2))
                + E_2d(1, ctx, d_logt=0) * (1/24
                                            - (t**2 - 28*t + 48)/(3*t**2*(t-4))
                                            - s2(beta)/beta**2 / t
                                            )
                + (
                    - Jbub(3,beta) * 2
                    + Jbub(2,beta) * 6
                    + Jbub(1,beta) * 3*(4-pi*pi)
                    ) / (24*t)
                + (pi*pi - 4*zeta(3) - 20) / (12*t)
                )

            # Corresponding to -(d/dt + 1/t)*E5
            dbdt = dbeta_dt(beta)
            G1 = -(dbdt*2*beta         + g1(beta)/t)
            G2 = -(dbdt*Jbub(1,beta)/2 + g2(beta)/t)

            if not approx_equal(G1, 0.):
                clogger.warn(f"G1 should be zero for E6bar, but got {G1} at t={t}")
            G1 = QuadError(0.)

        case _:
            raise NotImplementedError(f"{f'{ordinal(d_logt)} log(t)-derivative of ' if d_logt else ''}Ebar{n} not implemented")

    clogger.debug(f"E{n} boundary term = {boundary}")

    # The method that works for t=(∞,4] (sec. 4.1)
    # We may choose to still use the contour method for t=[0,4] or everywhere
    if (is_real(t)
        and Re(t) < (0 if ctx.timelike_contour else 4)
        and not ctx.all_with_contour):

        clogger.debug(f"Performing real integrals ({ctx.timelike_with_contour=}, {ctx.all_with_contour=})")

        ctx_s = ctx.but(t=Re(t))

        Erat = (0, int_Srat_gn(1, ctx_s), int_Srat_gn(2, ctx_s) if n != 6 else 0)
        EJ   = (0, int_SJ_gn  (1, ctx_s), int_SJ_gn  (2, ctx_s) if n != 6 else 0)
        EE1  = (0, int_E1_hn  (1, ctx_s), int_E1_hn  (2, ctx_s) if n != 6 else 0)

        clogger.debug(f"> Part 1:")
        clogger.debug(f"> > Erat = {Erat[1]}")
        clogger.debug(f"> > EJ   = {EJ  [1]}")
        clogger.debug(f"> > EE1  = {EE1 [1]}")
        clogger.debug(f"> Part 2:")
        clogger.debug(f"> > Erat = {Erat[2]}")
        clogger.debug(f"> > EJ   = {EJ  [2]}")
        clogger.debug(f"> > EE1  = {EE1 [2]}")
        clogger.debug(f"> G2(part 1) - G1(part 2):")
        clogger.debug(f"> > Erat = {G2*Erat[1] - G1*Erat[2]}")
        clogger.debug(f"> > EJ   = {G2*EJ  [1] - G1*EJ  [2]}")
        clogger.debug(f"> > EE1  = {G2*EE1 [1] - G1*EE1 [2]}")

        integral1 = +G2*(Erat[1] + EJ[1] + EE1[1])
        integral2 = -G1*(Erat[2] + EJ[2] + EE1[2])

        clogger.debug(f"> G1 = {G1}")
        clogger.debug(f"> G2 = {G2}")
        clogger.debug(f"> +G2 integral 1 = {integral1}")
        clogger.debug(f"> -G1 integral 2 = {integral2}")
        clogger.debug(f"> Integral total = {integral1+integral2}")
        clogger.debug(f"> Boundary       = {boundary}")
        clogger.debug(f"> EE1 + boundary = {G2*EE1[1] - G1*EE1[2] + boundary}")

        result = (
            + boundary
            + integral1 + integral2
            )

        if ctx.err_log:
            norm = abs(QuadError.get_error(result))**2
            print(f"""{t}{TAB}{TAB.join(str(abs(QuadError.get_error(x))**2 / norm)
                            for x in (
                                boundary,
                                G2*Erat[1] - G1*Erat[2],
                                G2*EJ[1] - G1*EJ[2],
                                0,
                                0,
                                G2*EE1[1] - G1*EE1[2],
                                0,
                                0,
                                ))}""", file=ctx.err_log)

        return result

    # The method that works everywhere (sec. 4.3)
    else:
        clogger.debug(f"Using contour integrals ({ctx.timelike_with_contour=}, {ctx.all_with_contour=})")

        # The parts of the definite integrals that are done along the real axis
        # (very expensive, but can be stored and reused)
        cached = get_star_point_integrals(ctx)

        # The parts of the definite integrals that are done as contour integrals
        # (moderately expensive)
        with workprec(int(mp.prec*(ctx.contour_lower_precision if ctx.contour_lower_precision else 1))):
            clogger.debug(f"Using {mp.prec}-bit precision contours")
            cont1_J, cont1_E = E5_contour_integral(1, ctx)
            cont2_J, cont2_E = E5_contour_integral(2, ctx) if n != 6 else (0,0)

        # The indefinite integrals
        # (very cheap)
        Ediv = {1: int_Hdiv   (   ctx),            2: 0                             }
        Erat = {1: int_Srat_gn(1, ctx),            2: int_Srat_gn(2, ctx)           }

        # Reassembled integrals
        EJ   = {1: cached['J'][1]+cont1_J,         2: cached['J'][2]+cont2_J        }
        EE1  = {1: cached['E'][1]+cont1_E+Ediv[1], 2: cached['E'][2]+cont2_E+Ediv[2]}

        clogger.debug(f"> [contour] Part 1:")
        clogger.debug(f"> > EJ   = {cont1_J}")
        clogger.debug(f"> > EE1  = {cont1_E}")
        clogger.debug(f"> [contour] Part 2:")
        clogger.debug(f"> > EJ   = {cont2_J}")
        clogger.debug(f"> > EE1  = {cont2_E}")

        clogger.debug(f"> [star] Part 1:")
        clogger.debug(f"> > EJ   = {cached['J'][1]}")
        clogger.debug(f"> > EE1  = {cached['E'][1]}")
        clogger.debug(f"> [star] Part 2:")
        clogger.debug(f"> > EJ   = {cached['J'][2]}")
        clogger.debug(f"> > EE1  = {cached['E'][2]}")

        clogger.debug(f"> Part 1:")
        clogger.debug(f"> > Erat = {Erat[1]}")
        clogger.debug(f"> > EJ   = {EJ  [1]}")
        clogger.debug(f"> > EE1  = {EE1 [1]}")
        clogger.debug(f"> > Ediv = {Ediv[1]}") # NOTE this is part of EE1 here
        clogger.debug(f"> Part 2:")
        clogger.debug(f"> > Erat = {Erat[2]}")
        clogger.debug(f"> > EJ   = {EJ  [2]}")
        clogger.debug(f"> > EE1  = {EE1 [2]}")
        clogger.debug(f"> G2(part 1) - G1(part 2):")
        clogger.debug(f"> > Erat = {G2*Erat[1] - G1*Erat[2]}")
        clogger.debug(f"> > EJ   = {G2*EJ  [1] - G1*EJ  [2]}")
        clogger.debug(f"> > EE1  = {G2*EE1 [1] - G1*EE1 [2]}")

        integral1 = +G2*(Erat[1] + EJ[1] + EE1[1])
        integral2 = -G1*(Erat[2] + EJ[2] + EE1[2])

        # clogger.debug(f"> +G2 contour 1 = {cont1}")
        # clogger.debug(f"> -G1 contour 2 = {cont2}")
        # clogger.debug(f"> Contour total = {cont1+cont2}")
        # clogger.debug(f"> +G2 realint 1 = {real1}")
        # clogger.debug(f"> -G1 realint 2 = {real2}")
        # clogger.debug(f"> Realint total = {real1+real2}")
        # clogger.debug(f"> +G2 all int 1 = {cont1+real1}")
        # clogger.debug(f"> -G1 all int 2 = {cont2+real2}")
        # clogger.debug(f"> All int total = {cont1+cont2+real1+real2}")
        # clogger.debug(f"> Boundary      = {boundary}")
        # clogger.debug(f"> Total         = {boundary+cont1+cont2+real1+real2}")

        clogger.debug(f"> G1 = {G1}")
        clogger.debug(f"> G2 = {G2}")
        clogger.debug(f"> +G2 integral 1 = {integral1}")
        clogger.debug(f"> -G1 integral 2 = {integral2}")
        clogger.debug(f"> Integral total = {integral1+integral2}")
        clogger.debug(f"> Boundary       = {boundary}")
        clogger.debug(f"> EE1 + boundary = {G2*EE1[1] - G1*EE1[2] + boundary}")

        result = (
            + boundary
            + integral1 + integral2
            )

        if ctx.err_log:
            norm = abs(QuadError.get_error(result))**2
            print(f"""{t}{TAB}{TAB.join(str(abs(QuadError.get_error(x))**2 / norm)
                            for x in (
                                boundary,
                                G2*Erat[1] - G1*Erat[2],
                                G2*EJ[1] - G1*EJ[2],
                                G2*cached['J'][1] - G1*cached['J'][2],
                                G2*cont1_J - G1*cont2_J,
                                G2*EE1[1] - G1*EE1[2],
                                G2*cached['E'][1] - G1*cached['E'][2],
                                G2*cont1_E - G1*cont2_E,
                                ))}""", file=ctx.err_log)

        return result

def Ebar_divergence(n, order, t, beta=None):

    if beta is None and n > 3:
        beta = t_to_beta(t)

    match (n,order):
        case (0,1):
            return (pi*pi + 35)/2
        case (0,2):
            return 23/3
        case (0,3):
            return 2

        case (1,1):
            return (t**2 - 54*t + 18*pi*pi + 630)/36
        case (1,2):
            return (23 - t)/3
        case (1,3):
            return 2

        case (2,1):
            return (24 - t + 2*pi*pi)/8
        case (2,2):
            return (28 - t)/12
        case (2,3):
            return 1

        case (3,1):
            return (t - 20)/8
        case (3,2):
            return -5/6

        case (4,1):
            return (t + 36 + 3*pi*pi - 51*Jbub(1,beta) + 9*Jbub(2,beta))/12
        case (4,2):
            return (14 - 9*Jbub(1,beta))/6
        case (4,3):
            return 1

        case (5,1):
            Jbub1 = Jbub(1,beta)
            return Jbub1**2 - Jbub1 + Jbub(2,t)/2 + pi*pi/12 + 1/3
        case (5,2):
            return (1 - Jbub(1,beta))/3
        case (5,3):
            return 1/3

        case (6,1):
            Jbub1 = Jbub(1,beta)
            return (4*Jbub1**2 + 10*Jbub1 + Jbub(2,beta) - 4)/(4*(t-4))
        case (6,2):
            return (2 + Jbub(1,beta))/(2*(t-4))

        case _:
            return 0

if __name__ == '__main__':
    from .clogging import init_clogging
    init_clogging(logging.DEBUG)

    import sys
    beta = float(sys.argv[1])
    L    = float(sys.argv[2])

    ctx = IntegrationContext(None,None,beta, method=Method.EISENSTEIN)

    # @tabulate_return
    def F1_integrand(xi):
        return h1(xi)*E_2d(1, ctx.but(beta=xi))

    F1_1 = QuadError.from_quad(F1_integrand, [(beta, L)])
    F1_2 = (
        - int_Hreg(ctx)
        + hybrid_integral(
            lambda xi: h1(xi) * (E_2d(1, ctx.but(beta=xi)) - E1_div(xi)), Hreg_series(),
            limit=L, xover=1000, int_real=True, var_real=True,
            **ctx.options)
        + QuadError.from_quad(lambda xi: h1(xi)*E1_div(xi), [(beta, L)])
        )

    print(f"[1] F1(beta, L) = {F1_1}")
    print(f"[2] F1(beta, L) = {F1_2}")

