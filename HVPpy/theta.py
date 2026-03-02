# The parametrization of the real t line in terms of the tau-friendly variable theta

from mpmath import mpf, mpc, mp
from mpmath import findroot
from mpmath import sqrt, log, exp, pi, inf
from scipy.optimize import newton

from .clogging import clogger
from .utilities import *
from .elliptics import eisen, tau_to_t, psi_t
from .integration import QuadError, evaluate_series
from .method import Method
from .t_tau_beta import beta_to_t, beta_to_tau, tau_to_beta, dt_dtau, tau_to_t

sqrt3 = sqrt(mpf(3))

def theta_to_tau(theta, deriv=0):
    # Propagate uncertainty
    if isinstance(theta, QuadError):
        val, err = theta()
        return QuadError(theta_to_tau(val, deriv), err*theta_to_tau(val, deriv+1))
    if not isinstance(theta, mpf):
        theta = mpf(theta)

    if deriv == 0:
        if theta <= 0:
            return -1j*theta/16
        elif theta <= 4:
            return (1 + exp(1j*pi*(1-theta/6)))/mpf(6)
        elif theta <= 16:
            return 1/2 + sqrt3/6 * exp(1j*pi*(34-theta)/36)
        else:
            return 1/2 + 1j*theta*sqrt3/96
    else:
        if any(theta == x for x in {0,4,16}):
            return float("nan")
        elif theta < 0:
            return -1j/16 if deriv == 1 else 0
        elif theta < 4:
            return 1/mpf(6) * exp(1j*pi*(1-theta/6)) * (-1j*pi/6)**deriv
        elif theta < 16:
            return sqrt3/6 * exp(1j*pi*(34-theta)/36) * (-1j*pi/36)**deriv
        else:
            return 1j*sqrt3/96 if deriv == 1 else 0

def theta_to_t(theta, method=None, error=True):
    if any(theta == x for x in {0,4,16}):
        return theta
    return Re(tau_to_t(theta_to_tau(theta), method=method, error=error))
def theta_to_beta(theta, method=None, error=True):
    return tau_to_beta(theta_to_tau(theta), method=method, error=error)

def dt_dtheta(theta, method=None):
    tau = theta_to_tau(theta)
    q = exp(2j*pi*tau)
    return Re(dt_dtau(tau_to_t(tau), tau, method)
            * theta_to_tau(theta, deriv=1))
def d2t_dtheta2(theta):
    q = exp(2j*pi*theta_to_tau(theta))
    return Re(
        -2j*pi*deriv(theta)
            * (1 + 6*eisen(0,2, q, psi_t))
            * theta_to_tau(theta, deriv=1)
        -2j*pi*theta_to_t(theta)
            * (1 + 6*eisen(0,2, q, psi_t))
            * theta_to_tau(theta, deriv=2)
        +24*pi**2*theta_to_t(theta)
            * (
                2*eisen(-1,3, q, psi_t)
                - eisen(-1,2, q, psi_t)
                )
            * theta_to_tau(theta, deriv=1)**2
            )

def dbeta_dtheta(beta, theta):
    q = exp(2j*pi*theta_to_tau(theta))
    return (pi*1j*(beta**2 - 1)/beta
            * (1 + 6*eisen(0,2, q, psi_t))
            * theta_to_tau(theta, deriv=1))

def t_to_theta(t, tol=None, error=True):
    t = Re(QuadError.decay(t))

    # By design, t=theta at these points
    if any(t == x for x in {0,4,16}):
        theta = QuadError(t,0) if error else t
    else:
        # Find theta using Newton's method
        # (We have the second derivative, but it's too expensive to make Halley worth it)
        # For very large |t|, theta(t) has nice asymptotics that can be used
        #  as initial guesses
        # For very small |t|, theta(t) is close to a step function and is best
        #  approximated as such (|t| < 3 is chosen rather arbitrarily)
        # For intermediate |t|, we wing it
        if t < -10:
            initial = -8/pi * log(-t)
        elif t > +20:
            initial = +48/(sqrt3*pi) * log(+t)
        elif t < 0:
            initial = min(t, -3)
        else:
            initial = max(t, +3)

        # mpmath doesn't provide nice metadata like scipy,
        # so we capture its verbose output and use that as a workaround
        captured_output = CapturedOutput()
        with captured_output:
            theta = findroot(solver='newton',
                f = lambda z: Re(QuadError.decay(theta_to_t(z) - t)),
                df= lambda z: Re(QuadError.decay(dt_dtheta(z))),
                x0=initial, maxstep=mp.prec,
                tol=tol if tol else mp_eps(), verify=error, verbose=error)

        if error:
            theta = QuadError(theta, float(captured_output.output[-1][len('error: '):]))
            clogger.debug(f"theta({t}) = {theta} found in {len(captured_output.output)//2} iterations using Newton's method")
        else:
            clogger.debug(f"theta({t}) = {theta}")

    if error and theta.value() and abs(theta.error() / theta.value()) > 1e-3:
        clogger.warning(f"Lagre relative error in theta: {theta}")

    return theta

def beta_to_theta(beta, tol=None, verify=True):
    return t_to_theta(beta_to_t(beta), tol, verify)

# Theta version of hybrid_integral
@print_return
def theta_integral(integrand, series, limit_beta, xover_beta,
                   method=None, **kwargs):

    if "hybrid_error_log" in kwargs:
        t0 = time.time_ns()

    limit_theta, limit_err = QuadError.val_err(beta_to_theta(limit_beta))

    # Adaptive crossover
    adaptive_xover = kwargs.get("adaptive_xover", False)
    if adaptive_xover:
        xover_theta = limit_theta
        xover_scale = 1/(1 + kwargs.get("xover_step", .1))
        while True:
            xover_tau = theta_to_tau(xover_theta)
            xover_beta = tau_to_beta(xover_tau)
            if abs(xover_beta) > kwargs.get("max_xover", inf):
                break
            head_at_xover = integrand(tau_to_beta(xover_tau), xover_tau)
            tail_at_xover = evaluate_series(series, var=xover_beta, error= -1)
            tail_relerr = abs(QuadError.decay((head_at_xover - tail_at_xover) / head_at_xover))
            if xover_theta != limit_theta and tail_relerr > prev:
                xover_theta /= xover_scale
                xover_tau = theta_to_tau(xover_theta)
                xover_beta = tau_to_beta(xover_tau)
                break
            xover_theta *= xover_scale
            prev = tail_relerr
        xover_err = 0
    else:
        xover_theta, xover_err = QuadError.val_err(beta_to_theta(xover_beta))
        xover_tau = theta_to_tau(xover_theta)

    head_at_xover = integrand(xover_beta, xover_tau)
    real_integrand = is_real(head_at_xover * dbeta_dtheta(xover_beta, xover_theta))

    def full_integrand(theta):
        tau = theta_to_tau(theta)
        beta = QuadError.decay(tau_to_beta(tau, method=method))
        result = QuadError.decay(integrand(beta, tau) * dbeta_dtheta(beta, theta))
        return Re(result, strict=True) if real_integrand else Im(result, strict=True)

    head = None
    if limit_theta < xover_theta < 0:
        limits = [(limit_theta, xover_theta)]
        sign = -1
    elif 0 < xover_theta < limit_theta < 4:
        limits = [(xover_theta, limit_theta)]
        sign = +1
    elif xover_theta == limit_theta:
        head = QuadError(0)
    else:
        raise ValueError(f"Invalid theta range: [{limit_theta}, {xover_theta}]")

    clogger.debug(f"Hybrid theta integral ({'adaptive ' if adaptive_xover else ''}xover at θ={xover_theta}, β={xover_beta}):")

    if not real_integrand:
        sign *= 1j
    if head is None:
        head = sign*QuadError.from_quad(full_integrand, limits, **kwargs)

    if xover_beta != inf:
        tail          = evaluate_series(series, var=xover_beta, deriv=-1, error=-1)
        tail_at_xover = evaluate_series(series, var=xover_beta, deriv= 0, error=-1)
        tail_relerr = abs(QuadError.decay((head_at_xover - tail_at_xover) / head_at_xover))
    else:
        tail = QuadError(0)
        tail_at_xover = 0
        tail_relerr = 0

    clogger.debug(f" > limit theta:   {limit_theta}")
    clogger.debug(f" > xover theta:   {xover_theta}")
    clogger.debug(f" > head at xover: {head_at_xover}")
    limit_err *= full_integrand(limit_theta)
    clogger.debug(f" > > limit error: {limit_err}")
    xover_err *= full_integrand(xover_theta)
    clogger.debug(f" > > xover error: {xover_err}")
    clogger.debug(f" > tail at xover: {tail_at_xover}")
    clogger.debug(f" > >  rel. error: {tail_relerr}")
    if tail_relerr > 1e-1:
        clogger.warning(f"Large tail relative error: {tail_relerr}")

    clogger.debug(f" > head-integral: {head}")
    head.add_error(limit_err)
    head.add_error(xover_err)
    clogger.debug(f" > > (+ lim.err): {head}")
    clogger.debug(f" > tail-integral: {tail}")
    tail.add_error(tail.value()*tail_relerr)
    clogger.debug(f" > > (+ rel.err): {tail}")

    result = head+tail

    if "hybrid_error_log" in kwargs:
        t1 = time.time_ns()
        print(kwargs.get('source', '???') + '\t' + '\t'.join((str(x.real) for x in (
            xover_beta,
            result.value(), abs(result.error()),
            head.value(), abs(head.error()),
            tail.value(), abs(tail.error()),
            QuadError.decay(head_at_xover),
            QuadError.decay(tail_at_xover),
            abs(tail_relerr),
            t1-t0
            ))), file=kwargs["hybrid_error_log"])


    return result


