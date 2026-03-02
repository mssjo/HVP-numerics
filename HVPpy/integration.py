from collections.abc import Sequence, Mapping
from copy import deepcopy as copy
from math import prod

import mpmath
from mpmath import quad, quadsubdiv, quadosc, nsum
from mpmath import sqrt, log, factorial, ceil
from mpmath import nstr, mpc, inf

from .clogging import clogger
from .method import Method
from .utilities import *

# The IntegrationContext holds all the variables necessary to carry out E_2d and Ebar.
# Without it, both functions become a mess of deeply nested inner functions.
# This also allows a more smooth implementation of on-demand variables:
# Only one of t,tau,beta and theta need be defined, and the others are computed from them as needed.
# Contexts can also be passed when the methods call each other instead of duplicating all the variables.
# If a variable is set to False, attempting to compute it raises an error.
class IntegrationContext:
    """
    Class for holding all variables necessary to carry out E_2d and Ebar.

    The reasons for havig this class are theefold.
    Firstly, it reduces the number of arguments needed in the auxiliary methods
    to E_2d and Ebar: now, only a context and possibly some index needs to be passed.
    Secondly, it allows flexibility in which of t, tau and beta is supplied:
    the appropriate variable is requested from the context, and if it is missing,
    it is computed on the fly from whichever variable has been provided. It is
    then stored, which is important since some variable conversions are expensive.
    Thirdly, it flexibly carries all miscellaneous options and reports missing
    ones as false, making the logic of using them cleaner
    (context.option rather than kwargs.get("option", False)).
    """

    def __init__(self, t,tau,beta, *, context=None, **options):
        """
        Create an integration context.

        Parameters:
        t, tau, beta - the kinematic variables.
            If a variable is set to None, its value is lazily determined from that of the others.
            If a variable is set to False, it is disabled for debugging, e.g. "assert t is not used"
        context - if provided, all options from this context are copied here.
        options - any number of options to be set for this context.
            Options present here override those set in context.
        """

        if all((x is None or x is False) for x in (t,tau,beta)):
            raise ValueError("At least one of t, tau or beta must be supplied")

        self._t         = t
        self._tau       = tau
        self._beta      = beta

        if context is not None:
            for opt, val in context.options.items():
                setattr(self, opt, val)

        for opt, val in options.items():
            setattr(self, opt, val)

    def but(self, *, t=None, tau=None, beta=None, **options):
        """
        Make a modified copy of this context.

        Parameters:
        t, tau and beta - new values for the kinematic variables.
            If neither is provided, the original ones are retained.
            If at least one is provided, the others are lazily determined from it.
        options - any number of options to be set.
            Original options are retained unless overridden here.
        """
        if all(x is None for x in (t,tau,beta)):
            t = self._t
            tau = self._tau
            beta = self._beta

        return IntegrationContext(
            t, tau, beta,
            context=self,
            cache = {}, # Cache is invaliated so it is cleared
            **options
            )

    def get(self, attr, default=False, type=None):
        if attr.startswith('_'):
            raise AttributeError(f"IntegrationContext has no attribute '{attr}'")
        if type is None:
            return self.__dict__.get(attr, default)
        else:
            return type(self.__dict__.get(attr, default))

    def __getattr__(self, attr):
        """
        Provide lazily determined t, tau and beta, and default missing attributes to False.
        """
        match attr:
            case 't':
                if self._t is False:
                    raise AttributeError(f"t is disabled for this context")
                if self._t is None:
                    if self._beta is not None:
                        from .t_tau_beta import beta_to_t
                        self._t = beta_to_t(self._beta)
                    else:
                        from .t_tau_beta import tau_to_t
                        self._t = tau_to_t(self.tau)
                return self._t

            case 'tau':
                if self._tau is False:
                    raise AttributeError(f"tau is disabled for this context")
                if self._tau is None:
                    if self.use_theta and (not is_complex(self.t) or self.t.imag == 0):
                        from .theta import theta_to_tau, t_to_theta
                        self._tau = QuadError.decay(theta_to_tau(t_to_theta(self.t)))
                    else:
                        from .t_tau_beta import t_to_tau
                        self._tau = t_to_tau(self.t)
                return self._tau

            case 'beta':
                if self._beta is False:
                    raise AttributeError(f"beta is disabled for this context")
                if self._beta is None:
                    from .t_tau_beta import t_to_beta
                    self._beta = t_to_beta(self.t)
                return self._beta

            case 'cache':
                self.cache = {}
                return self.cache # Cache defaults to empty

            case 'options':
                return {attr: val for attr, val in self.__dict__.items() if not attr.startswith('_')}

            case _:
                if attr.startswith('_'): # Hide private attributes
                    raise AttributeError(f"IntegrationContext has no attribute '{attr}'")
                else:
                    return False # Options default to False

    def __bool__(self):
        """ Always return True in contrast to None """
        return True

    def __str__(self):
        """ Write a summary of what is stored in the context """
        items = (
            f"t={self._t}" if self._t is not None else None,
            f"τ={self._tau}" if self._tau is not None else None,
            f"β={self._beta}" if self._beta is not None else None,
            f"{self.method.value}" if "method" in self.__dict__ else None,
            "θ" if self.use_theta else None,
            )
        return f"({','.join(item for item in items if item is not None)})"

class QuadError:
    """
    Class for error propagation in quadrature.

    QuadError objects are meant as drop-in replacements for exact numbers,
    but carry error information that is propagated through all arithmetic operations.
    For efficiency, the total error is lazily evaluated when needed.
    """

    def __init__(self, val=0, err=0):
        """
        Create a QuadError object.

        Parameters:
        val - the central value (default: zero)
        err - the standard deviation (default: zero)

        Thus, QuadError(val) is an exact value and QuadError() is zero.
        """
        if val or err:
            if isinstance(val, QuadError):
                self._values = copy(val._values)
            else:
                self._values = [(val, err)]
        else:
            self._values = []

    # Implementation of arithmetic operations with error propagation
    def __add__(self, x):
        res = QuadError(self)
        res += x
        return res
    def __radd__(self, x):
        return self + x
    def __iadd__(self, x):
        if isinstance(x, QuadError):
            self._values.extend(x._values)
        elif isinstance(x, Sequence):
            self._values.append((QuadError.decay(x[0]), QuadError.decay(x[1]) if len(x)>1 else 0.))
        else:
            self._values.append((QuadError.decay(x),0.))
        return self
    def __sub__(self, x):
        res = QuadError(self)
        res -= x
        return res
    def __rsub__(self, x):
        return -(self - x)
    def __isub__(self, x):
        if isinstance(x, QuadError):
            self._values.extend((-val,err) for val,err in x._values)
        elif isinstance(x, Sequence):
            self._values.append((-QuadError.decay(x[0]), QuadError.decay(x[1]) if len(x)>1 else 0.))
        else:
            self._values.append((-QuadError.decay(x),0.))
        return self
    def __neg__(self):
        return QuadError() - self
    def __pos__(self):
        return QuadError(self)
    def __imul__(self, x):
        self._values = [(val*QuadError.decay(x), (QuadError.decay(x))*err) for val,err in self._values]
        return self
    def __mul__(self, x):
        res = QuadError(self)
        res *= x
        return res
    def __rmul__(self, x):
        return self * x
    def __truediv__(self, x):
        return self * (1/x);
    def __rtruediv__(self, x):
        val,err = self()
        return QuadError(QuadError.decay(x) / val, (QuadError.decay(x))*err)
    def __pow__(self, x):
        val,err = self()
        exp = QuadError.decay(x)
        powm1 = val**(exp-1)
        return QuadError(val*powm1, (err*powm1*exp))
    def __abs__(self):
        val, err = self()
        return QuadError(abs(val), abs(err))
    @staticmethod
    def exp(obj):
        from mpmath import exp

        if isinstance(obj, QuadError):
            val, err = obj()
            res = exp(val)
            return QuadError(res, res*err)
        else:
            return exp(obj)

    def __call__(self):
        """ Force evaluation of total error and return the resulting value,error pair """
        if len(self._values) == 1:
            return self._values[0]
        else:
            self._values = [(
                sum(val for val,_ in self._values),
                ( # NOTE: real and imaginary parts get independent errors
                    + sqrt(sum((err.real if is_complex(err) else err)**2 for _,err in self._values))
                    + 1j*sqrt(sum((err.imag if is_complex(err) else 0)**2 for _,err in self._values))
                )
                )]
            return self._values[0]
    def value(self):
        """ Return the central value """
        return self()[0]
    def error(self):
        """ Return the standard deviation """
        return self()[1]

    def add_error(self, err):
        """ Add an extra source of error without affecting the central value """
        self._values.append((0, QuadError.decay(err)))
        return self

    def __getattr__(self, attr):
        """
        Allow self.real and self.imag to work correctly.
        Delegate other unknown atributes to self.value(), ignoring the error.
        """
        if attr in {'real', 'imag'}:
            val, err = self()
            return QuadError(getattr(val, attr), getattr(err, attr, err))
        elif attr.startswith('_'):
            raise AttributeError(f"Attempting to obtain nonexistent private attribute '{attr}'")
        else:
            return getattr(self.value(), attr)

    def is_real(self, epsrel, epsabs):
        """
        Determine whether self is approximately real.

        Parameters:
        epsrel - value is deemed real if abs(imag/real) is less than this
        epsabs - value is deemed real if abs(imag) is less than either this or the error
        """
        val, err = self()
        return is_real(val, epsrel, max(epsabs, abs(err.imag)))
    def is_imag(self, epsrel, epsabs):
        """
        Determine whether self is approximately purely imaginary.

        Parameters:
        epsrel - value is deemed imaginary if abs(real/imag) is less than this
        epsabs - value is deemed imaginary if abs(real) is less than either this or the error
        """
        val, err = self()
        return is_imag(val, epsrel, max(epsabs, abs(err.real)))

    def __str__(self):
        val, err = self()
        return f"{val} ± {nstr(err,3)}" # The error doesn't need much precision
    def __repr__(self):
        val, err = self()
        return f"QuadError({val}, {err})"

    # Methods that treat non-QuadError objects as QuadError objects with zero error.
    # Useful in contexts where it is unclear whether something will be a QuadError or not.
    @staticmethod
    def decay(obj):
        """ Strip away error, or just return obj if it has no error. """
        return QuadError.get_value(obj)
    @staticmethod
    def get_value(obj):
        """ Return the central value, or just the value itself if already exact. """
        if isinstance(obj, QuadError):
            return obj.value()
        else:
            return obj
    @staticmethod
    def get_error(obj):
        """ Return the error, or zero itself if exact. """
        if isinstance(obj, QuadError):
            return obj.error()
        else:
            return 0.
    @staticmethod
    def is_imprecise(obj, tol=.5):
        """
        Test if the error of an object is large compared to its value.

        Parameters:
        obj - any value. A non-QuadError always results in False.
        tol - the smallest error/value ratio that is considered "large"
        """
        if isinstance(obj, QuadError):
            return abs(obj.error()) > abs(obj.value())*tol
        else:
            return False
    @staticmethod
    def val_err(obj):
        """ Return a (value,error) pair, with zero error for precise values. """
        if isinstance(obj, QuadError):
            return (obj.value(), obj.error())
        else:
            return (obj, 0)

    @staticmethod
    def from_nsum(summand, limits, **kwargs):
        """
        Get a QuadError from a numerically evaluated sum

        Parameters:
        summand - the summand, an integer-to-exact-value function
        limits - the summation limits, an ordered pair
        kwargs - options to be passed on to mpmath.nsum

        Returns the sum with the numerical tolerance as its error.
        """
        return QuadError(nsum(summand, limits, **kwargs), kwargs.get('tolerance', mp_eps()))

    @staticmethod
    def from_quad(integrand, limits, *, subdiv=False, osc=False,
                  integrand_error=False, **kwargs):
        """
        Get a QuadError from a numerically evaluated, possibly multidimentional integral.

        Parameters:
        integrand - the integrand function
        limits - the integration limits, a list of ordered pairs
        subdiv - use mpmath.quadsubdiv, suitable e.g. for discontinuous integrands
        osc - use mpmath.quadosc, suitable for oscillatory integrals
        estimate_integrand_error - normally, the error of the integrand is discarded.
            This option instead enables rough propagation of it by adding the integral of the error
            as another source of error of the integral.
        kwargs - options to be passed on to the integration method (normally mpmath.quad).

        Returns the integral with its numerically estimated error.
        """

        if "method" in kwargs and isinstance(kwargs["method"], Method):
            kwargs.pop("method") # because this interferes with quad's method
        if "quad_method" in kwargs:
            kwargs["method"] = kwargs["quad_method"]

        def trim(limits):
            # quad itself is not good at reporting these errors
            if len(limits) < 2:
                raise ValueError(f"Integration limits must contain at least a lower and and upper limit, got {limits}")
            if limits[0] >= limits[-1]:
                raise ValueError(f"Lower integration limit must be below upper, got {limits}")

            # Trim midpoints that fall outside the integration range
            return [limits[0]] + [limit for limit in limits[1:-1] if limit > limits[0] and limit < limits[-1]] + [limits[-1]]

        if subdiv:
            func = quadsubdiv
        elif osc:
            func = quadosc
        else:
            func = quad

        integral = QuadError(*func(
                        lambda var: QuadError.get_value(integrand(var)),
                        *(trim(limit) for limit in limits),
                        error=True, **kwargs))

        if integrand_error:
            err = func( lambda var: QuadError.get_error(integrand(var)),
                        *(trim(limit) for limit in limits),
                        error=False, **kwargs)
            clogger.debug(f"[QuadError.from_quad] integrand error (est.): {err}")
            integral.add_error(err)

        return integral


    @staticmethod
    def from_line_contour(integrand, points, **kwargs):
        if len(points) < 2:
            raise ValueError("At least two points must be provided")


        clogger.debug(f"Contour-integrating from {' to '.join(str(p) for p in points)}")

        result = 0
        for start, end in zip(points, points[1:]):
            result += QuadError.from_quad(
                lambda z: integrand(start*(1-z) + end*z) * (end-start),
                [(0,1)], **kwargs)
        return result

    @staticmethod
    def from_arc_contour(integrand, points, **kwargs):
        if len(points) != 2:
            raise ValueError("Exactly two points must be provided")

        start = points[0] if is_complex(points[0]) else complex(points[0],0)
        end   = points[1] if is_complex(points[1]) else complex(points[1],0)

        center = (abs(start) - abs(end)) / (2*(start.real - start.real))
        radius = hypot(start.imag, start.real - center)

        assert approx_equal(radius, hypot(end.imag, end.real - center)), "Inconsistent radius"

        arg_start = arg(start-center)
        arg_end = arg(end-center)

        def contour_func(z):
            return center + radius*exp(1j*z)
        def contour_deriv(z):
            return 1j*radius*exp(1j*z)

        clogger.debug(f"Integrating over arc centered at {center:.3f} with radius {radius:.3f} from {start:.3f} (arg={arg_start:.3f}) to {end:.3f} (arg={arg_end:.3f})")

        if arg_start > arg_end:
            arg_start,arg_end = arg_end,arg_start

        return QuadError.from_quad(
            lambda z: integrand(contour_func(z)) * contour_deriv(z),
            [(arg_start, arg_end)])

    def highlight(value, *, prefix="", suffix="", reference=False):
        """
        Highlight a value according to its precision using ANSI escapes.

        Parameters:
        value - the value, may or may not be QuadError
        prefix - string prepended to the result
        suffix - string appended to the result
        reference - a value to compare against, or True to highlight this as the
            reference for other values (default: False)

        Without a reference, significant digits are printed in white
        and non-significant ones (those falling within the error band) in gray.
        With a reference, significant digits matching those of the reference
        are printed in green, and those not matching it in red.
        The reference itself has its significant digits printed in cyan.

        The prefix, suffix, exponents and errors are printed in the "main" color
        of the value:
            green if there is at least one green digit,
            red if all digits are red,
            gray if no digit is significant (the error is comparable or larger
                than the value),
            cyan for the reference (unless gray),
            white otherwise.

        Real and complex parts are printed separately, with the prefix colored
        according to the real part and the suffix according to the imaginary.
        """
        from .clogging import ColorFormatter as cf

        has_ref = (reference is not True and reference is not False)

        # Recurse on real and imaginary parts
        if is_complex(value) and QuadError.decay(value.imag) != 0:
            return f"""{
                QuadError.highlight(value.real,
                                    prefix=f'{prefix}(', suffix=')',
                                    reference=(reference.real if has_ref else reference))
                } + {
                QuadError.highlight(value.imag,
                                    prefix='(', suffix=f')i{suffix}',
                                    reference=(reference.imag if has_ref else reference))
                }"""

        if has_ref and is_complex(reference):
            ref = reference.real
        else:
            ref = reference

        # clogger.debug(f"value: {value}, reference: {reference}")

        val, err = QuadError.val_err(value)
        val = val.real
        err = err.real

        vstr = str(val)
        estr = f" ± {nstr(err, 3)}" if isinstance(value, QuadError) else ""

        if has_ref:
            rval, rerr = QuadError.val_err(ref)
            rstr = str(rval)
            err = max(abs(err), abs(rerr))
        else:
            err = abs(err)

        # Find number of dignificant digits
        if err != 0:
            # All this because converting mpf to int is hard
            sigdig = 0
            pow10 = err
            while pow10 < abs(val):
                sigdig += 1
                pow10 *= 10
            # clogger.debug(f"{sigdig} significant digits")
        else:
            sigdig = len(vstr)
            # clogger.debug(f"All digits significant")

        # Compare exponents, when applicable
        epos = vstr.find('e') if 'e' in vstr else len(vstr)
        if has_ref:
            repos = rstr.find('e') if 'e' in rstr else len(rstr)
            if vstr[epos:] == rstr[repos:]:
                # clogger.debug(f"Value has same power ({vstr[epos:]}) as reference ({rstr[repos:]})")
                rpos = -1
            else:
                # clogger.debug(f"Value has different power ({vstr[epos:]}) than   reference ({rstr[repos:]})")
                rpos = 0

        # Explore string representation of value
        spos = -1
        ndig = 0
        for i,c in enumerate(vstr[:epos]):
            if not c.isdigit():
                continue
            if ndig or c != '0':
                ndig += 1 # Increment digit count, excluding leading zeroes
            if ndig > sigdig and spos < 0:
                spos = i  # Find point where significant digits end
                break
            if has_ref and i < len(rstr) and c != rstr[i] and rpos < 0:
                rpos = i  # Find point where value diverges from reference
                if spos >= 0:
                    break
        if has_ref and rpos < 0:
            rpos = epos
        if spos < 0:
            spos = epos

        # Completely imprecise value
        if sigdig == 0:
            return f"{cf.BRIGHTCOLOR%cf.BLACK}{prefix}{vstr}{estr}{suffix}{cf.RESET}"
        # Completely wrong value
        if has_ref and rpos == 0:
            return f"{cf.BOLD}{cf.COLOR%cf.RED}{prefix}{vstr}{estr}{suffix}{cf.RESET}"

        # Normal case
        if has_ref:
            maincolor = cf.COLOR%cf.GREEN
            return f"{cf.BOLD}{maincolor}{prefix}{vstr[0:min(spos,rpos)]}{f'{cf.COLOR%cf.RED}{vstr[rpos:spos]}' if rpos < spos else ''}{cf.BRIGHTCOLOR%cf.BLACK}{vstr[spos:epos]}{maincolor}{vstr[epos:]}{estr}{suffix}{cf.RESET}"
        else:
            maincolor = (cf.BRIGHTCOLOR%cf.CYAN if ref else cf.COLOR%cf.WHITE)
            return f"{cf.BOLD}{maincolor}{prefix}{vstr[0:spos]}{cf.BRIGHTCOLOR%cf.BLACK}{vstr[spos:epos]}{maincolor}{vstr[epos:]}{estr}{suffix}{cf.RESET}"


# Semi-numerical integral of the integrand from infinity to the limit
#  - the order of limits must be reversed for quad
#  - the integral and integration variable may be either real or purely imaginary,
#    and this must be accounted for so that only real integration is needed
#  - at the xover point (always given as a real number),
#    the numerical integral of the exact integrand is replaced by
#    the exact integral of the Taylor series of the integrand around inifinity
def hybrid_integral(integrand, series, limit, *, xover,
                         int_real, var_real,
                         other_limits=[], midpoints=[],
                         **kwargs):

    if "hybrid_error_log" in kwargs:
        t0 = time.time_ns()

    if any(coeff for n,coeff in series.items() if n >= -1):
        raise ValueError("Integral of series-expanded integrand does not converge")

    # Calculate head (finite part of integration range) numerically
    if var_real:
        limits = [Re(limit), *midpoints, xover]
        actual_var = lambda z: z
        var_prefactor = 1
    else:
        limits = [Im(limit), *midpoints, xover]
        actual_var = lambda z: 1j*z
        var_prefactor = 1j

    if int_real:
        real_integrand = lambda z: Re(z, strict=True)
        int_prefactor = 1
    else:
        real_integrand = lambda z: Im(z, strict=True)
        int_prefactor = 1j

    # Adaptive crossover
    adaptive_xover = kwargs.get("adaptive_xover", False)
    if adaptive_xover:
        xover = limit
        xover_scale = 1 + kwargs.get("xover_step", .1)
        while abs(xover) < kwargs.get("max_xover", 1000):
            head_at_xover = integrand(actual_var(xover))
            tail_at_xover = evaluate_series(series, var=actual_var(xover), error=-1)
            tail_relerr = abs(QuadError.decay((head_at_xover - tail_at_xover) / head_at_xover))
            if xover != limit and tail_relerr > prev:
                xover /= xover_scale
                break
            xover *= xover_scale
            prev = tail_relerr

    clogger.debug(f"Hybrid integral ({'adaptive ' if adaptive_xover else ''}xover at {xover}):")

    limits.append(xover)

    if limits[0] < xover:
        # -1 for swapping the integration limits
        head = -1 * var_prefactor * int_prefactor * QuadError.from_quad(
            lambda z: real_integrand(integrand(actual_var(z))),
            [[limit for limit in limits if limit <= xover], *other_limits],
            **kwargs
            )
    else:
        # This is for when the entire integral is above the crossover
        head = QuadError(0)
        xover = limits[0]

    if xover != inf:
        # Calculate tail (infinite part) from series
        # (over)estimate relative error of series approximation of integral
        # using relative error of integrands at the xover
        head_at_xover = integrand(actual_var(xover))
        tail          = evaluate_series(series, var=actual_var(xover), deriv=-1, error=-1)
        tail_at_xover = evaluate_series(series, var=actual_var(xover), error=-1)

        tail_relerr = abs(QuadError.decay((head_at_xover - tail_at_xover) / head_at_xover))
        if tail_relerr > 1e-1:
            clogger.warning(f"Large tail relative error: {tail_relerr}")
        tail = QuadError(tail.value(),
                        abs(tail.value()*tail_relerr)
                        # abs(tail.value())              # TEMP
                        )
    else:
        head_at_xover = 0
        tail_at_xover = 0
        tail_relerr = 0
        tail = QuadError(0)

    clogger.debug(f" > head at xover: {head_at_xover}")
    clogger.debug(f" > tail at xover: {tail_at_xover}")
    clogger.debug(f" > >  rel. error: {tail_relerr}")
    clogger.debug(f" > head-integral: {head}")
    clogger.debug(f" > tail-integral: {tail}")

    result = head+tail

    if "hybrid_error_log" in kwargs:
        t1 = time.time_ns()
        print(kwargs.get('source', '???') + '\t' + '\t'.join((str(x.real) for x in (
            xover,
            result.value(), abs(result.error()),
            head.value(), abs(head.error()),
            tail.value(), abs(tail.error()),
            QuadError.decay(head_at_xover),
            QuadError.decay(tail_at_xover),
            abs(tail_relerr),
            t1-t0,
            ))), file=kwargs["hybrid_error_log"])

    return result

def hybrid_error_log_header():
    return '\t'.join((
        'source',
        'chi',
        'result', 'dresult',
        'head', 'dhead',
        'tail', 'dtail',
        'head_at_chi',
        'tail_at_chi',
        'tail_relerr',
        'time'))

def evaluate_series(coeffs, var,log_var=None, deriv=0,log_deriv=0, *,
                    error=+1, max_order=False):
    """
    Evaluate a power series, or its derivative, or its logarithmic derivative

    Arguments:
    coeffs - a sequence or mapping such that coeffs[n] is the coefficient of
            var**n in the series. However, if coeffs[n] is itself a sequence or
            mapping, then coeffs[n][m] is the coefficient of var**n * log_var**m.
    var - the variable of the series
    log_var - the logarithm of the variable. If not provided, it is computed if needed.
    deriv - the number of derivatives with respect to var.
            Negative numbers give integrals, ignoring integration constants.
    log_deriv - the number of derivatives (or integrals if negative) with respect to log_var.
            Mutually exclusive with deriv.
    error - if +1 (or any truthy value except -1),
            the error of the series is estimated as the value of the
            last term (most positive power), and a QuadError object is returned
            If -1, the error is estimated as the most negative power instead.
            If false, no error estimation is made and the result is returned as-is.
    max_order - if truthy, truncate the series at this order.
            This is done as if coeffs[n] did not exist for abs(n) > max_order:
            the unused terms are not used to improve any error estimates.
    """

    # NOTE: allowing both deriv and log_deriv would work,
    #  with all derivatives taken before any log-derivatives.
    # That would be confusing, though, so we disallow it
    if deriv and log_deriv:
        raise ValueError("Simultaneous derivatives and log-derivatives are not allowed")

    # Standardize the layout of the coefficients
    if not isinstance(coeffs, Mapping):
        return evaluate_series({n: coeff for n,coeff in enumerate(coeffs)}, var,log_var, deriv,log_deriv, error=error, max_order=max_order)
    # Standardize the error
    if error and error != -1:
        error = +1

    # Inner function for evaluating a single term
    def series_term(var_pow=0,log_pow=1, deriv=0,log_deriv=0):
        # Evaluate log if needed
        nonlocal log_var
        if log_var is None:
            if (log_pow or deriv <= var_pow < 0 or (log_deriv < 0 and var_pow == 0)):
                log_var = log(var)

        if deriv == 0 and log_deriv == 0:
            if log_pow:
                return var**var_pow * log_var**log_pow
            else:
                return var**var_pow
        elif deriv > 0:
            # NOTE: short-circuiting zero paths to hopefully speed up higher derivatives
            return (
                + (var_pow * series_term(var_pow-1,log_pow,   deriv=deriv-1) if var_pow else 0)
                + (log_pow * series_term(var_pow-1,log_pow-1, deriv=deriv-1) if log_pow else 0)
                )
        elif log_deriv > 0:
            return (
                + (var_pow * series_term(var_pow,log_pow,   log_deriv=log_deriv-1) if var_pow else 0)
                + (log_pow * series_term(var_pow,log_pow-1, log_deriv=log_deriv-1) if log_pow else 0)
                )
        elif deriv < 0: # integral
            if var_pow == -1:
                return log_var**(log_pow+1) / (log_pow+1)
            else:
                # clogger.debug(f"int log(x)^{log_pow}x^{var_pow} dx =")
                # for k in range(0,log_pow+1):
                #     clogger.debug(f" {'-' if k%2 else '+'} log(x)^{log_pow-k}*x^{var_pow+1} * {factorial(log_pow)/factorial(log_pow-k)} / {(var_pow+1)**(k+1)}")
                return sum(
                    series_term(var_pow+1,log_pow-k, deriv=deriv+1)
                        * (-1)**(k)
                        * factorial(log_pow)/factorial(log_pow-k)
                        / (var_pow+1)**(k+1)
                    for k in range(0,log_pow+1)
                    )
        elif log_deriv < 0: # log-integral
            # Simply change integration variable to var
            # return series_term(log_pow, var_pow-1, deriv=-1, log_deriv=log_deriv+1)
            if var_pow == 0:
                return series_term(0,log_pow+1, log_deriv=log_deriv+1) / (log_pow+1)
            return sum(
                series_term(var_pow,log_pow-k, log_deriv=log_deriv+1)
                    * (-1)**(k)
                    * factorial(log_pow)/factorial(log_pow-k)
                    / (var_pow)**(k+1)
                for k in range(0,log_pow+1)
                )

    # Evaluate the series
    series = 0
    last_n = 0
    err = 0
    for n, coeff in coeffs.items():
        # Skip terms higher than max_order
        if max_order and abs(n) > max_order:
            continue

        if isinstance(coeff, Sequence):
            term = sum(subcoeff * series_term(n,m, deriv,log_deriv) for m,subcoeff in enumerate(coeff))
        else:
            term = coeff * series_term(n,0, deriv,log_deriv)

        series += term

        if (error == -1 and n < last_n) or (error and n > last_n):
            last_n = n
            err = copy(term)

    if error:
        # Error is crudely estimated by giving the last term 100% uncertainty
        return QuadError(series, err)
    else:
        return series



    # def series(term, ncoeff):
    #     series = sum(term(n,coeff) for n,coeff in ncoeff)
    #     if error:
    #         # Error is extimated as the last term in the series
    #         return QuadError(series, term(*ncoeff[-1]))
    #     else:
    #         return series
    #
    # if deriv and log_deriv:
    #     raise ValueError("Simultaneous derivatives and log-derivatives are not allowed")
    # elif log_deriv:
    #     if inv_var:
    #         return series(
    #             lambda n,coeff: var**-n * (-n)**log_deriv * coeff,
    #             list(enumerate(coeff))[1:]
    #             )
    #     else:
    #         return series(
    #             lambda n,coeff: var**+n * (+n)**log_deriv * coeff,
    #             list(enumerate(coeffs))[1:]
    #             )
    # elif deriv > 0:
    #     if inv_var:
    #         return series(
    #             lambda n,coeff: var**(-n+deriv) * prod(-x for x in range(n, n+deriv)) * coeff,
    #             list(enumerate(coeffs))[1:]
    #             )
    #     else:
    #         return series(
    #             lambda n,coeff: var**(+n-deriv) * prod(range(1+n-deriv, 1+n)) * coeff,
    #             list(enumerate(coeffs))[deriv:]
    #             )
    # elif deriv < 0:
    #     if inv_var:
    #         clogger.debug(' + '.join(f"x^{-n-deriv}*c_{n}/{prod(-x for x in range(n+deriv, n))}" for n,_ in list(enumerate(coeffs))[1-deriv:]) + " + logs")
    #         return (
    #             + series(
    #                 lambda n,coeff: var**(-n+deriv) / prod(-x for x in range(n+deriv, n)) * coeff,
    #                 list(enumerate(coeffs))[1-deriv:]
    #                 )
    #             + series(
    #                 lambda n,coeff: log_int(var, -(n+deriv)) / prod(-x for x in range(1, n)) * coeff,
    #                 list(enumerate(coeffs))[:1-deriv]
    #                 )
    #             )
    #     else:
    #         return series(
    #             lambda n,coeff: var**(n-deriv) / prod(range(1+n, 1+n-deriv)) * coeff,
    #             list(enumerate(coeffs))
    #             )
    # else:
    #     if inv_var:
    #         return series(lambda n,coeff: var**(-n) * coeff, list(enumerate(coeffs)))
    #     else:
    #         return series(lambda n,coeff: var**(+n) * coeff, list(enumerate(coeffs)))


# def evaluate_log_series(var, coeffs,log_var=None, deriv=0, log_deiv=0, error=True):
#
#     if log_var is None:
#         log_var = log(var)
#
#     def series(term, mncoeff):
#         def inner_series(m,ncoeff):
#             return sum(term(m,n,coeff) for n,coeff in coeffs)
#         series = sum(inner_series(m,ncoeff) for m,ncoeff in mncoeff)
#         if error:
#             # Error is extimated as the last term in the series,
#             # assuming log expansion to be exact
#             return QuadError(serie, inner_series(*ncoeff[-1]))
#         else:
#             return series
#
#
#     if deriv and log_deriv:
#         raise ValueError("Simultaneous derivatives and log-derivatives are not allowed")
#     elif log_deriv > 0:
#         return series(lambda m,n,coeff: n*
#
#     elif deriv > 0:
#         return series(lambda n,coeff: var**(n-deriv) * prod(range(1+n-deriv, 1+n)) * coeff, list(enumerate(coeffs[deriv:], start=deriv)))
#     elif deriv < 0:
#         return series(lambda n,coeff: var**(n-deriv) / prod(range(1+n, 1+n-deriv)) * coeff, list(enumerate(coeffs)))
#     else:
#         return series(
#             lambda m,n,coeff: var**m * log_var**n * coeff,
#             list(enumerate(list(enumerate(c)) for c in coeffs))
#             )

if __name__ == '__main__':
    from .series_expansion import *

    t = .1
    print(f"t = {t}")
    print(f"E1 = {evaluate_series(E_2d_series(1), var=t)}")
    print(f"dE1/dt = {evaluate_series(E_2d_series(1), var=t, deriv=1)}")
    print(f"d2E1/dt2 = {evaluate_series(E_2d_series(1), var=t, deriv=2)}")
    print(f"dE1/dlogt = {evaluate_series(E_2d_series(1), var=t, log_deriv=1)}")
    print(f"d2E1/dlogt2 = {evaluate_series(E_2d_series(1), var=t, log_deriv=2)}")
    print(f"int E1 dt = {evaluate_series(E_2d_series(1), var=t, deriv=-1)}")
    print(f"int E1 dt2 = {evaluate_series(E_2d_series(1), var=t, deriv=-2)}")
    print(f"int E1 dlogt = {evaluate_series(E_2d_series(1), var=t, log_deriv=-1)}")
    print(f"int E1 dlogt2 = {evaluate_series(E_2d_series(1), var=t, log_deriv=-2)}")
