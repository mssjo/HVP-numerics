
import functools
import logging
import threading
import sys
import time

from io import StringIO

from collections.abc import Sequence, Mapping, Set
from inspect import signature
from mpmath import mpc

from .clogging import clogger

epsilon = 1e-12
tolerance = epsilon

# So that these can be used in f-strings
NL = '\n'
TAB = '\t'

class CapturedOutput:
    def __init__(self, allow_stderr=False):
        self._allow_stderr = allow_stderr
        self.output = []
        self.errors = []

    def __enter__(self):
        self._stdout = sys.stdout
        self._stderr = sys.stderr
        sys.stdout = self._outfile = StringIO()
        sys.stderr = self._errfile = StringIO()

        return self

    def __exit__(self, exc_type, exc, traceback):
        self.output = self._outfile.getvalue().splitlines()
        self._outfile.close()
        self.errors = self._errfile.getvalue().splitlines()
        self._errfile.close()

        sys.stdout = self._stdout
        sys.stderr = self._stderr

        if not self._allow_stderr and self.errors:
            raise IOError("Captured stderr:\n" + '\n'.join(self.errors))

    def stdout(self):
        return self._outfile;
    def stderr(self):
        return self._errfile;

def abbreviate(obj, length=2):
    full_length = 2*length + 1
    if callable(obj):
        return "<func>"
    elif isinstance(obj, Sequence):
        return "<list>"
    elif isinstance(obj, Mapping):
        return "<dict>"
    elif isinstance(obj, Set):
        return "<set>"
    elif not isinstance(obj, str):
        obj = str(obj)
    return obj if len(obj) <= full_length else f"{obj[:length]}…{obj[-length:]}"

def log_errors(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as err:
            clogger.error(f"{type(err).__name__} in {func.__name__}({', '.join(f'{k}={v}' for k,v in signature(func).bind(*args, **kwargs).arguments.items())})")
            raise err
    return wrapper

def add_QuadError(func):
    from .integration import QuadError

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        result = func(*args, **kwargs)
        if isinstance(result, QuadError):
            return result
        else:
            return QuadError(result, 0)
    return wrapper


# Function decorator to cancel evaluation after a timeout
# If interrupted by the timeout or manually, the return value is default
def timeout(seconds, default=None):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            def interrupt():
                # print to stderr, unbuffered in Python 2.
                clogger.info(f"Function {func.__name__} interrupted after {seconds} seconds")
                raise KeyboardInterrupt

            timer = threading.Timer(seconds, interrupt)
            timer.start()
            try:
                result = func(*args, **kwargs)
            except KeyboardInterrupt:
                return default
            finally:
                timer.cancel()
            return result
        return wrapper
    return decorator

def timed(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        t0 = time.time_ns()
        result = func(*args, **kwargs)
        t1 = time.time_ns()
        return result, (t1-t0)
    return wrapper

# Function decorator to print return values
def print_return(func, level=logging.DEBUG):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        arguments = signature(func).bind(*args, **kwargs).arguments
        clogger.log(level, f"Calling function {func.__name__}({', '.join(f'{k}={v}' for k,v in arguments.items())})")
        retval = func(*args, **kwargs)
        clogger.log(level, f"        function {func.__name__}({','.join(abbreviate(v) for v in arguments.values())}) returns {retval}")
        return retval
    return wrapper
# Similar, but tabulates them in a file every time the function is called
class tabulate_return():

    def __init__(self, func, dir='./tabulate_return', level=logging.INFO, print_level=logging.DEBUG):

        self._func = func
        self._file = open(f"{dir}/{func.__name__}.dat", 'w')
        self._level = level
        clogger.log(print_level, f"Return values of {func.__name__} will be printed to {self._file.name}")
        if clogger.level <= self._level:
            print('\t'.join(signature(func).parameters.keys()) + '\tretval', file=self._file)

    def __del__(self):
        self._file.close()

    def __call__(self, *args, **kwargs):
        from .integration import QuadError # import here to avoid circular dependency
        def flatten(z):
            if isinstance(z, QuadError):
                return flatten(z.value())
            elif is_complex(z):
                return str(z.real) + '\t' + str(z.imag)
            else:
                return str(z)

        retval = self._func(*args, **kwargs)
        if clogger.level <= self._level:
            print('\t'.join(flatten(z) for z in signature(self._func).bind(*args, **kwargs).arguments.values()) + '\t' + flatten(retval), file=self._file)
        return retval

# Check if z is a complex number type, regardless of whether it represents a real number
def is_complex(z):
    return hasattr(z, 'imag')

def approx_equal(val_1, val_2, comp=None, epsabs=1e-6, epsrel=1e-4):
    if hasattr(val_1, 'approx_equal'):
        return val_1.approx_equal(val_2, comp, epsabs, epsrel)
    if hasattr(val_2, 'approx_equal'):
        return val_2.approx_equal(val_1, comp, epsabs, epsrel)
    if comp is None:
        comp = (abs(val_1) + abs(val_2))/2
    diff = abs(val_1 - val_2)

    if diff < epsabs or (abs(comp) > epsabs and abs(diff/comp) < epsrel):
        return True
    else:
        # clogger.debug(f"approx_equal failed: {diff} >= {epsabs} (abs), {f'no rel (comp={abs(comp)})' if abs(comp) <= epsabs else f'{abs(diff/comp)} >= {epsrel} (rel)'}")
        return False


# Check for real or purely imaginary within some tolerance
def is_real(z, epsabs=1e-6, epsrel=1e-4):
    if hasattr(z, 'is_real'):
        return z.is_real(epsrel, epsabs)
    return (not is_complex(z)) or approx_equal(z.imag, 0., z.real, epsabs,epsrel)
def is_imag(z, epsabs=1e-6, epsrel=1e-4):
    if hasattr(z, 'is_imag'):
        return z.is_imag(epsrel, epsabs)
    return is_complex(z) and approx_equal(z.real, 0., z.imag, epsabs,epsrel)

# Real imaginary parts, expecting the number to be almost purely real/imaginary
def Re(z, epsabs=1e-6, epsrel=1e-4, strict=False):
    if not is_complex(z):
        return z
    if not is_real(z, epsabs, epsrel):
        if strict:
            raise ValueError(f"Expecting (approximately) real value, but re={z.real}, im={z.imag}")
        else:
            clogger.warning(f"Expecting (approximately) real value, but re={z.real}, im={z.imag}")
    return z.real
def Im(z, epsabs=1e-6, epsrel=1e-4, strict=False):
    if not is_complex(z):
        return 0
    if not is_imag(z, epsabs, epsrel):
        if strict:
            raise ValueError(f"Expecting (approximately) imaginary value, but re={z.real}, im={z.imag}")
        else:
            clogger.warning(f"Expecting (approximately) imaginary value, but re={z.real}, im={z.imag}")
    return z.imag

# Take real part of result if all arguments are (close enough to) real
def real_in_real_out(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        result = func(*args, **kwargs)
        if all((not is_complex(arg)) or is_real(arg) for arg in list(args) + list(kwargs.values())):
            try:
                result = Re(result)
            except ValueError:
                clogger.warning(f"Real-in-real-out failed: {result} not real")
        return result
    return wrapper

def verify_inverse(ifun, *, epsabs=1e-6, epsrel=1e-4):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            assert len(args) == 1, "verify_inverse not implemented for multi-argument functions"
            retval = func(*args, **kwargs)
            invval = ifun(retval, **kwargs)

            if not approx_equal(invval, args[0], epsabs=epsabs, epsrel=epsrel):
                clogger.warning(f"Inverse verification failed: {ifun.__name__}({func.__name__}({', '.join(f'{k}={v}' for k,v in signature(func).bind(*args, **kwargs).arguments.items())})) = {args[0]}")
            return retval
        return wrapper
    return decorator

# Ordinal numbers
def ordinal(n):
    if n < 0:
        raise ValueError(f"Ordinal number should be possitive, got {n}")
    match n % 10:
        case 1 if n != 11:
            return f"{n}st"
        case 2 if n != 12:
            return f"{n}nd"
        case 3 if n != 13:
            return f"{n}rd"
        case _:
            return f"{n}th"

