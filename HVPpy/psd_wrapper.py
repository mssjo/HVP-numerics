import contextlib
import os
import logging
import subprocess
import sympy as sp
import pySecDec as psd

from mpmath import mpc

from .utilities import CapturedOutput, print_return
from .clogging import clogger
from .integration import QuadError

# The propagator basis
props_3loop = [
    "l1**2 - mp2",
    "l2**2 - mp2",
    "l3**2 - mp2",
    "(l1-q)**2 - mp2",
    "(l2-q)**2 - mp2",
    "(l3-q)**2 - mp2",
    "(l1+l3)**2 - mp2",
    "(l2+l3)**2 - mp2",
    "(l1+l2)**2 - mp2"
    ]
props_nloop = {
    1 : [props_3loop[i-1] for i in (1,4)],
    2 : [props_3loop[i-1] for i in (1,2,4,5,9)],
    3 : props_3loop
    }
def get_n_loop(nu):
    return {len(props) : nloop for nloop,props in props_nloop.items()}[len(nu)]

# nu arrays for masters E0,...,E6
nu_master = [
    [1,1,0, 0,0,0, 1,1,0],
    [1,0,0, 0,1,0, 1,1,0],
    [2,0,0, 0,1,0, 1,1,0],
    [3,0,0, 0,1,0, 1,1,0],
    [1,1,0, 1,0,0, 1,1,0],
    [1,1,0, 1,1,0, 1,1,0],
    [2,1,0, 1,1,0, 1,1,0]
    ]

# Convert (nu,d) to and from a compact string representation
def encode(nu, dim='4-2*eps', maxeps=None, name='I'):
    return f"{name}_{dim.split('-')[0]}d_{f'{maxeps}eps_' if maxeps is not None else ''}{''.join(str(i) if i>0 else '0' for i in nu)}_{''.join(str(-i) if i<0 else '0' for i in nu)}"
def decode(key):
    den, num = key.rsplit('_', maxsplit=1)
    nu = [int(pos) if pos != '0' else -int(neg) for pos,neg in zip(den[-len(num):], num)]
    print(f"Decode {key} -> {nu}")
    return nu

# Wrapper around pySecDec.
# The constructor creates a integral library if it does not already exist
#   (warning: time-consuming!). Creating anew can be forced.
# It works as a context manager, and the integral is available for evaluation
#  (the integral library is loaded) inside the context.
# Evaluating the integral is done using function call syntax.
class pySecDec:

    def __init__(self, n, *, dim="4-2*eps", maxeps=0, force=False, TSIL=False, **kwargs):
        nu = nu_master[n]
        clogger.debug(f"Processing {nu=}{' (forced)' if force else ''}")

        # Setup the key used to refer to the integral
        n_loop = get_n_loop(nu)
        self.key = encode(nu, dim, maxeps, name='TSIL' if TSIL else 'I')
        clogger.debug(f"({n_loop=}, {self.key=})")

        # Clear existing directory if forced or if a previous calculation failed
        # if not os.path.exists(self.path()):
            # clogger.info(f"[pySecDec] Creating directory {self.path()}")
            # os.makedirs(self.path(), exist_ok=True)
        if not force and not os.path.isfile(f"./{self.path()}/disteval.done"):
        # if not force and not os.path.isfile(f"./{self.key}_pylink.so"):
            clogger.info(f"[pySecDec] Integral exists but incomplete - deleting and starting anew")
            force = True
        if force:
            if subprocess.run(["rm", "-rf", f"{self.path()}"]).returncode:
                raise IOError(f"[pySecDec] Failed to remove directory for integral {self.key}")

        # Generate the integal library
        try:
            massdim = sp.sympify(f'{sum(nu)} - {n_loop}*({dim})/2')
            clogger.debug(f"[pySecDec] Massdim = {massdim}")
            if TSIL:
                prefactor = f"{'-' if sum(nu)%2 else '+'}(4*pi)**({n_loop}*eps) * mp2**({massdim})"
            else:
                prefactor = f"exp({n_loop}*eps*EulerGamma)*mp2**({massdim})"

            clogger.debug(f"[pySecDec] Prefactor = {prefactor}")

            with contextlib.chdir("./secdec/"):
                integral = psd.loop_package(
                    name = self.key,
                    loop_integral = psd.LoopIntegralFromPropagators(
                        propagators = props_nloop[n_loop],
                        powerlist = nu,
                        loop_momenta = [f'l{i+1}' for i in range(n_loop)],
                        external_momenta = ['q'],
                        replacement_rules = [ ('q*q', 't * mp2') ],
                        dimensionality = dim
                        ),
                    real_parameters = ['mp2'],
                    complex_parameters = ['t'],
                    additional_prefactor = prefactor,
                    requested_orders = [0 if maxeps is None else maxeps],
                    form_work_space = '8G',
                    form_threads = 4,
                    form_optimization_level = 1
                    )

            process = subprocess.run(
                ["make", "-j", "4", "disteval"],
                cwd=self.path(),
                capture_output=True, text=True
                )
            if process.returncode:
                raise IOError(f"[pySecDec] Failed to make integral library for integral {self.key}")

            for line in process.stderr.splitlines():
                clogger.warning(f"[pySecDec] {line}")
            for line in process.stdout.splitlines():
                clogger.debug(f"[pySecDec] {line}")

        except FileExistsError as err:
            if force:
                raise err
            clogger.debug(f"[pySecDec] (Integral already exists - skipping)")


    def __enter__(self):
        clogger.debug(f"[pySecDec] Loading disteval for {self.key}...")
        self.integral = psd.integral_interface.DistevalLibrary(
                            f'{self.path()}/disteval/{self.key}.json',
                            verbose = (clogger.level < logging.DEBUG)
                        )
        # self.integral = psd.integral_interface.IntegralLibrary(
        #                     f'{self.key}/{self.key}_pylink.so'
        #                 )
        # self.integral.use_Qmc()

        return self

    def __exit__(self, exc_type, exc, traceback):
        self.integral.close()
        self.integral = None

    def __call__(self, t, *, epsrel=1e-3, epsabs=1e-10):
        if self.integral is None:
            raise RuntimeError("Disteval library not open!")
        captured_output = CapturedOutput(allow_stderr=True)
        with captured_output:
            result = self.integral(
                format = 'json',
                parameters = {'t': t, 'mp2': 1},
                epsrel = epsrel, epsabs = epsabs)
        for line in captured_output.errors:
            clogger.debug(f"[pySecDec] {line}") # psd sends debug to stderr
        for line in captured_output.output:
            clogger.debug(f"[pySecDec] {line}")

        clogger.debug(f"[pySecDec] {self.key}:")
        clogger.debug(f"    {result}")

        # The return format is a bit awkward.
        # This changes it into an {order: (value, error)} dict
        #  with (value, error) as a QuadError object
        # NOTE: overall sign due to convention difference
        return {order: QuadError(-mpc(val),mpc(err)) for (order,),(val,err) in result["sums"][self.key].items()}

    def path(self):
        return f"./secdec/{self.key}/"
