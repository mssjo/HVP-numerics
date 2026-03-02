
import contextlib
import os
import subprocess

from mpmath import mp, mpf

from .psd_wrapper import encode, decode, nu_master, get_n_loop
from .utilities import CapturedOutput
from .clogging import clogger
from .integration import QuadError

class AMFlow:

    def __init__(self, n, *, dim="4", force=False, **kwargs):
        self.nu = nu_master[n]
        clogger.debug(f"Processing {self.nu=}{' (forced)' if force else ''}")

        # Due to the need to embed t in integral.wl we must defer a lot of things
        #  to __call__, but we keep the psd_wrapper semantics
        self.dim = dim
        self.kwargs = kwargs

        # Setup the key used to refer to the integral
        n_loop = get_n_loop(self.nu)
        self.key = encode(self.nu, dim, 0, name='I')
        clogger.debug(f"({n_loop=}, {self.key=})")

        # Clear existing directory if forced
        if force:
            if subprocess.run(["rm", "-rf", f"./{self.path()}"]).returncode:
                raise IOError(f"[AMFlow] Failed to remove directory for integral {self.key}")
        if not os.path.exists(self.path()):
            clogger.info(f"[AMFlow] Creating directory {self.path()}")
            os.makedirs(self.path(), exist_ok=True)

    def __enter__(self):
        clogger.debug(f"[AMFlow] Preparing run of {self.key}...")
        return self

    def __exit__(self, exc_type, exc, traceback):
        self.integral = None

    def __call__(self, t, *, prec=None):

        if prec is None:
            prec = mp.dps

        with open(f"{self.path()}integral.wl", 'w') as math:
            print('\n'.join((
                f'SetDirectory["{os.getcwd()}/{self.path()}"];',
                f'Get["../AMFlow_defs.m"];',
                f'SetAMFOptions["D0" -> {self.dim}];',
                f'integral = j[HVP{get_n_loop(self.nu)}loop, {','.join((str(n) for n in self.nu))}];',
                f'Compute[integral, {t}, {prec}];',
                f'Quit[];')),
                file=math);

            clogger.debug(f"Generated {math.name}")


        process = subprocess.run(
                ["wolframscript", "-file", "integral.wl"],
                cwd=self.path(),
                # capture_output=True, text=True
            )
        clogger.debug(f"[AMFlow] Completed run of {self.key}.")

        if process.returncode or not os.path.isfile(f"{self.path()}integral.dat"):
            raise IOError(f"[AMFlow] Failed to compute {self.key} at t={t}")

        # for line in process.stderr.splitlines():
        #     clogger.warning(f"[AMFlow] {line}")
        # for line in process.stdout.splitlines():
        #     clogger.debug(f"[AMFlow] {line}")

        with open(f"{self.path()}integral.dat", 'r') as raw:
            for line in raw:
                result = {}
                sign = +1
                for token in line.split():
                    if token == '+':
                        sign = +1
                        continue
                    elif token == '-':
                        sign = -1
                        continue

                    if (mark := token.find("*eps^")) != -1:
                        eps = +int(token[mark+len("*eps^"):])
                    elif (mark := token.find("*eps")) != -1:
                        eps = +1
                    elif (mark := token.find("/eps^")) != -1:
                        eps = -int(token[mark+len("/eps^"):])
                    elif (mark := token.find("/eps")) != -1:
                        eps = -1
                    else:
                        eps = 0

                    token = token[:mark]

                    # NOTE: opposite sign convention
                    if (mark := token.find('`')) != -1:
                        result[eps] = -QuadError(mpf(token[:mark]), mpf(f'1e-{token[mark+1:]}'))
                    else:
                        result[eps] = -QuadError(mpf(token))

        clogger.debug(f"[AMFlow] {self.key}:")
        clogger.debug(f"    {result}")

        return result;

    def path(self):
        return f"./amflow/{self.key}/"
