
import contextlib
import os
import subprocess

from mpmath import mp, mpf

from .psd_wrapper import encode, decode, nu_master, get_n_loop
from .utilities import CapturedOutput
from .clogging import clogger
from .integration import QuadError

class AMFlow:

    def __init__(self, nu, *, dim="4", force=False, verbose=False, **kwargs):
        self.verbose = verbose
        if self.verbose:
            print(f"Processing {nu=}{' (forced)' if force else ''}", end=' ')

        # Setup the key used to refer to the integral
        n_loop = get_n_loop(nu)
        self.key = encode(nu, dim, 0, name='I')
        if verbose:
            print(f"({n_loop=}, {self.key=})")

        # Clear existing directory if forced
        if not os.path.exists(self.path()):
            clogger.info(f"[AMFlow] Creating directory {self.path()}")
            os.makedirs(self.path(), exist_ok=True)
        if force:
            if subprocess.run(["rm", "-rf", f"./{self.path()}"]).returncode:
                raise IOError(f"[AMFlow] Failed to remove directory for integral {self.key}")

        with open(f"{self.path()}integral.wl", 'w') as math:
            print(f'SetDirectory["{os.getcwd()}/{self.path()}"];', file=math);
            print(f'Get["../AMFlow_defs.m"];', file=math);
            print(f'SetAMFOptions["D0" -> {dim}];', file=math);
            print(f'integral = j[HVP{n_loop}loop, {','.join((str(n) for n in nu))}];', file=math);
            print(f'Compute[integral];', file=math);

            clogger.debug(f"Generated {math.name}")


    def __enter__(self):
        if self.verbose:
            clogger.info(f"[AMFlow] Preparing run of {self.key}...")
        return self

    def __exit__(self, exc_type, exc, traceback):
        self.integral = None

    def __call__(self, t, *, prec=mp.prec):

        # This is a workaround since wolframscript drops command-line arguments
        with contextlib.chdir(self.path()):
            with open("t.m", 'w') as out:
                print(f"{t}", file=out)
            with open("prec.m", 'w') as out:
                print(f"{prec}", file=out)

        process = subprocess.run(
                ["wolframscript", "-file", "integral.wl"],
                cwd=self.path(),
                capture_output=True, text=True
            )

        if process.returncode or not os.path.isfile(f"{self.path()}integral.dat"):
            raise IOError(f"[AMFlow] Failed to compute {self.key} at t={t}")

        for line in process.stderr.splitlines():
            clogger.warning(f"[AMFlow] {line}")
        for line in process.stdout.splitlines():
            if self.verbose:
                clogger.info(f"[AMFlow] {line}")
            else:
                clogger.debug(f"[AMFlow] {line}")

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

                    if (mark := token.find('`')) != -1:
                        result[eps] = QuadError(mpf(token[:mark]), mpf(f'1e-{token[mark+1:]}'))
                    else:
                        result[eps] = QuadError(mpf(token))

        if self.verbose:
            clogger.debug(f"[AMFlow] {self.key}:")
            clogger.debug(f"    {result}")

        return result;

    def path(self):
        return f"./amflow/{self.key}/"
