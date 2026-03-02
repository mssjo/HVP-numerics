
import contextlib
import json
import os
import subprocess

from mpmath import mp, mpf, mpc

from .psd_wrapper import encode, decode, nu_master
from .utilities import Re, Im
from .clogging import clogger
from .integration import QuadError

graph_master = [
    [((0,1), 1), ((0,1), 1), ((0,1), 1), ((0,1), 1)],
    [((0,1), 1), ((0,1), 1), ((0,1), 1), ((0,1), 1)],
    [((0,1), 2), ((0,1), 1), ((0,1), 1), ((0,1), 1)],
    [((0,1), 3), ((0,1), 1), ((0,1), 1), ((0,1), 1)],
    [((0,1), 1), ((0,2), 1), ((0,2), 1), ((0,2), 1), ((2,1), 1)],
    [((0,2), 1), ((0,3), 1), ((2,3), 1), ((2,3), 1), ((2,1), 1), ((3,1), 1)],
    [((0,2), 2), ((0,3), 1), ((2,3), 1), ((2,3), 1), ((2,1), 1), ((3,1), 1)],
    ]
def scalar_prods_master(n, t):
    pt = float(str(+Re(t, strict=True)))
    mt = float(str(-Re(t, strict=True)))
    # pt = str(+t))
    # mt = str(-t))

    if n == 0:
        return [[0, 0], [0, 0]]

    match n:
        case 1 | 2 | 3:
            int_vert = 0
        case 4:
            int_vert = 1
        case 5 | 6:
            int_vert = 2

    return [
        [pt, mt] + [0]*int_vert,
        [mt, pt] + [0]*int_vert
        ] + [[0]*(int_vert+2)] * int_vert


class FeynTrop:

    def __init__(self, n, *, dim=4, maxeps=0, force=False, **kwargs):
        nu = nu_master[n]
        clogger.debug(f"Processing {nu=}{' (forced)' if force else ''}")

        # Due to the need to embed t in params.json we must defer a lot of things
        #  to __call__, but we keep the psd_wrapper semantics
        self.n = n
        self.dim = int(dim)
        self.maxeps = maxeps
        self.kwargs = kwargs

        # Setup the key used to refer to the integral
        self.key = encode(nu_master[n], str(dim), maxeps, name='I')
        clogger.debug(f"({self.key=})")

        # Clear existing directory if forced
        if force:
            if subprocess.run(["rm", "-rf", f"./{self.path()}"]).returncode:
                raise IOError(f"[FeynTrop] Failed to remove directory for integral {self.key}")
        if not os.path.exists(self.path()):
            clogger.info(f"[FeynTrop] Creating directory {self.path()}")
            os.makedirs(self.path(), exist_ok=True)

    def __enter__(self):
        clogger.debug(f"[FeynTrop] Preparing run of {self.key}...")
        return self

    def __exit__(self, exc_type, exc, traceback):
        pass

    def __call__(self, t):

        params = {
            "graph": graph_master[self.n],
            "dimension": self.dim,
            "scalarproducts": scalar_prods_master(self.n, t),
            "masses_sqr": [1] * len(graph_master[self.n]),
            "num_eps_terms": self.maxeps+1,
            "lambda": self.kwargs.get("lambda", 0),
            "N": self.kwargs.get("samples", 1_000_000),
            "seed": 0
            }

        json_params = json.dumps(params)
        clogger.debug(f"[FeynTrop] {json_params}")

        process = subprocess.Popen(["feyntrop"],
                             cwd=self.path(),
                             stdout=subprocess.PIPE, stdin=subprocess.PIPE,
                             encoding='utf8')
        out, err = process.communicate(json_params)

        if process.returncode:
            clogger.debug(f"[FeynTrop] {err}")
            raise IOError(f"[FeynTrop] Failed to compute {self.key} at t={t}")
        clogger.debug(f"[FeynTrop] Completed run of {self.key}.")

        result = {
            # NOTE: sign convention
            order : -(-1)**sum(nu for nu in nu_master[self.n]) * QuadError(mpc(re_val, im_val), mpc(re_err, im_err))
                for order, ((re_val, re_err), (im_val, im_err))
                in enumerate(json.loads(out)["integral"])
            }


        clogger.debug(f"[FeynTrop] {self.key}:")
        clogger.debug(f"    {result}")

        return result

    def path(self):
        return f"./feyntrop/{self.key}/"
