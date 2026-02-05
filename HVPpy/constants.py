import numpy as np
from .integration import QuadError

pi16 = 1/(16*np.pi*np.pi)
pi0_PDG25_MeV = QuadError(134.9768, 0.0005)
piC_PDF25_MeV = QuadError(139.57039,  0.00017)
mu_MeV = 770
fpi_MeV = 92.2

def xi_f(mpi, fpi):
    return pi16 * mpi * mpi / (fpi * fpi)

def Lpi_f(mpi, mu):
    return 2*np.log(mpi/mu)

Lpi_c = Lpi_f(pi0_PDG25_MeV.value(), mu_MeV)
xi_c = xi_f(pi0_PDG25_MeV.value(), fpi_MeV)

# Literature (lit) values
lq1_lit = QuadError(-0.65, 0.10)
lq2_lit = QuadError(0.273, 0.033)
lq4_lit = QuadError(0.92,0.20)
lq6_lit = QuadError(-1.96,0.07)
lq_pilComb_216 = lq2_lit -2*lq1_lit - 2*lq6_lit
lq_pilComb_126 = 2*lq1_lit - lq2_lit - 2*lq6_lit
rqV1_lit = QuadError(-5.0, 8.0)
rqV2_lit = QuadError(4.0,1.2)

class LiteratureLEC:
    def __init__(
        self,
    ):
        self.lq1 = lq1_lit
        self.lq2 = lq2_lit
        self.lq4 = lq4_lit
        self.lq6 = lq6_lit
        self.rqV1 = rqV1_lit
        self.rqV2 = rqV2_lit
