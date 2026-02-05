# Calculation methods

from enum import Enum

class Method(Enum):
    DOUBLE_BESSEL   = 'B'
    BESSEL          = 'b'
    ELLIPTIC        = 'l'
    EISENSTEIN      = 'e'
    POISSON         = 'p'
    SECDEC          = 's'
    EXPANSION_0     = '0'
    EXPANSION_INF   = '∞'
    # TSIL            = 't'
    # KALMAN          = 'k'
    AMFLOW          = 'a'

    def needs_IBP_derivs(self):
        return self == self.SECDEC or self == self.AMFLOW
