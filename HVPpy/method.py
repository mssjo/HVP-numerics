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
    FEYNTROP        = 'f'
    GRID            = 'g'
    AUTO            = '*'
    DUMMY           = ''

    def needs_IBP_derivs(self):
        return self.value in "saf"
    def needs_i_epsilon(self):
        return self.value not in "sf"
    def is_concrete(self):
        return self.value not in "*"
    def has_epsilon(self):
        return self.value in "sfa"

