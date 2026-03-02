from mpmath import sqrt, besseli, besselj, besselk, inf

from .utilities import *

def bessel_real(v,t,x):
    # A variant of the modified Bessel function appearing in E1
    # that is manifestly real for real arguments
    rt = t.real
    if rt < 0:
        arg = x*sqrt(-rt)
        return besselj(v, arg) * (-arg)**v
    else:
        arg = x*sqrt(+rt)
        return besseli(v, arg) * (+arg)**v

# @tabulate_return
def bessel_integrand(n,t,x, d_logt=0):

    if x == 0:
        return 0

    match (n,d_logt):
        case (1,0):
            # k0 goes to zero faster than i0 blows up, but we have to insert the limit manually
            i0 = bessel_real(0, t,x)
            if i0 == inf:
                return 0
            else:
                return -8*x * i0 * besselk(0,x)**4

        case (1,1):
            i1 = bessel_real(1, t,x)
            if i1 == inf:
                return 0
            else:
                return -4*x * i1 * besselk(0,x)**4

        case (1,2):
            i0 = bessel_real(0, t,x) * t * x**2
            i1 = bessel_real(1, t,x) * 2
            i2 = bessel_real(2, t,x)
            if i0 == inf or i1 == inf or i2 == inf:
                return 0
            else:
                return -x * (i0 + i1 + i2) * besselk(0,x)**4

        # TODO 3rd derivative

        case (1,_):
            raise ValueError(f"{ordinal(d_logt)} log(t)-derivative not implemented")
        case (_,d) if d != 0:
            raise ValueError("log(t)-derivatives implemented only for n=1")

        case (2,0):
            i0 = bessel_real(0, t,x)
            i1 = bessel_real(1, t,x)/2
            if i0 == inf or i1 == inf:
                return 0
            else:
                return 2*x * (i0 + i1) * besselk(0,x)**4

        case (3,0):
            i0 = bessel_real(0, t,x)*(2+x**2)
            i1 = bessel_real(1, t,x)*(1+2/t)
            i2 = bessel_real(2, t,x)/t
            if i0 == inf or i1 == inf or i2 == inf:
                return 0
            else:
                return -x/2 * (i0 + i1 + i2) * besselk(0,x)**4

        case _:
            raise ValueError(f"Not implemented: E{n} in 2d (Bessel)")
