from mpmath import log, polylog, inf, pi, zeta

def Jbub(n,beta):
    if beta == inf:
        return 0
    if beta == 0:
        match n:
            case 0:
                return 1.
            case 1:
                return -2.
            case 2:
                return 8.
            case 3:
                return -48.
            case _:
                raise IndexError(f"Jbub({n},beta) not implemented")

    bp = (beta + 1)/(2*beta)
    bm = (beta - 1)/(2*beta)
    Li = polylog

    match n:
        case 0:
            return 1
        case 1:
            return beta * (log(bp) - log(bm)) - 2
        case 2:
            return (8
                    - 2*beta*(Li(2,bp)+(1/2)*log(bp)**2+2*log(bp))
                    + 2*beta*(Li(2,bm)+(1/2)*log(bm)**2+2*log(bm))
                    )
        case 3:
            return (-48
                    + 12*beta*(Li(3,bp)+ Li(2,bp)+(1/12)*log(bp)**3+(1/2)*log(bp)**2+(2-((pi**2)/(12)))*log(bp))
                    - 12*beta*(Li(3,bm)+ Li(2,bm)+(1/12)*log(bm)**3+(1/2)*log(bm)**2+(2-((pi**2)/(12)))*log(bm))
                    + 3*beta*log(bp)*log(bm)*(log(bp)-log(bm))
                    )
        case _:
            raise IndexError(f"Jbub({n},beta) not implemented")

def Tab(a,b, beta, eps=0):
    if not all(isinstance(x, int) and x >= 0 for x in (a,b)):
        raise ValueError(f"Invalid index for Tab: {a}, {b}")
    if a+b != 3:
        raise NotImplementedError(f"Not implemented: T{a}{b}")

    match (a,b, eps):
        case (3,0, 2):
            return 21 + 5*pi**2/2 + 19*pi**4/160 - (6 + pi**2/4)*zeta(3) - 3/5*zeta(5)
        case (3,0, 1):
            return 15 + 3*pi**2/2 + 19*pi**4/480 - 3*zeta(3)
        case (3,0, 0):
            return 10 + 3*pi**2/4 - zeta(3)
        case (3,0, -1):
            return 24 + pi**2/4
        case (3,0, -2):
            return 3
        case (3,0, -3):
            return 1

        case (2,1, 0):
            return 4 + pi**2/2 - zeta(3) - (3 + pi**2/4)*Jbub(1,beta) + Jbub(2,beta) - Jbub(3,beta)/6
        case (2,1, -1):
            return (12 + pi**2 - 8*Jbub(1,beta) + 2*Jbub(2,beta))/4
        case (2,1, -2):
            return 2 - Jbub(1,beta)
        case (2,1, -3):
            return 1

        case (1,2, 0):
            Jbub1 = Jbub(1,beta)
            return 1 + pi**2/4 - zeta(3) - (2+pi**2/2)*Jbub1 + Jbub1**2 + (1-Jbub1)*Jbub(2,beta) - 1/3*Jbub(3,beta)
        case (1,2, -1):
            Jbub1 = Jbub(1,beta)
            return (4 + pi**2 - 8*Jbub1 + 4*Jbub1**2 + 4*Jbub(2,beta))/4
        case (1,2, -2):
            return 1 - 2*Jbub(1,beta)
        case (1,2, -3):
            return 1

        case (0,3, 0):
            Jbub1 = Jbub(1,beta)
            return -(zeta(3) + Jbub1**3 + 3*(pi**2/4 - Jbub(2,beta))*Jbub1 + 1/2*Jbub(3,beta))
        case (0,3, -1):
            return (pi**2 + 12*Jbub(1,beta)**2 + 6*Jbub(2,beta))/4
        case (0,3, -2):
            return -3*Jbub(1,beta)
        case (0,3, -3):
            return 1

        case _ if eps < 0:
            return 0
        case _:
            raise NotImplementedError(f"Not implemented: eps^{eps} term of T{a}{b}")

