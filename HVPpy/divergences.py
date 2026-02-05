def Ebar_divergence(n, order, t, beta=None):

    if beta is None and n > 3:
        beta = t_to_beta(t)

    match (n,order):
        case (0,1):
            return (pi**2 + 35)/2
        case (0,2):
            return 23/3
        case (0,3):
            return 2

        case (1,1):
            return (t**2 - 54*t + 18*pi**2 + 630)/36
        case (1,2):
            return (23 - t)/3
        case (1,3):
            return 2

        case (2,1):
            return (24 - t + 2*pi**2)/8
        case (2,2):
            return (28 - t)/12
        case (2,3):
            return 1

        case (3,1):
            return (t - 20)/8
        case (3,2):
            return -5/6

        case (4,1):
            return (t + 36 + 3*pi**2 - 51*Jbub(1,beta) + 9*Jbub(2,beta))/12
        case (4,2):
            return (14 - 9*Jbub(1,beta))/6
        case (4,3):
            return 1

        case (5,1):
            Jbub1 = Jbub(1,beta)
            return Jbub1**2 - Jbub1 + Jbub(2,t)/2 + pi**2/12 + 1/3
        case (5,2):
            return (1 - Jbub(1,beta))/3
        case (5,3):
            return 1/3

        case (6,1):
            Jbub1 = Jbub(1,beta)
            return (4*Jbub1**2 + 10*Jbub1 + Jbub(2,beta) - 4)/(4*(t-4))
        case (6,2):
            return (2 + Jbub(1,beta))/(2*(t-4))

        case _:
            return 0

