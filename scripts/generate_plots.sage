
# This SAGE program generates the domain-colored figures for the plots
# using SAGE's superior facilities for that.
# Call as sage generate_plots.sage TARGET[:RES] ...
#  which produces plots/TARGET.png (or variants thereof) at resolution RES
#  (each TARGET has a different default RES)
# The special TARGET "all" expands to all available targets
# It must be given access to parts of the main integral implementation.
# trim_plots.sh must be called afterward.

import sys

import HVPpy

from mpmath import mp, mpf
mp.prec = 128

from HVPpy import clogger, init_clogging
init_clogging()

verbose = True;
directory = '../plots'

def _safe(func, z, safe):
    global verbose;
    if verbose:
        if verbose % (res**2/64) == 0:
            print('.', end='', flush=True)
        verbose += 1
    if not safe:
        return complex(func(complex(z)))
    try:
        return complex(func(complex(z)))
    except Exception as err:
        clogger.warning(f"Ignoring {type(err)}: {err}")
        return NaN;

def _plot(func, re, im, name, *, plot_opt={}, pic_opt={}, safe=True):
    global verbose;

    clogger.info(f"Generating '{name}.png' at resolution {res}")
    if verbose:
        verbose = 1;
        print('[', end='', flush=True);

    save(complex_plot(lambda z: _safe(func, z, safe), re, im, contoured=True, plot_points=res, **plot_opt), f'{directory}/{name}.png', axes=False, **pic_opt)
    if(verbose):
        print(']', flush=True)

res = None
def _get_res(default):
    global res
    if not res:
        res = default

ctx = HVPpy.IntegrationContext(-2,None,None, method=HVPpy.Method.EISENSTEIN); # dummy context

targets = ["tau", "q", "tsun", "rho2", "rho4", "rhoF", "rhoF_inset", "tau_t_tau", "Hbball0", "Hbball1", "Hbball2", "Hbball3", "E1q", "E2q", "E3q", "w1", "Li0", "Li1", "Li2", "Li3"]

def main(argv):
    global res
    for target in argv:

        if (mark := target.find(':')) != -1:
            res = int(target[mark+1:])
            target = target[:mark]
        else:
            res = None

        match target:

            case "all":
                main([f"{t}:{res}" for t in targets] if res else targets)

            case "tau":
                _get_res(2048)
                eps = .5/res
                _plot(HVPpy.tau_to_t, (-.5,.5), (eps,.5), 'tau')

            case "q":
                _get_res(2048)
                _plot(lambda q: (HVPpy.tau_to_t(log(q)/(2j*pi)) if abs(q)<1 else NaN), (-1,1), (-1,1), 'q', pic_opt={'aspect_ratio':1})

            case "tsun":
                _get_res(2048)
                _plot(lambda ts: -64*ts/((ts-9)*(ts-1)), (-8,16), (-5,5), 'tsun')

            case "rho2":
                _get_res(1024)
                _plot(lambda z: sqrt((z-4)*(z-16)), (-4,24), (-6,6), 'rho2_before')
                _plot(lambda z: HVPpy.rho2_corrected(z, HVPpy.rho2_region(z)), (-4,24), (-6,6), 'rho2_after')

            case "rho4":
                _get_res(1024)
                _plot(lambda z: (z**3 - 9*z**2 + 3*z - 3)**(1/4), (-3,12), (-3,3), 'rho4_before')
                _plot(lambda z: HVPpy.rho4_corrected(z, HVPpy.rho4_region(z)), (-3,12), (-3,3), 'rho4_after')

            case "rhoF":
                _get_res(2048)
                _plot(lambda z: HVPpy.F1(1728/HVPpy.j_sun(z))-HVPpy.F1_at_1, (-3,12), (-3,0), 'rhoF_before')
                _plot(lambda z: HVPpy.rhoF_corrected(1728/j_sun(z), rhoF_region(z))-F1_at_1, (-3,12), (0,3), 'rhoF_after')

            case "rhoF_inset":
                _get_res(256)
                _plot(lambda z: HVPpy.F1(1728/HVPpy.j_sun(z))-HVPpy.F1_at_1, (8.5,9.25), (-.25,0), 'rhoF_before_inset')
                _plot(lambda z: HVPpy.rhoF_corrected(1728/HVPpy.j_sun(z), HVPpy.rhoF_region(z))-HVPpy.F1_at_1, (8.5,9.25), (0,.25), 'rhoF_after_inset')

            case "tau_t_tau":
                _get_res(2048)
                eps = .5/res
                _plot(lambda tau: (2**60*(HVPpy.t_to_tau(HVPpy.tau_to_t(tau))-tau) if (abs(tau-1/2)-eps > sqrt(3)/6 and abs(tau-1/6)-eps > 1/6) else float('nan')), (eps, .5-eps), (eps, .5), 'tau_t_tau')

            case "Hbball0" | "Hbball1" | "Hbball2" | "Hbball3":
                _get_res(512)
                mp.prec = 64
                r = int(target[-1])
                _plot(lambda q: (HVPpy.QuadError.decay(H_bball(r, q,log(q), ctx)) if abs(q)<1 else NaN), (-1,1), (-1,1), f'Hbball{r}', pic_opt={'aspect_ratio':1})

            case "E1q" | "E2q" | "E3q":
                _get_res(512)
                mp.prec = 64
                r = int(target[1])
                E0 = {1: -7*zeta(3), 2: 7*zeta(3)/4, 3: (6 - 35*zeta(3))/32}
                _plot(lambda q: (HVPpy.QuadError.decay(HVPpy.E_2d(r, tau=log(q)/(2j*pi), method=ctx.method) - E0[r]) if abs(q)<1 else NaN), (-1,1), (-1,1), f'E{r}q', pic_opt={'aspect_ratio':1})

            case "E4q":
                _get_res(512)
                mp.prec = 64
                r = int(target[1])
                E0 = -5/mpf(3) + (7*pi**2)/12 - zeta(3)
                _plot(lambda q: (HVPpy.QuadError.decay(HVPpy.Ebar(4, tau=log(q)/(2j*pi), method=ctx.method) - E0) if abs(q)<1 else NaN), (-1,1), (-1,1), f'E4q', pic_opt={'aspect_ratio':1})

            case "w1":
                _get_res(512)
                mp.prec = 64
                _plot(lambda q: (HVPpy.QuadError.decay(HVPpy.varpi_1(log(q)/(2j*pi), method=ctx.method)) if abs(q)<1 else NaN), (-1,1), (-1,1), f'w1', pic_opt={'aspect_ratio':1})

            case "Li0" | "Li1" | "Li2" | "Li3":
                _get_res(512)
                mp.prec = 64
                r = int(target[-1])

                def Li_theta(r, argz, theta):
                    from theta import theta_to_tau
                    logq = 2j*pi*HVPpy.theta_to_tau(theta)
                    z = exp(1j*pi*argz)

                    return HVPpy.QuadError.decay(HVPpy.Li_elliptic(r, exp(logq),logq, z))

                _plot(lambda zt: Li_theta(r, zt.imag, zt.real), (-4,20), (0,2), f'Li{r}', safe=False)

if __name__ == '__main__':
    main(sys.argv[1:])
    clogger.info(f"All requested plots generated in directory '{directory}'.")



# DEPRECATED by join_contours.py
# clogger.info(f"Generating 'abs_jsun_geq1.png'...")
# save(complex_plot(lambda z: (1 if abs(1728/j_sun(z)) > 1 else 0), (-3,12), (-3,3), plot_points=res), f'{directory}/abs_jsun_geq1.png', axes=False)
