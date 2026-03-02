
from heapq import merge
from itertools import product, groupby
from contextlib import nullcontext

from mpmath import mpf, mpc, mp
from mpmath import log, zeta, pi, inf, exp, factorial
from numpy import linspace

from .clogging import clogger, ColorFormatter
from .E_2d import E_2d
from .Ebar import Ebar
from .Jbub import Jbub
from .utilities import *
from .method import Method
from .t_tau_beta import t_to_beta, t_to_tau, tau_to_t
from .theta import t_to_theta, theta_to_tau, theta_to_t
from .integration import QuadError, IntegrationContext
from .points import *

def plot_linspace(x_min, x_max, x_res,
                  refine_pts=[], refine_lvl=3, refine_size=1, refine_pow=2,
                  over_theta=False):
    if refine_lvl <= 0 or not refine_pts:
        yield from linspace(x_min, x_max, int(1 + x_res*(x_max-x_min)))
        return

    # Multi-way merge of sorted iterables, collapsing identical elements to one
    def _merge_unique(*iterables):
        def _next(iterator):
            try:
                return next(iterator)
            except StopIteration:
                return None

        iterators = {i:v for i,v in {it: _next(it) for it in (iter(iterable) for iterable in iterables)}.items() if v is not None}

        while iterators:
            min_val = min(val for val in iterators.values() if val is not None)
            yield min_val

            iterators = {i:v for i,v in {it: (_next(it) if val == min_val else val) for it,val in iterators.items()}.items() if v is not None}

    # The whole range at standard resolution, plus each point ± refine_size
    #  at resoluion times refine_pow (and recursively so up to refine_lvl).
    yield from _merge_unique(
        plot_linspace(x_min, x_max, x_res),
        *(
            plot_linspace(point-refine_size, point+refine_size, x_res*refine_pow,
                          [point], refine_lvl-1, refine_size/refine_pow, refine_pow)
            if point != inf else
            sorted(1/x if x != 0. else inf # Special treatment for point at infinity
                for x in plot_linspace(-refine_size, +refine_size, x_res,
                                        [0.], refine_lvl, refine_size/refine_pow, refine_pow))
            for point in refine_pts if not ((point < x_min or point > x_max) or (over_theta and point == 0))
            )
        )

def plot_dir():
    return f"plots/prec{mp.prec}"
def plot_filename(name, method, suffix="", **options):
    return f"""{plot_dir()}/{name}{method.value}{
            options.get("max_series_order", '')
        }{
            'e' if "elliptic_e4" in options else ''
        }{suffix}.dat"""

def evaluate_masters(tt, masters, methods, **options):
    over_theta = options.get("theta_point", False)
    if over_theta:
        t = theta_to_t(tt)
        clogger.debug(f"{masters} theta={tt}, t={t}, tau={theta_to_tau(tt)}, beta={t_to_beta(t)}{', '.join(f'{o}={v}' for o,v in options.items())}")
    else:
        t = tt
        clogger.debug(f"{masters} t={t}, tau={t_to_tau(t)}, beta={t_to_beta(t)}{', '.join(f'{o}={v}' for o,v in options.items())}")

    for n in masters:
        if (
            (not "keep_problematic" in options)
            and
            (t in (problematic_with_theta if options.get("use_theta", False) else problematic)[n])
            ):
            clogger.info(f"Skipping problematic t value for E{n}: {t}")
            continue

        values = []
        for method in methods:

            if d_logt := options.get('d_logt', 0):
                prefix = f"(d/dlogt){f'^{d_logt}' if d_logt > 1 else ''} "
            else:
                prefix = ""

            if n >= 4 or "bar" in options:
                Efun, Ename = Ebar, f"{prefix}E{n}bar"
            else:
                Efun, Ename = E_2d, f"{prefix}E{n}_2d"

            if over_theta:
                var = theta_to_tau(tt)
                values.append((f"{Ename},{method.value}", Efun(n, tau=var, method=method, **options)))
            else:
                if is_real(t) and t.real > threshold[n] and method.needs_i_epsilon():
                    var = t + 1j*mp_eps() * (-1 if options.get("below_cut", False) else +1)
                else:
                    var = t
                values.append((f"{Ename},{method.value}", Efun(n, t=var, method=method, **options)))

        print()
        for i, (expr,val) in enumerate(values):
            print(QuadError.highlight(val,
                                      prefix=f"{expr} = ",
                                      reference=((values[0][1] if i else True) if options.get("compare", False) else False)))
    print()

def print_data_point(t, out):
    tsun = t_sun(t)
    rtsun = r_map(tsun)
    tau = t_to_tau(t)
    q = exp(2j*pi*tau)
    print('\t'.join(f"{x}" for x in (t.real, t.imag, tsun.real, tsun.imag, rtsun.real, rtsun.imag, tau.real, tau.imag, q.real, q.imag)), file=out)

def err_log_header():
    return TAB.join(("t", "boundary", "Irat", "IJ", "IJstar", "IJcont", "IE", "IEstar", "IEcont"))

def print_data_header(out):
    print('\t'.join(("Re(t)", "Im(t)", "Re(tsun)", "Im(tsun)", "Re(rtsun)", "Im(rtsun)", "Re(tau)", "Im(tau)", "Re(q)", "Im(q)")), file=out)

def plot_Jbub_data(t_min,t_max,t_res, ns):
    for n in ns:
        clogger.info(f"Plotting Jbub{n}/{n}!...")
        with open(f"{plot_dir()}/Jbub{n}.dat", 'w') as out:
            print('\t'.join(("t", f"ReJ{n}", f"ImJ{n}")), file=out)
            clogger.debug(f"Plotting Jbub{n}...")

            for t in plot_linspace(t_min, t_max, t_res, refine_pts=[4], refine_lvl=5):

                # if t == 4:
                #     print('', file=out)
                #     continue

                if t < 4:
                    beta = t_to_beta(t)
                else:
                    beta = t_to_beta(t+mp_eps())

                J = Jbub(n, beta) / factorial(n)

                print('\t'.join(str(x) for x in (t, J.real, J.imag)), file=out)
                clogger.debug(f"{t}: {J}")

def _quotient(num, den, where=""):
    if QuadError.decay(den) != 0:
        return num/den

    # clogger.debug(f"Elided division by zero{f' ({where})' if where else ''}: {num} / {den}")

    if QuadError.decay(num) == 0:
        return 0
    else:
        return 1

def _complex_val(value):
    val, _ = QuadError.val_err(value)
    return mpc(val.real, val.imag)
def _complex_abs(value):
    _, err = QuadError.val_err(value)
    return mpc(abs(err.real), abs(err.imag))
def _complex_rel(value):
    val, err = QuadError.val_err(value)
    return mpc(
        abs(_quotient(err.real, val.real)),
        abs(_quotient(err.imag, val.imag))
        )
def _complex_min(value):
    val, err = QuadError.val_err(value)
    return mpc(val.real - abs(err.real), val.imag - abs(err.imag))
def _complex_max(value):
    val, err = QuadError.val_err(value)
    return mpc(val.real + abs(err.real), val.imag + abs(err.imag))

@log_errors
def plot_E_data(t_min,t_max,t_res, ns, methods, err_log=False, **options):

    reference = {}
    use_theta = options.get("use_theta", False)
    over_theta = options.get("theta_point", False)

    header = f'''{
            TAB.join(f"{norm}{err}{reim}"
                    for norm,reim,err in product(
                        ("", "Norm"),
                        ("Re", "Im", "ReX"),
                        ("", "Abs", "Rel", "Min", "Max")
                        )
                    )
        }{TAB}{"Time"}'''

    for n in ns:
        for meth in methods:
            E_func = E_2d if n < 4 else Ebar

            with (
                open(plot_filename(f"E{n}", meth, **options), 'w') as out,
                (open(plot_filename(f"E{n}", meth, "_err", **options), 'w') if err_log else nullcontext()) as err
                ):

                print(f"t{TAB}{header}", file=out)
                if err:
                    print(err_log_header(), file=err)

                clogger.info(f"Plotting E{n} ({meth.value})...")

                problem_points = (problematic_with_theta if use_theta else problematic)[n]

                for tt in plot_linspace(t_min, t_max, t_res, refine_pts=problem_points, over_theta=over_theta):

                    if tt in problem_points:
                        clogger.debug(f"Skipping problematic t value for E{n}: {tt}")
                        print('', file=out)
                        continue

                    if n > 4 or options.get("bar", False):
                        Efun = Ebar
                    else:
                        Efun = E_2d

                    if over_theta:
                        tau = theta_to_tau(tt)
                        t = tau_to_t(tau)
                        clogger.debug(f"Plotting θ = {tt} (t = {t})")

                        result = timed(E_func)(n, t=t, tau=tau, method=meth, err_log=err, **options)
                    else:
                        if tt > threshold[n] and meth.needs_i_epsilon():
                            t = tt + 1j*mp_eps() * (-1 if options.get("below_cut", False) else +1)
                        else:
                            t = tt
                        clogger.debug(f"Plotting t = {t}")
                        result = timed(E_func)(n, t=t, method=meth, err_log=err, **options)

                    value, time = result

                    if meth == methods[0]:
                        reference[(n,tt)] = value.value()

                    norm_value = (
                        +    (_quotient(value.real, reference[(n,tt)].real, f'{tt},re') - 1.)
                        + 1j*(_quotient(value.imag, reference[(n,tt)].imag, f'{tt},im') - 1.)
                        )

                    clogger.debug(f"values:{NL} > {value}")
                    clogger.debug(f"(norm):{NL} > {value}")

                    # Note the parallellism with how the header is constructed!
                    line = f'''{
                            TAB.join(str(Re(err(reim(norm))))
                                for norm,reim,err in product(
                                    (value, norm_value),
                                    (lambda x: x.real, lambda x: x.imag, lambda x: x.real.add_error(x.imag) if tt < threshold[n] else x.real),
                                    (lambda x: x.value(), _complex_abs, _complex_rel, _complex_min, _complex_max)
                                    )
                                )
                        }{TAB}{time}'''

                    print(f"{t.real}{TAB}{line}", file=out)

                clogger.info(f"Wrote {out.name}")
    clogger.debug(f"Columns are {header.replace('\t', ', ')}")


def plot_Pi_data(t_min,t_max,t_res, methods, **options):

    ctx = IntegrationContext(1,None,None, **options)
    Lp = 2*log(135/770)

    def PiJEZ(JEZ, ctx):
        clogger.debug(ctx)
        t = ctx.t
        terms = {}

        match JEZ:
            case 'J':
                beta = ctx.beta
                terms[3] = -(t**2/648 + t/324 - 1/108 - 1/(9*t))*Jbub(3,beta)
                terms[21] = -(t**2/54 - 8*t/27 + 4/3 - 2/(3*t))*Jbub(2,beta)*Jbub(1,beta)
                terms[111] = -(t**2/54 - 2*t/9 + 8/9 - 32/(27*t))* Jbub(1,beta)**3
                terms[2] = -(5*t**2/216 - 103*t/108 + 151/36 + 1/(3*t))*Jbub(2,beta)
                terms[11] = -((1+ Lp)*t**2/27 - (63+4*Lp)*t/54 + (275-64*Lp)/54 + (10 + 96*Lp)/(27*t))*Jbub(1,beta)**2
                terms[1] = +((97 - 9*pi**2)*t**2/3888 - (2285 + 9*pi**2)*t/1944 + 407/81+pi**2/72 - (4-pi**2)/(6*t))*Jbub(1,beta) - (t**2*Lp/27 + (23 + 8*Lp)*t/27 - (92 + 48*Lp)/27)*Lp*Jbub(1,beta)
                terms[0] = - (4-pi**2)/36

            case 'E':
                terms[1] = (107*t**3/46656 - 233*t**2/3888 - 47*t/243 + 17327/1458 - 20813/(486*t) + 1952/(81*t**2))*E_2d(1,ctx)
                terms[2] = -(143*t**4/46656 - 257*t**3/1944 + 2209*t**2/1944 + 32287*t/2916 - 32126/243 + 72214/(243*t) - 7808/(81*t**2))*E_2d(2,ctx)
                terms[3] = -(143*t**4/23328 - 166*t**3/729 + 316*t**2/243 + 29413*t/1458 - 99722/729 + 5216/(27*t))*E_2d(3,ctx)
                terms[5] = -(t/3 - 13/9 - 2/(3*t))*Ebar(5,ctx)
                terms[6] = (t**3/54 - 19*t**2/27 + 143*t/27 - 100/9)*Ebar(6,ctx)
                terms[0] = +(1216/27+2*zeta(3)/9-pi**2/18)*1/t - (26885/972-2131*zeta(3)/108+91*pi**2/108)

            case 'Z':
                terms[0] =  ( - t*zeta(3)/9 - (t**2 + 18*t)/81*Lp**3 - 20*t/27*Lp**2 - (17*t**2/432 - 845*t/648)*Lp - (t**2/216 - 5*t/27)*pi**2 + 20245*t**2/139968 - 55511*t/23328 )

        return {'': sum(term for term in terms.values())} | terms


    for JEZ in 'JEZ':
        terms0 = PiJEZ(JEZ, ctx.but(method=Method.EXPANSION_0))

        with open(f"{plot_dir()}/Pi{JEZ}.dat", 'w') as out:
            header = f'''{
                    TAB.join(f"{meth.value}{err}{reim}Pi{JEZ}{term}"
                            for meth,reim,err,term in product(
                                methods,
                                ("Re", "Im"),
                                ("", "Abs", "Rel", "Min", "Max"),
                                terms0.keys()
                                )
                            )
                }{TAB}{
                    TAB.join(f"{meth.value}Time" for meth in methods)
                }'''

            print(f"t{TAB}{header}", file=out)
            clogger.info(f"Plotting Pi{JEZ}...")

            for t in plot_linspace(t_min, t_max, t_res, refine_pts=[0,4,16], refine_lvl=5):

                clogger.debug(f"Plotting t = {t}")

                if t in {0,4,16}:
                    clogger.debug(f"Skipping problematic t value for Pi{JEZ}: {t}")
                    print('', file=out)
                    continue

                    if tt > 4 and meth.needs_i_epsilon():
                        t += 1j*mp_eps() * (-1 if options.get("below_cut", False) else +1)

                results = [timed(PiJEZ)(JEZ, ctx.but(t=t, method=meth)) for meth in methods]
                values = [val for val,_ in results]
                times = [time for _,time in results]

                def _quotient(num, den):
                    return num/den if den != 0 else 0 if num == 0 else float('nan')

                clogger.debug(f"values:{NL}{NL.join(f' > {v['']}' for v in values)}")

                # Note the parallellism with how the header is constructed!
                line = f'''{
                        TAB.join(str(Re(err(reim(values[meth][key]))))
                            for meth,reim,err,key in product(
                                range(len(methods)),
                                (lambda x: x.real, lambda x: x.imag),
                                (lambda x: QuadError.get_value(x), _complex_abs, _complex_rel, _complex_min, _complex_max),
                                terms0.keys()
                                )
                            )
                    }{TAB}{
                        TAB.join(str(times[meth]) for meth in range(len(methods)))
                    }'''

                print(f"{t.real}{TAB}{line}", file=out)

        clogger.info(f"Wrote {out.name}")
    clogger.debug(f"Columns are {header.replace('\t', ', ')}")

def plot_theta_data(t_min,t_max,t_res, tol=1e-9):
    with open(f"{plot_dir()}/theta.dat", 'w') as out:
        print('\t'.join(("t", "theta", "err")), file=out)
        for t in plot_linspace(t_min, t_max, t_res):
            theta = t_to_theta(t, tol)
            print('\t'.join(str(x) for x in (t, theta.value(), theta.error())), file=out)


def plot_theta_t_err(t_min, t_max, theta_res, **options):
    theta_min = t_to_theta(t_min, error=False)
    theta_max = t_to_theta(t_max, error=False)

    with open(f"{plot_dir()}/theta_t_err.dat", 'w') as out:
        print('\t'.join(("theta", "Retau", "Imtau", "Ret", "AbsImt", "eRet", "eAbsImt", "eReErrt", "eImErrt")), file=out)
        for theta in plot_linspace(theta_min, theta_max, theta_res):
            if theta == 0:
                continue
            tau = theta_to_tau(theta)
            try:
                t = tau_to_t(tau, error=False) # NOTE theta_to_t discards the imaginary part
            except ValueError:
                t = mpc('nan', 'nan')
            try:
                et, det = tau_to_t(tau, method=Method.EISENSTEIN, error=True)()
            except ValueError:
                et = det = mpc('nan', 'nan')
            print('\t'.join(str(x) for x in (theta, tau.real, tau.imag, t.real, abs(t.imag), et.real, abs(et.imag), abs(det.real), abs(det.imag))), file=out)
        clogger.info(f"Wrote {out.name}")


def plot_t_tau_t_err(t_min, t_max, t_res, **options):
    from math import sqrt,atan2

    imsize = abs(t_max - t_min)/2

    with open(f"{plot_dir()}/t_tau_t_err.dat", 'w') as out:
        print('\t'.join(("Ret", "Imt", "Retau", "Imtau", "ReErr", "ImErr", "AbsErr", "ArgErr")), file=out)
        for Re_t in plot_linspace(t_min, t_max, t_res):
            for Im_t in plot_linspace(-imsize, +imsize, t_res):
                t = mpc(Re_t, Im_t)
                tau = t_to_tau(t)
                err = tau_to_t(tau) - t
                # try:
                #     tau = t_to_tau(t)
                #
                #     try:
                #         err = tau_to_t(tau) - t
                #     except ValueError:
                #         err = mpc('nan', 'nan')
                # except ValueError:
                #     tau = err = mpc('nan', 'nan')

                print('\t'.join(str(x) for x in (t.real, t.imag, tau.real, tau.imag, err.real, err.imag, abs(err), atan2(err.imag, err.real))), file=out)
            print('', file=out)
        clogger.info(f"Wrote {out.name}")

def all_about_tau_header():
    return [
        "Retau", "Imtau",
        "Ret", "Imt",
        "Retau", "Imtau",
        "Req", "Imq",
        "Rezc", "Imzc",
        "Rezr", "Imzr",
        "Rezc_alt", "Imzc_alt",
        "Rezr_alt", "Imzr_alt",
        "reg2",
        "reg4c", "regFc",
        "reg4r", "regFr"]

def all_about_tau(tau):
    from .t_tau_beta import t_sun,r_map, rho2_region,rho4_region,rhoF_region

    t = tau_to_t(tau)
    q = exp(2j*pi*tau)

    zc = t_sun(t)
    zc_alt = t_sun(t, alternative=True)

    zr = r_map(zc)
    zr_alt = r_map(zc_alt)

    return [
        tau.real, tau.imag,
        t.real, t.imag,
        tau.real, tau.imag,
        q.real, q.imag,
        zc.real, zc.imag,
        zr.real, zr.imag,
        zc_alt.real, zc_alt.imag,
        zr_alt.real, zr_alt.imag,
        rho2_region(t),
        rho4_region(zc), rhoF_region(zc),
        rho4_region(zr), rhoF_region(zr)]

def plot_tau_t_tau_err(_, tau_max, tau_res, **options):
    from math import sqrt,atan2

    eps = .5/tau_res

    with open(f"{plot_dir()}/tau_t_tau_err.dat", 'w') as out:
        print('\t'.join(all_about_tau_header() + ['ReErr', 'ImErr', 'AbsErr', 'ArgErr']), file=out)
        for Im_tau in plot_linspace(0.+eps, tau_max, tau_res * .5/tau_max):
            for Re_tau in plot_linspace(0.+eps, .5-eps, tau_res):
                tau = mpc(Re_tau, Im_tau)
                # Terminate scanline if we hit the upper circular arc
                if abs(tau - 1/2)-eps < sqrt(3)/6:
                    break
                # Skip values under the lower circular one
                if abs(tau - 1/6)-eps <= 1/6:
                    data = [mpf('nan')] * len(all_about_tau_header())
                    err = mpc('nan', 'nan')
                else:
                    data = all_about_tau(tau)
                    t = mpc(data[2], data[3])
                    try:
                        err = t_to_tau(t) - tau
                    except (ValueError, ArithmeticError) as err:
                        clogger.warning(f"Invalid result at tau={tau}: {err}")
                        err = mpc('nan', 'nan')
                print('\t'.join(str(x) for x in data + [err.real, err.imag, abs(err), atan2(err.imag, err.real)]), file=out)
            print('', file=out)
        clogger.info(f"Wrote {out.name}")
    clogger.debug(f"Columns are {all_about_tau_header.replace('\t', ', ')}")

def plot_t_tau_images(Re_min, Re_max, Re_res, Ims, **options):

    for Im_t in Ims:
        with open(f"{plot_dir()}/image_Imt{Im_t.real}.dat", 'w') as out:
            print('\t'.join(all_about_tau_header()), file=out)
            for Re_t in plot_linspace(Re_min, Re_max, Re_res, [0,4,16], refine_lvl=10):
                if Re_t == 10: # Avoid the exact transition for rho2 since there is roundoff error
                    Re_t += .1/Re_res

                t = mpc(Re_t, Im_t.real)
                tau = t_to_tau(t)
                data = all_about_tau(tau)
                data[2] = Re_t
                data[3] = Im_t.realq

                print('\t'.join(str(x) for x in data), file=out)
            clogger.info(f"Wrote {out.name}")
    clogger.debug(f"Columns are {all_about_tau_header.replace('\t', ', ')}")

def plot_hybrid_err(log_xover_min, log_xover_max, log_xover_res, points, **options):
    from .integration import hybrid_error_log_header


    for point in points:
        filename = f"{plot_dir()}/hybrid_err_t{point.real}"
        with open(f"{filename}.dat", 'w') as out:
            print(hybrid_error_log_header(), file=out)
            for log_xover in plot_linspace(max(log_xover_min, log(t_to_beta(point), 10)), log_xover_max, log_xover_res):
                xover = 10**Re(log_xover)
                clogger.debug(f"t={point}, chi={xover}")
                Ebar(5, t=point,
                     hybrid_error_log=out,
                     SJ_xover=xover, E1h2_xover=xover, Hreg_xover=xover,
                     **options)
            clogger.info(f"Wrote {out.name}")

        with open(f"{filename}_ref.dat", 'w') as out:
            Ebar(5, t=point,
                hybrid_error_log=out,
                adaptive_xover=True,
                **options)
            clogger.info(f"Wrote {out.name} (temporary file)")

        with open(f"{filename}_ref.dat", 'r') as ref:
            reference = {}
            for line in ref:
                cols = line.strip().split('\t')
                reference[cols[0]] = mpf(cols[2])
        outfile = {}
        with (
            open(f"{filename}.dat", 'r') as infile,
            open(f"{filename}_SJ1.dat", 'w') as outfile["SJ1"],
            open(f"{filename}_SJ2.dat", 'w') as outfile["SJ2"],
            open(f"{filename}_Hreg.dat", 'w') as outfile["Hreg"],
            open(f"{filename}_E1h2.dat", 'w') as outfile["E1h2"]
            ):
            for line in infile:
                cols = line.strip().split('\t')
                key = cols[0]
                if key == "source":
                    for key in ("SJ1", "SJ2", "Hreg", "E1h2"):
                        print(f"{cols[1]}{TAB}reldiff{TAB}{TAB.join(cols[2:])}", file=outfile[key])
                else:
                    print(f"{cols[1]}{TAB}{mpf(cols[2])/reference[key] - 1}{TAB}{TAB.join(cols[2:])}", file=outfile[key])
            for out in outfile.values():
                clogger.info(f"Wrote {out.name}")

    clogger.debug(f"Columns are {hybrid_error_log_header().replace('\t', ', ')}")

def plot_diffeq(t_min, t_max, t_res, ns, methods, **options):

    for n in ns:
        match n:
            case 1:
                order = 3
                problem_points = {0,4,16}
            case 5:
                order = 2
                problem_points = {0,4}
            case _:
                raise NotImplementedError(f"Not implemented: diffeq for E{n}")

        header = f"t{TAB}" + TAB.join((f'{reim}{err}{val}'
            for reim,err,val in product(
                ("Re", "Im"),
                ("", "Abs", "Rel", "Min", "Max"),
                ["Res"] + [f'D{d}' for d in range(order+1)]
            )))
        for meth in methods:
            with open(f"{plot_dir()}/diffeqE{n}{meth.value}.dat", 'w') as out:
                print(header, file=out)

                clogger.info(f"Evaluating diffeq for E{n} ({meth.value})...")

                prev_tt = t_min
                for tt in plot_linspace(t_min, t_max, t_res, refine_pts=problem_points, refine_lvl=5):

                    if is_real(tt) and tt.real > threshold[n]:
                        t = tt + 1j*mp_eps() * (-1 if options.get("below_cut", False) else +1)
                    else:
                        t = tt

                    if tt in problem_points:
                        if t == threshold[n]:
                            # If threshold is problematic, it is actually a singularity and should be skipped
                            clogger.debug(f"Skipping problematic t value for E{n}: {t}")
                            print('', file=out)
                            continue
                        else:
                            # Otherwise, dodge the removable singularity
                            clogger.debug(f"Avoiding problematic t value for E{n}: {t}")
                            t -= (tt - prev_tt)/2


                    clogger.debug(f"Evaluating at t = {t}")

                    ctx = IntegrationContext(t,None,None, method=meth, **options)

                    match n:
                        case 1:
                            # Use caching to avoid duplicate labor
                            if meth.needs_IBP_derivs():
                                ctx.cache[(2,0)] = E_2d(2, ctx)
                                ctx.cache[(3,0)] = E_2d(3, ctx)
                            for d in (0,1,2,3):
                                ctx.cache[(1,d)] = E_2d(1, ctx, d_logt=d)

                            res = (
                                + ctx.cache[(1,3)] * (t-16)*(t-4)*t**2
                                + ctx.cache[(1,2)] * 6*t*(t**2 - 15*t + 32)
                                + ctx.cache[(1,1)] * (7*t**2 - 68*t + 64)
                                + ctx.cache[(1,0)] * (t - 4)
                                - 24)

                        case 5:
                            ctx.cache[(1,0)] = E_2d(2, ctx)
                            if meth.need_IBP_derivs():
                                ctx.cache[(2,0)] = E_2d(2, ctx)
                                ctx.cache[(3,0)] = E_2d(3, ctx)
                            for d in (0,1,2):
                                ctx.cache[(5,d)] = Ebar(1, ctx, d_logt=d)

                            # TODO

                    line = f"{tt}{TAB}" + TAB.join((str(reim(err(val)))
                        for reim,err,val in product(
                            (lambda x: x.real, lambda x: x.imag),
                            (_complex_val, _complex_abs, _complex_rel, _complex_min, _complex_max),
                            [res] + [ctx.cache[(n,d)] for d in range(order+1)]
                        )))
                    print(line, file=out)

                    prev_tt = tt;

            clogger.info(f"Wrote {out.name}")
    clogger.debug(f"Columns are {header.replace('\t', ', ')}")

# DEPRECATED everything below here

# The default tmin/tmax values correspond to Im(tau)=1/2, Im(t)=0
def plot_t_tau_data(eps, tmin=-17.7355421803587, tmax=29.8564064605510, npoints=256, refine=0, refine_radius=.1, filename=None):
    # print(f"tmin={tmin}, tmax={tmax}")
    if not filename:
        filename = f"{plot_dir()}/eps{eps:.5f}.dat"
    with open(filename, 'w') as out:
        print_data_header(out)

        Ret = tmin
        step = (tmax-tmin)/npoints
        while Ret < tmax:
            print_data_point(Ret + 1j*eps, out)
            scale = 1
            while scale <= refine and min(abs(Ret), abs(Ret-4), abs(Ret-16)) < refine_radius**scale:
                scale += 1

            Ret += step / (2**scale)

    print(f"Generated plot data to '{filename}'")

def plot_bridge_data(t, spacelike=False, eps=1e-6, npoints=256, filename=None):
    if not filename:
        filename = f"{plot_dir()}/t{t:.5f}_{'s' if spacelike else 't'}bridge.dat"
    with open(filename, 'w') as out:
        print_data_header(out)
        t += 1j*eps

        if spacelike:
            origin = 1
            rad = np.absolute(t-1)
            argmax = np.angle(t-1)
            argmin = pi-argmax
        else:
            origin = 0
            rad = np.absolute(t)
            argmin = np.angle(t)
            argmax = pi/2

        for arg in np.linspace(argmin, argmax, npoints):
            print_data_point(beta_to_t(origin + rad*exp(1j*arg)), out)

    print(f"Generated {'spacelike' if spacelike else 'timelike'} bridge data to '{filename}'")

def plot_t0_data(t=-1, scale=.9, steps=100, eps=1e-6, filename=None):
    if not filename:
        filename = f"{plot_dir()}/t0.dat"

    with open(filename, 'w') as out:
        print('\t'.join(("t", "tau", "b", "l", "e", "p")), file=out, end='')

        t += 1j*eps
        for _ in range(steps):
            tau = t_to_tau(t)
            print(f"{t:.5f}" + '\t' + f"{tau:.5f}", file=out, end='')
            for method in "blep":
                try:
                    print("\t" + f"{E_2d(1,tau, method):.5f}", file=out)
                except Exception:
                    print("\tnan", file=out, end='')

            t *= scale

