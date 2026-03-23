"""
File: numerics.py
Author: Mattias Sjö

This is the main file, which implements the command-line interface.
It has no use when using this as a library.
"""

import re
import sys, os
import traceback
from datetime import datetime

# logging, including fallback to standard logging module
import logging
from HVPpy import clogger, init_clogging

from mpmath import mp, mpc, mpf

def get_with_default(args, defaults):
    """
    Parse a list of arguments, replacing missing ones by the provided defaults.

    Arguments:
    args - a list of arguments.
    defaults - a list of defaults values.

    Returns a list of values of length len(defaults). The first len(args) are
    the elements of args, cast to the type of the corresponding default. The
    remaining ones are taken directly from defaults. Elements of args beyond
    len(defaults) are ignored.
    """
    return (type(default)(args[i]) if i < len(args) else default for i, default in enumerate(defaults))

def match_prefix(key, pattern):
    if key == ''.join(s[0] for s in pattern.split('_')):
        return True
    ks = key.split('_')
    ps = pattern.split('_')
    if len(ks) != len(ps):
        return False
    for k, p in zip(ks, ps):
        if not p.startswith(k):
            return False
    return True

def get_via_prefix(key, options, name, *, case=lambda s: s, default=None):
    """
    Interpret an abbreviated name.

    Arguments:
    key - a string.
    options - a set of strings.
    name - a description of what the strings represent, used for error reporting.
    case - a function for casefolding and similar (default: identity function)
    default - a value to return if nothing matches. If callable, default(key) is returned instead. (default: raise an error)

    If case(key) is a prefix of exactly one element in options, return that element.
    Otherwise return default if applicable, otherwise raise a KeyError.
    """
    matches = [match for match in options if match_prefix(case(key), case(match))]
    if not matches:
        if default is not None:
            if callable(default):
                return default(key)
            else:
                return default
        else:
            raise KeyError(f"Unrecognized {name}: '{key}'")
    if len(matches) > 1:
        raise ValueError(f"Ambiguous {name} '{key}': could mean {', '.join(matches)}")
    return matches[0]

def find_via_prefix(pattern, strings, options, *, key=lambda x: x, case=lambda s: s):
    for string in strings:
        if case(pattern).startswith(case(key(string))):
            alts = [opt for opt in options if opt != pattern and match_prefix(case(key(string)), case(opt))]
            if alts:
                raise ValueError(f"""Ambiguous match: '{key}' could match '{pattern}' but also {', '.join(f"'{alt}'" for alt in alts)}""")
            return string
    return False


def get_from_enum(key, enum, name=None, case=lambda s: s.upper(), default=None):
    """
    Like get_via_prefix, but from the set of values of an enum.

    Arguments:
    key - a string.
    enum - an enum.
    name - a description of what the strings represent, used for error reporting (default: the lowercase form of the enum's class name).
    case - a function for casefolding and similar (default: identity function)
    default - a value to return if nothing matches. If callable, default(key) is returned instead. (default: raise an error)

    If key is the value of an element of the enum, return that element.
    Otherwise, see if it is a prefix of the name of any enum element as in get_via_prefix.
    """
    if name is None:
        name = enum.__name__.lower()
    for item in enum:
        if key == item.value:
            return item
    else:
        get_via_prefix(key, enum, name, case, default)

def parse_point(point):
    realnum = r"[0-9]*(?:\.[0-9]*)?(?:[eE][-+]?[0-9]+)?"
    if match := re.fullmatch(f"({realnum})([-+]{realnum})j", point):
        return mpc(mpf(match[1]), mpf(match[2]))
    elif match := re.fullmatch(f"([-+]?{realnum})j", point):
        return mpc(0, mpf(match[1]))
    elif match := re.fullmatch(f"[-+]?{realnum}", point):
        return mpf(match[0])
    else:
        raise ValueError(f"Malformed number: '{point}'")

points_default = [-2]
indices_default = [1,2,3,4,5,6]

def print_help():
    """ Print a helpful usage message. """

    from HVPpy import Method
    NL='\n'

    print(NL.join([
        "=== HVP, by Mattias Sjö 2026 ===",
        "",
        "Arguments may be given in any order, separated by whitepsace.",
        "Many arguments take one or more comma-separated values following '='.",
        "For instance, the arguments",
        "    index=1,2,3 method=0,e,s",
        "indicate that E1, E2 and E3 will be computed using the methods 0, e and s (see below).",
        "Whitespace following a comma is ignored.",
        "",
        "Arguments may be abbreviated as long as the abbreviation is unambiguous.",
        "Each part of a multi-word argument may be abbreviated separately, or may",
        "be entirely replaced by an acronym. Thus, timelike_with_contour may be",
        "expressed as either of time_w_cont, timelike, t_w_c, twc, etc.",
        "",
        "The main behavior of the program is controlled with the following arguments:",
        "  index=i,j,...   Set the list of indices, normally indicating Ei, Ej,...",
       f"                  (if omitted, defaults to 'index {','.join(str(x) for x in indices_default)}')",
        "  indices=i,j,... Alias of index",
        "  points=t,u,...  Set the points at which the integrals are to be computed.",
       f"                  (if omitted, defaults to 'point {','.join(str(x) for x in points_default)}')",
        "  methods=a,b,... Set the computation methods to be used.",
        "                  The methods, listed as '[short] long', are:"]
        + [f"                    [{method.value}] {method.name}" for method in Method]
        + [
        "                  Either the long (case-insensitive, may be abbreviated) or short name may be used.",
        "                  (if omitted, all methods are used)",
        "  plot=what,tmin,tmax,tres",
        "                  Produce data for plotting the values, from tmin to tmax in steps of 1/tres.",
       f"                 What (case-insensitive, may be abbreviated) may be any of{NL}{NL.join(f'{' '*20}{name:13s} {(NL+" "*35).join(descr)}' for name,descr in plotlist.items())}",
       f"                 The arguments default to {', '.join(str(x) for x in plot_defaults)}, respectively.",
        "  precision=n     Set the delimal precision of the calculations.",
        "  bits=n          Set the number of bits of precision in the calculations",
        "                  (bits and precision overwrite each other).",
        "  debug           Print much more information about the calculations.",
        "  profile         Run the program under cProfile to analyze its performance.",
        "                  It is recommended to pipe the output to a file for inspection.",
        "  arguments=file  Add the contents of file to the list of arguments.",
        "                  The profile and debug arguments don't work if loaded this way.",
        "  cache           Display the cached star-point integral values and terminate.",
        "  help            Print this message and terminate.",
        "",
        "There are many additional options, which affect the behavior of specific",
        "methods and may be implemented on the fly to add new features.",
        "An incomplete list of these is:",
       f"{NL.join(f'{'  '}{name:16s} {(NL+" "*19).join(descr)}' for name,descr in optionlist.items())}",
        "The options are also passed on to the mpmath library calls; see the",
        "documentation of quad and nsum in particular. Options not in this list",
        "will be accepted but may not be abbreviated.",
        "",
        "The plot argument causes the corresponding action(s) to be executed.",
        "Otherwise, the given integrals Ei are evaluated using the given methods at the given points."]))

def print_cache():
    import json
    from mpmath import workprec
    from HVPpy.Ebar import STAR_POINT_FILE
    from HVPpy.method import Method
    from HVPpy.clogging import ColorFormatter as cf
    from HVPpy.integration import QuadError
    from HVPpy.utilities import json_to_QE

    try:
        with open(STAR_POINT_FILE, 'r') as cache:
            cached = json.load(cache)
    except IOError:
        print(f"{cf.BOLD}{cf.COLOR%cf.RED}(no cache - exiting){cf.RESET}")
        return

    for m, mdata in cached.items():
        print(f"{cf.BOLD}{get_from_enum(m, Method).name}:{cf.RESET}")
        max_prec = max((int(p) for p in mdata.keys())) if mdata else 0
        for p, pdata in mdata.items():
            print(f"{' '*4}{cf.BOLD}{p}-bit precision:{cf.RESET}")
            for x in "JE":
                for n in "12":
                    with workprec(p):
                        print(f"{' '*8}{QuadError.highlight(
                            json_to_QE(pdata[x][n]),
                            prefix=f'I{x}({n}) = ',
                            reference=json_to_QE(mdata[str(max_prec)][x][n])) if p!=max_prec else True}")



plot_defaults = ("masters", -8., +32., 8.)
plotlist = {
    "masters":          ["Plot the selected masters (E_n for n in indices) over the given t range."],
    "Jbub":             ["Plot the selected Jbub functions (Jbub(n) for n in indices) over the given t range."],
    "Pi":               ["Plot the HVP components PiE, PiJ and PiZ over the given t range."],
    "t_theta":          ["Plot the mapping t->theta over the given t range."],
    "theta_t":          ["Plot the deviation from real t in the mapping theta->t over the given t range."],
    "tau_t_tau":        ["Plot the deviation from the identity when mapping tau->t and back,",
                            "with Im(tau) in the given range and Re(tau) between 0 and 1/2."],
    "t_tau_t":          ["Like tau_t_tau_err, but over the square in the complex t plane symmetric",
                            "across the real line and of the given real extent."],
    "images":           ["Plot the images of lines of constant Im(t) (one for each value in points)",
                            "under the various coordinate mappings involved in t->tau, over the given",
                            "range in Re(t)."],
    "hybrid":           ["Plot the various hybrid integrals involved in E5bar over the given range",
                            "in log10(chi), where chi is the crossover."],
    "diffeq":           ["Plot the deviation from satisfying the chosen integral's differential equation."],
    "errors":           ["Like masters, but also output a file detailing the distribution of the error budget",
                         "for E5 and E6."],
    }
optionlist = {
    "bar":              ["Compute Ebar also when index=1, 2 or 3 (default is to compute E_2d here)"],
    "2d":               ["Compute E_2d also when index=4, 5 or 6 (default is to compute Ebar here)"],
    "compare":          ["When computing master integrals with different methods, highlight",
                            "the output to show the degree of agreement with the result of the",
                            "first method."],
    "use_theta":        ["Perform computations in theta space when applicable."],
    "theta_point":      ["Interpret the point input, and the plot ranges when applicable,",
                            "as theta values rather than t values."],
    "timelike_with_contour": ["Use the contour integration approach also in the timelike subthreshold",
                            "(0 < t < 4) region."],
    "all_with_contour": ["Always use the contour integration approach."],
    "imaginary_error":  ["Treat the imaginary part as an error on functions known to be real"],
    "subtract_0":       ["Subtract the value at t=0 from the master integrals."],
    "IJ_xover":         ["Set the crossover (chi) for the IJ integral."],
    "IE_xover":         ["Set the crossover (chi) for the IE integral."],
    "adaptive_xover":   ["Adapt the crossover (chi) to the situation.",
                            "(Warning: does not work very well!)"],
    "xover_step":       ["Set the stepsize used with adaptive_xover"],
    "below_cut":        ["Obtain values below the branch cut, not above it.",
                            "Not possible for all methods."],
    "d_logt":           ["Set the number of log(t)-derivatives to apply to the integrals."],
    "no_cached":        ["Ignore precomputed integral values; always recompute everything."],
    "overwrite_cache":  ["Recompute precomputed integral values the first time they are needed."],
    "integrand_error":  ["Perform rough error propagation from integrands to integrals.",
                            "Makes affected integrals more expensive."],
    "tolerance":        ["Set the tolerance of certain library calls independently of the working precision."],
    "quad_method":      ["Set the integration method used by the quad library call."],
    "samples":          ["Set the number of samples used by the FeynTrop method."],
    "lambda":           ["Set the deformation parameter used by the FeynTrop method.", "This is needed above threshold and must neither be too large nor too small."],
    "max_series_order": ["Limit the order of series expansions."],
    "contour_lower_precision": ["Set a factor by which to lower the working precision",
                            "while computing contour integrals."],
    "timeout":          ["Set the maximum time (in seconds) spent computing any value.",
                            "Overtime computations are interrupted and return NaN."],
    }
arglist = ["index", "indices", "points",  "methods", "plot", "help", "profile", "debug", "precision", "bits", "cache", "arguments"] + list(optionlist.keys())


def get_args(argv):
    for item in re.sub(r",\s*", ",", ' '.join(argv)).split():
        if '=' in item:
            arg, vals = item.split('=', 1)
            yield arg, [v.strip() for v in vals.split(',')]
        else:
            yield item, []

def main():
    """
    Implement the command-line interface of the numerics.

    Read sys.argv as described by print_help() (accessed when "help" in sys.argv).
    """

    NL = '\n' # for use in f-strings

    # This is done as early as possible to catch all potential logging outputs
    if arg := find_via_prefix('debug', get_args(sys.argv), arglist,
                                 key=lambda x: x[0], case=lambda s: s.lower()):
        if len(arg[1]) > 1:
            raise ValueError(f"Argument 'debug' takes at most one value, got {','.join(arg[1])}")
        if arg[1]:
            init_clogging(arg[1][0])
        else:
            init_clogging(logging.DEBUG)
    else:
        init_clogging(logging.INFO)

    if len(sys.argv) < 2:
        print_help()
        print("\n(Help message printed since no arguments were given.)")
        return

    from HVPpy import Method

    indices = []
    points = []
    methods = []
    options = {}
    actions = {}

    def parse_args(argv):
        # Read all arguments and take the corresponding actions
        nonlocal indices, points, methods, options, actions
        for arg, vals in get_args(sys.argv):

            # Attempt to interpret abbreviations, otherwise just casefold
            try:
                arg = get_via_prefix(arg, arglist, "argument", case=lambda s: s.lower())
            except KeyError:
                arg = arg.lower()

            # Match the special arguments
            match arg:

                case "arguments":
                    for val in vals:
                        with open(val, 'r') as infile:
                            if parse_args(infile.read()):
                                return True

                case "index" | "indices":
                    try:
                        indices += [int(x) for x in vals]
                    except ValueError:
                        raise ValueError(f"Argument '{arg}' should be followed by a comma-separated list of integrers (indices of Ei, Jbub(i), ...)")

                case "points":
                    try:
                        points += [parse_point(x) for x in vals]
                    except ValueError:
                        raise ValueError(f"Argument '{arg}' should be followed by a comma-separated list of complex numbers")

                case "precision":
                    if len(vals) != 1:
                        raise ValueError(f"Argument {arg} takes exactly one value, got {','.join(vals)}")
                    mp.dps = int(vals[0])
                case "bits":
                    if len(vals) != 1:
                        raise ValueError(f"Argument {arg} takes exactly one value, got {','.join(vals)}")
                    mp.prec = int(vals[0])

                case "methods":
                    methods += [get_from_enum(x, Method) for x in vals]

                case "plot":
                    try:
                        what,tmin,tmax,tres = get_with_default(vals, plot_defaults)
                    except ValueError:
                        raise ValueError(f"Argument 'plot' should be followed by what,t_min,t_max,t_res (defaulting to {', '.join((str(x) for x in plot_defaults))} respectively)")

                    name = "plot_" + get_via_prefix(what, plotlist, "plot option", case=lambda s: s.lower())
                    actions[name] = (tmin,tmax,tres)

                    clogger.debug(f"Added action {name}={tmin},{tmax},{tres}")

                case "cache":
                    print_cache()
                    return True

                case "help":
                    print_help()
                    return True

                # These are already handled
                case "profile":
                    pass
                case "debug":
                    pass

                # Miscellaneous options
                case _:
                    if len(vals) > 1:
                        raise ValueError(f"Argument {arg} takes at most one value, got {','.join(vals)}")

                    if vals:
                        # Values are cast to int if possible, else float, else string
                        try:
                            options[arg] = int(vals[0])
                        except ValueError:
                            try:
                                options[arg] = float(vals[0])
                            except ValueError:
                                options[arg] = vals[0]
                    else:
                        options[arg] = True

                    clogger.debug(f"Added option {arg}={options[arg]}")
        return False

    if parse_args(sys.argv[1:]):
        return

    # Defaults in case of empty arguments
    if not points:
        points = points_default
    if not indices:
        indices = indices_default
    if not methods:
        methods = Method

    clogger.info(f"{points=}")
    clogger.info(f"{indices=}")
    clogger.info(f"{methods=}")
    clogger.info(f"options=[{', '.join(f'{k}={v}' for k,v in options.items())}]")
    clogger.info(mp)

    # Import most functionality AFTER setting all settings
    from HVPpy import output

    # Make output directory if needed
    if any(action.startswith("plot") for action in actions):
        os.makedirs(output.plot_dir(), exist_ok=True)

    # Default action
    if not actions:
        for t in points:
            output.evaluate_masters(t, indices, methods, **options)

    for action, args in actions.items():
        match action:
            case "diffeq":
                for t in points:
                    diffeq.run_diffeq_suite(t, *args, **options)

            case "plot_masters":
                output.plot_E_data(*args, indices, methods, **options)
            case "plot_Jbub":
                # This reinterprets indices as Jbub indices
                output.plot_Jbub_data(*args, indices)
            case "plot_Pi":
                output.plot_Pi_data(*args, methods, **options)
            case "plot_t_theta":
                output.plot_theta_data(*args, **options)
            case "plot_theta_t":
                output.plot_theta_t_err(*args, **options)
            case "plot_tau_t_tau":
                output.plot_tau_t_tau_err(*args, **options)
            case "plot_t_tau_t":
                output.plot_t_tau_t_err(*args, **options)
            case "plot_images":
                output.plot_t_tau_images(*args, points, **options)
            case "plot_hybrid":
                output.plot_hybrid_err(*args, points, **options)
            case "plot_diffeq":
                output.plot_diffeq(*args, indices, methods, **options)
            case "plot_errors":
                output.plot_E_data(*args, indices, methods, err_log=True, **options)

            case _:
                clogger.error(f"Unknown action: {action}")
                return


# Run the main method, possibly under a profiler, and highlight its errors
if __name__ == '__main__':
    with open(f"HVP.log", 'a') as log:
        print(f"[{datetime.now().strftime("%c")}] {' '.join(sys.argv)}", file=log)

    try:
        if arg := find_via_prefix('profile', get_args(sys.argv), arglist,
                           key=lambda x: x[0], case=lambda s: s.lower()):
            if arg[1]:
                raise ValueError(f"Argument 'profile' takes no values")
            import cProfile
            cProfile.run('main()')
        else:
            main()
    except Exception:
        logging.error(traceback.format_exc())
    except SystemExit:
        logging.error(traceback.format_exc())
        clogger.error("Process interrupted")
    except KeyboardInterrupt:
        logging.error(traceback.format_exc())
        clogger.error("Process interrupted by user")
