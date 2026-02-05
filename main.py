"""
File: numerics.py
Author: Mattias Sjö

This is the main file, which implements the command-line interface.
It has no use when using this as a library.
"""

import sys
import traceback

# logging, including fallback to standard logging module
import logging
from HVPpy import clogger, init_clogging

from mpmath import mp

def main():
    """
    Implement the command-line interface of the numerics.

    Read sys.argv as described by print_help() (accessed when "help" in sys.argv).
    """

    global epsilon, tolerance

    NL = '\n' # for use in f-strings

    # This is done as early as possible to catch all potential logging outputs
    init_clogging(logging.DEBUG if "debug" in sys.argv else logging.INFO)


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

    def get_via_prefix(key, options, name, case=lambda s: s, default=None):
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
        matches = [match for match in options if case(match).startswith(case(key))]
        if not matches:
            if default is not None:
                if callable(default):
                    return default(key)
                else:
                    return default
            else:
                raise KeyError(f"Unrecognized {name}: '{key}'")
        if len(matches) > 1:
            raise KeyError(f"Ambiguous {name} '{key}': could mean {', '.join(matches)}")
        return matches[0]

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

    def print_help():
        """ Print a helpful usage message. """

        print('\n'.join([
            "Arguments may be given in any order, with most followed by comma-separated values.",
            "For instance, the arguments",
            "    index 1,2,3 method 0,e,s",
            "indicate that E1, E2 and E3 will be computed using the methods 0, e and s (see below).",
            "Arguments are case-insensitive and may be abbreviated."
            "",
            "The following arguments are supported:",
            "  index i,j,...   Set the list of indices, normally indicating Ei, Ej,...",
            f"                  (if omitted, defaults to 'index {','.join(str(x) for x in indices)}')",
            "  indices i,j,... Alias of index",
            "  points t,u,...  Set the points at which the integrals are to be computed.",
            f"                  (if omitted, defaults to 'point {','.join(str(x) for x in points)}')",
            "  methods a,b,... Set the computation methods to be used.",
            "                  The methods, listed as '[short] long', are:"]
            + [f"                     [{method.value}] {method.name}" for method in Method]
            + [
            "                  Either the long (case-insensitive, may be abbreviated) or short name may be used.",
            "                   (if omitted, all methods are used)",
            # " contour a,b,... Sets the contour to be used for E5 or E6.",
            # "                 The contours, listed as '[short] long', are:"]
            # + [f"                     [{contour.value}] {contour.name}" for contour in Contour]
            # + [
            # "                 Any contour may be followed by '=h' to set the contour height (only relevant for some contours)",
            # f"                  (if omitted, contour {Contour.FIXED_POINT.value} is used)",
            "  options a,b,... Set any other options. Each one may be a name=value pair, or just a name,",
            "                   in which it is interpreted as name=True",
            "                  These options are passed as **kwargs to most methods, making it easy to add new ones on the fly.",
            "                  An incomplete list is:",
            "                     use_theta               Use the theta-parametrization for integrals of E1 along the real line",
            "                     timelike_with_contour   Use contour integration also for 0 < t < 4",
            "                     all_with_contour        Use contour integration for all t",
            "                     no_cached               Recompute precomputed integral values",
            "                  Methods are taken in pairs, the first for the left-hand side (must support derivatives)",
            "                    and the second for the right-hand side (no derivatives).",
            "  plot what,tmin,tmax,tres",
            "                  Produce data for plotting the values, from tmin to tmax in steps of 1/tres.",
            f"                  What (case-insensitive, may be abbreviated) may be any of{NL}{NL.join(f'{' '*20}{name:13s} {(NL+" "*35).join(descr)}' for name,descr in plots.items())}",
            f"                  The arguments default to {', '.join(str(x) for x in plot_defaults)}, respectively.",
            "  precision n     Set the delimal precision of the calculations",
            "  bits n          Set the number of bits of precision in the calculations",
            "                   (bits and precision overwrite each other)",
            "  profile         Run the program under cProfile to analyze its performance.",
            "                  It is recommended to pipe the output to a file for inspection.",
            "  debug           Print much more information about the calculations.",
            "  help            Print this message and terminate.",
            "",
            "The arguments diffeq and plot cause the corresponding actions to be executed.",
            "Otherwise, the given integrals Ei are evaluated using the given methods at the given points."]))

    if len(sys.argv) < 2:
        print_help()
        print("\n(Help message printed since no arguments were given.)")
        return

    from HVPpy import Method

    # Plotting options
    plot_defaults = ("E", -8., +32., 8.)
    plots = {"E":           ["Plot the selected masters (E_n for n in indices) over the given t range."],
             "Jbub":        ["Plot the selected Jbub functions (Jbub(n) for n in indices) over the given t range."],
             "Pi":          ["Plot the HVP components PiE, PiJ and PiZ over the given t range."],
             "t_theta":     ["Plot the mapping t->theta over the given t range."],
             "theta_t":     ["Plot the deviation from real t in the mapping theta->t over the given t range."],
             "tau_t_tau":   ["Plot the deviation from the identity when mapping tau->t and back,", "with Im(tau) in the given range and Re(tau) between 0 and 1/2."],
             "t_tau_t":     ["Like tau_t_tau_err, but over the square in the complex t plane symmetric", "across the real line and of the given real extent."],
             "images":      ["Plot the images of lines of constant Im(t) (one for each value in points)", "under the various coordinate mappings involved in t->tau, over the given range in Re(t)."],
             "hybrid":      ["Plot the various hybrid integrals involved in E5bar over the given range", "in log10(chi), where chi is the crossover."],
             "diffeq":      ["Plot the deviation from satisfying the chosen integral's differential equation."]}

    # Default values for various arguments
    indices = [1]
    points = [-1]
    methods = Method
    options = {}
    actions = {}
    mp.dps = 20
    tolerance = 10**(3-mp.dps)

    # Read all arguments and take the corresponding actions
    try:
        argv = iter(sys.argv)
        next(argv)
        while True:
            arguments = ["index", "indices", "points",  "methods", "options", "plot", "help", "profile", "debug", "precision", "bits"]
            arg = get_via_prefix(next(argv), arguments, "argument", lambda s: s.lower())

            match arg:
                case "debug":
                    pass # this is already handled
                case "index":
                    try:
                        indices = [int(x) for x in next(argv).split(',')]
                    except ValueError:
                        raise ValueError(f"Argument 'index' should be followed by a comma-separated list of numbers (indices of Ei, Jbub(i), ...)")

                case "points":
                    try:
                        points = [complex(x) for x in next(argv).split(',')]
                    except ValueError:
                        raise ValueError(f"Argument 'point' should be followed by a comma-separated list of t values")

                case "precision":
                    mp.dps = int(next(argv))
                    tolerance = 10**(3-mp.dps)
                case "bits":
                    mp.prec = int(next(argv))
                    tolerance = 10**(3-mp.dps)

                case "methods":
                    methods = [get_from_enum(x, Method) for x in next(argv).split(',')]

                case "options":
                    for opt in next(argv).split(','):
                        if '=' in opt:
                            name, val = opt.split('=', maxsplit=1)
                            if '=' in val:
                                raise ValueError(f"More than one '=' not permitted in option, got '{opt}'")
                        else:
                            name, val = opt, True

                        # Values are cast to int if possible
                        try:
                            options[name] = int(val)
                        except ValueError:
                            options[name] = val

                    clogger.debug(f"Added option {name}={val}")

                case "profile":
                    pass # implemented outside main()

                case "plot":
                    try:
                        what,tmin,tmax,tres = get_with_default(next(argv).split(','), plot_defaults)
                    except StopIteration:
                        what,tmin,tmax,tres = plot_defaults

                    except ValueError:
                        raise ValueError(f"Argument 'plot' should be followed by what,t_min,t_max,t_res (defaulting to {', '.join((str(x) for x in plot_defaults))} respectively)")

                    name = "plot_" + get_via_prefix(what, plots, "plot option", lambda s: s.lower())
                    actions[name] = (tmin,tmax,tres)

                    clogger.debug(f"Added action {name}={tmin},{tmax},{tres}")

                case "help":
                    print_help()
                    return

    # This signals the end of sys.argv
    except StopIteration:
        pass

    clogger.debug(mp)


    # Import most functionality AFTER setting all settings
    from HVPpy import output
    #from twoloop import twoloop

    # Default action
    if not actions:
        for t in points:
            output.evaluate_masters(t, indices, methods, **options)

    for action, args in actions.items():
        match action:
            case "diffeq":
                for t in points:
                    diffeq.run_diffeq_suite(t, *args, **options)

            case "plot_E":
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

            case _:
                clogger.error(f"Unknown action: {action}")
                return




# Run the main method, possibly under a profiler, and highlight its errors
if __name__ == '__main__':
    try:
        if 'profile' in sys.argv:
            import cProfile
            cProfile.run('main()')
        else:
            main()
    except Exception:
        logging.error(traceback.format_exc())
    except SystemExit:
        clogger.error("Process interrupted")
    except KeyboardInterrupt:
        clogger.error("Process interrupted by user")
