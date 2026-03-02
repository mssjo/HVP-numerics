
from mpmath import mpf

TAB = '\t'

with open("amflow/AMFlow_data.dat", 'r') as infile:
    indata = [line.strip().split() for line in infile]
    print(f"Read Mathematica output from {infile.name}")

err = 1e-12
def reldiff(val, ref):
    return str((val - ref)/ref) if ref != 0 else "nan"

for n in range(6):
    with open(f"plots/prec64/E{n+1}e.dat", 'r') as reffile, open(f"plots/prec64/E{n+1}a.dat", 'w') as outfile:
        print('\t'.join((
            "t",
            "Re",
            "MinRe",
            "MaxRe",
            "Im",
            "MinIm",
            "MaxIm",
            "NormRe",
            "NormMinRe",
            "NormMaxRe",
            "NormIm",
            "NormMinIm",
            "NormMaxIm",
            "Time"
            )),
            file=outfile)

        ref = iter(reffile)
        next(ref) # Skip header
        for line in indata:
            t = mpf(line[0])
            try:
                while True:
                    tokens = next(ref).strip().split()
                    if tokens and mpf(tokens[0]) == t:
                        # print(f"{t}:\t{tokens[1]}\t{tokens[6]}")
                        break
            except StopIteration:
                raise IOError(f"Reference file lacks entry for t={t}")

            ref_real = mpf(tokens[1])
            val_real = mpf(line[2*n+1]) * (-1 if n < 3 else 1)
            ref_imag = mpf(tokens[6])
            val_imag = mpf(line[2*n+2]) * (-1 if n < 3 else 1)

            print('\t'.join(str(x) for x in (
                t,
                val_real,
                val_real - err,
                val_real + err,
                val_imag,
                val_imag - err,
                val_imag + err,
                reldiff(val_real,       ref_real),
                reldiff(val_real - err, ref_real),
                reldiff(val_real + err, ref_real),
                reldiff(val_imag,       ref_imag),
                reldiff(val_imag - err, ref_imag),
                reldiff(val_imag + err, ref_imag),
                mpf(line[2*6+1] if n < 3 else line[2*6+2])/3,
                )),
                file=outfile)

    print(f"Read reference values from {reffile.name}")
    print(f"Wrote data to {outfile.name}")
