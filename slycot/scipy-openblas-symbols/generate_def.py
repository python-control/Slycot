#!/usr/bin/env python

"""Generate scipy_openblas_symbols.def from 3 .txt files

Result is checked-in to repo; only re-run if symbols change.

The source data is from Scipy - see headers of the source .txt files
and comment in lapack_exclusions.py for detail.

Symbols generated work with scipy_openblas32 0.3.31.22.1 wheels, from
Feb 2026, and have only been tested with Slycot; there's no guarantee
they'll work with other BLAS and LAPACK dependent Fortran sources.

"""

from lapack_exclusions import lapack_exclusions

# scipy doesn't know about this one, but SLICOT needs it; from netlib:
#   DLAMC3 is intended to force  A  and  B  to be stored prior to doing
#   the addition of  A  and  B,  for use in situations where optimizers
#   might hold one of these in a register.
own_symbols = ["dlamc3"]


def parse_line(line):
    """process line; return symbol as string, or None if none found"""

    # lines are "rettype identifier(args...)"
    # we want "identifier"

    pline = line.strip()
    if not pline or pline.startswith("#"):
        return None

    if pline[0] == " " or " " not in pline:
        raise ValueError(f"bad line {line}")

    _, rest = pline.split(" ", 1)

    if not rest or rest[0] == "(" or "(" not in rest:
        raise ValueError(f"bad line {line}")

    symbol, _ = rest.split("(", 1)

    return symbol


def get_symbols(filename):
    """Get all symbols from file"""
    with open(filename) as infile:
        lines = infile.readlines()

    symbols = [parse_line(line) for line in lines]

    return [s for s in symbols if s is not None]


def find_intersections(data):
    """Check for intersection between the various symbol sources"""
    # this w
    from itertools import combinations

    for p1, p2 in combinations(data, 2):
        s1 = set(data[p1])
        if not len(s1) == len(data[p1]):
            raise ValueError(f"{p1} has duplicates: ")
        s2 = set(data[p2])
        if not len(s2) == len(data[p2]):
            raise ValueError(f"{p2} has duplicates")

        if s1.issubset(s2):
            raise ValueError(f"{p1} is a subset of {p2}")

        if s2.issubset(s1):
            raise ValueError(f"{p2} is a subset of {p1}")

        both = s1.intersection(s2)

        if both:
            raise ValueError(f"{p1} and {p2} have non-empty intersection: {both}")


def main():
    """Get symbols and create output file"""
    blas = get_symbols("cython_blas_signatures.txt")
    lapack = get_symbols("cython_lapack_signatures.txt")

    find_intersections(
        {
            "blas": blas,
            "lapack": lapack,
            "lapack_exclusions": lapack_exclusions,
            "own_symbols": own_symbols,
        }
    )

    with open("scipy-openblas-symbols.def", "wt") as outfile:
        for symbol in blas + lapack + own_symbols + lapack_exclusions:
            outfile.write(f"-D{symbol.upper()}=SCIPY_{symbol.upper()}\n")


if __name__ == "__main__":
    main()
