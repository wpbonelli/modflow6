"""Compare first_step.f90 output with the Python evaluation of the same forms.

Reads the Fortran program's output (three values per line: integer exponent,
real exponent, expm1/log1p) from stdin. For each case, prints the error in ulps
against exact rational arithmetic and whether the Fortran value is bit-identical
to the Python one.

Run with: ./first_step | python3 -I first_step_cmp.py
The case list must match first_step.f90.
"""

import math
import sys
from fractions import Fraction as F

CASES = [
    (1.0, 10, 1.0001),
    (1.0, 10, 1.000001),
    (1.0, 3, 1 + 1e-9),
    (100.0, 50, 0.999),
    (10.0, 100, 1.005),
]


def powi(r, n):
    res, b, m = 1.0, r, n
    while m:
        if m & 1:
            res *= b
        m >>= 1
        if m:
            b *= b
    return res


lines = [line for line in sys.stdin.read().splitlines() if line.strip()]
if len(lines) != len(CASES):
    sys.exit(f"expected {len(CASES)} lines of Fortran output, got {len(lines)}")
for (p, n, r), line in zip(CASES, lines):
    fortran = list(map(float, line.split()))
    python = [
        p * (1 - r) / (1 - powi(r, n)),
        p * (1 - r) / (1 - r ** float(n)),
        p * (r - 1) / math.expm1(n * math.log1p(r - 1)),
    ]
    ex = F(p) * (1 - F(r)) / (1 - F(r) ** n)
    u = F(math.ulp(float(ex)))
    cols = []
    for name, fv, pv in zip(["int", "real", "expm1"], fortran, python):
        err = float(abs(F(fv) - ex) / u)
        cols.append(f"{name} {err:10.1f} ulp ({'same' if fv == pv else 'DIFF'})")
    print(f"{p:6g} {n:4d} {r!r:>14}: " + ", ".join(cols))
