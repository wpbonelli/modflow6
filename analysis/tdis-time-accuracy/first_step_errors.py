"""Error in the TDIS first time step length for a geometric period.

Evaluates the first-step formula from tdis_set_delt / check_tdis_timing

    delt = perlen * (1 - tsmult) / (1 - tsmult**nstp)

in double precision with several ways of computing the power, and compares each
with the exact value of the expression for the same double inputs (computed with
fractions.Fraction). Errors are in ulps of the exact result.

Forms:
  int    integer exponent by repeated squaring, like libgcc's __powidf2 (gfortran)
  real   real exponent, tsmult**real(nstp, DP), i.e. libm pow (#3204)
  expm1  perlen * (tsmult - 1) / expm1(nstp * log1p(tsmult - 1))
  hybrid expm1 form for |nstp * log1p(tsmult - 1)| < 1, real form otherwise

Produces the tables in analysis.md, section "First time step length".
Run with: python3 -I first_step_errors.py
"""

import math
import random
from fractions import Fraction as F


def powi(r, n):
    """r**n by repeated squaring, rounding at each multiply (__powidf2)."""
    res, b, m = 1.0, r, n
    while m:
        if m & 1:
            res *= b
        m >>= 1
        if m:
            b *= b
    return res


def exact(p, n, r):
    return F(p) * (1 - F(r)) / (1 - F(r) ** n)


def forms(p, n, r):
    y = n * math.log1p(r - 1)
    f = {
        "int": p * (1 - r) / (1 - powi(r, n)),
        "real": p * (1 - r) / (1 - r ** float(n)),
        "expm1": p * (r - 1) / math.expm1(y),
    }
    f["hybrid"] = f["expm1"] if abs(y) < 1 else f["real"]
    return y, f


def ulps(x, ex):
    return float(abs(F(x) - ex) / F(math.ulp(float(ex))))


def table(title, cases):
    print(f"\n{title}")
    print(
        f"{'perlen':>8} {'nstp':>5} {'tsmult':>12} {'y':>9}"
        f" {'int':>10} {'real':>10} {'expm1':>8} {'hybrid':>8}  (ulp)"
    )
    for p, n, r in cases:
        ex = exact(p, n, r)
        y, f = forms(p, n, r)
        e = {k: ulps(v, ex) for k, v in f.items()}
        print(
            f"{p:8g} {n:5d} {r!r:>12} {y:9.3g} {e['int']:10.1f} {e['real']:10.1f}"
            f" {e['expm1']:8.1f} {e['hybrid']:8.1f}"
        )


table(
    "Near 1 (cancellation dominates)",
    [
        (1.0, 10, 1.01),
        (1.0, 10, 1.001),
        (1.0, 5, 1.001),
        (1.0, 10, 1.0001),
        (1.0, 10, 1.000001),
        (1.0, 3, 1 + 1e-9),
        (100.0, 50, 0.999),
    ],
)
table(
    "Large |y| (expm1/log1p alone loses accuracy)",
    [
        (1.0, 100, 1.005),
        (1.0, 100, 1.05),
        (1.0, 100, 1.2),
        (1.0, 200, 1.2),
        (1.0, 1000, 1.5),
        (1.0, 100, 0.95),
        (1.0, 500, 0.9),
    ],
)
table(
    "Cases from #3204",
    [
        (10.0, 100, 1.005),
        (1095.75, 100, 1.05),
        (5000.0, 100, 1.02),
        (365.25, 6, 1.3),
    ],
)

# Random survey. The seed and the order of random calls reproduce the numbers in
# analysis.md; don't reorder them.
random.seed(1)
cases = []
for _ in range(20000):
    cases.append(
        (
            random.choice([1.0, 10.0, 365.25, 1e4]),
            random.randint(2, 1000),
            1 + random.choice([-1, 1]) * 10 ** random.uniform(-8, -0.3),
        )
    )
worst = {"real": 0.0, "expm1": 0.0, "hybrid": 0.0}
for p, n, r in cases:
    ex = exact(p, n, r)
    _, f = forms(p, n, r)
    for k in worst:
        worst[k] = max(worst[k], ulps(f[k], ex))
print("\nWorst over 20000 random cases (|tsmult - 1| in 1e-8..0.5, nstp 2..1000):")
for k, v in worst.items():
    print(f"  {k:>6}: {v:.1f} ulp")
