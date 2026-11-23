"""Middle-step and last-step error within one geometric TDIS period.

Replays tdis_set_delt / tdis_set_timestep for a single period starting at t = 0:

    kstp == 1      delt = first-step formula
    1 < kstp < n   delt = tsmult * delt
    kstp == n      delt = topertim + perlen - totim
    every step     totim = totim + delt

and compares delt and totim at each step with the exact geometric values for the
same double inputs. Run once with the current integer-exponent first step and once
with an accurate (real-exponent) first step, to separate first-step error from
compounding and last-step cancellation.

Produces the last-step table and middle-step figures in analysis.md, section
"Current behavior". Python only: not checked against a Fortran build.
Run with: python3 -I step_errors.py
"""

import math
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


def ulps(x, ex):
    return float(abs(F(x) - ex) / F(math.ulp(float(ex))))


def run(perlen, n, r, first):
    d1 = F(perlen) * (1 - F(r)) / (1 - F(r) ** n)  # exact first step
    topertim = 0.0
    t = 0.0
    max_mid = max_t = 0.0
    for k in range(1, n + 1):
        if k == 1:
            d = first(perlen, n, r)
        elif k == n:
            d = topertim + perlen - t
        else:
            d = r * d
        t = t + d
        if k < n:
            max_mid = max(max_mid, ulps(d, d1 * F(r) ** (k - 1)))
            max_t = max(max_t, ulps(t, d1 * (F(r) ** k - 1) / (F(r) - 1)))
    d_last = d1 * F(r) ** (n - 1)
    rel = float(abs(F(d) - d_last) / d_last)
    return max_mid, max_t, ulps(d, d_last), rel


def first_int(p, n, r):
    return p * (1 - r) / (1 - powi(r, n))


def first_real(p, n, r):
    return p * (1 - r) / (1 - r ** float(n))


print(
    f"{'perlen':>7} {'nstp':>5} {'tsmult':>7} {'first step':>11}"
    f" {'max mid delt':>13} {'max totim':>10} {'last delt':>10} {'last rel':>9}"
)
for perlen, n, r in [(10.0, 100, 1.005), (10.0, 100, 0.95), (1000.0, 200, 0.97)]:
    for name, first in [("integer", first_int), ("real", first_real)]:
        mid, tt, last, rel = run(perlen, n, r, first)
        print(
            f"{perlen:7g} {n:5d} {r:7g} {name:>11} {mid:10.0f} ulp {tt:6.0f} ulp"
            f" {last:6.0f} ulp {rel:9.1e}"
        )
