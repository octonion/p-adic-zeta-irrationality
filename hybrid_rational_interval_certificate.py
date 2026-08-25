#!/usr/bin/env python3
"""Rigorous rational fixed-point certificate for the 22 hybrid margins.

The program uses only Python's standard-library arbitrary-precision integers and
``fractions.Fraction``.  Every real interval is represented by two integers
``(lo, hi)`` meaning [lo/SCALE, hi/SCALE], with all arithmetic rounded outwards.

Transcendental enclosures are obtained from elementary identities with explicit
remainder bounds:

* pi = 16 atan(1/5) - 4 atan(1/239) (Machin);
* log(p) = 2 atanh((p-1)/(p+1));
* atan(z) uses its alternating power series, with the next term as remainder;
* atanh(z) uses its positive power series and a geometric tail bound;
* sqrt(a/b) is enclosed by integer square root at the fixed scale;
* acos(x) = 2 atan(sqrt((1-x)/(1+x))) for 0 < x < 1.

Thus no binary floating point, mpmath, MPFI, Arb, or external numerical library
is used in the certification step.
"""

from __future__ import annotations

from fractions import Fraction
from math import isqrt
from typing import Dict, Iterable, Tuple

PRECISION_DIGITS = 90
PRINT_DIGITS = 60
SCALE = 10**PRECISION_DIGITS
# The power-series loops stop once the rigorous next-term/tail enclosure is at
# most this many units of 1/SCALE.  The final interval widths are still far
# below 10^(-70), while the smallest margin is > 10^(-1).
SERIES_STOP_UNITS = 1000

Interval = Tuple[int, int]


def floor_div(a: int, b: int) -> int:
    assert b > 0
    return a // b


def ceil_div(a: int, b: int) -> int:
    assert b > 0
    return -((-a) // b)


def interval_fraction(q: Fraction) -> Interval:
    return (
        floor_div(SCALE * q.numerator, q.denominator),
        ceil_div(SCALE * q.numerator, q.denominator),
    )


def interval_add(x: Interval, y: Interval) -> Interval:
    return x[0] + y[0], x[1] + y[1]


def interval_sub(x: Interval, y: Interval) -> Interval:
    return x[0] - y[1], x[1] - y[0]


def interval_mul(x: Interval, y: Interval) -> Interval:
    products = (
        x[0] * y[0],
        x[0] * y[1],
        x[1] * y[0],
        x[1] * y[1],
    )
    return floor_div(min(products), SCALE), ceil_div(max(products), SCALE)


def interval_mul_int(x: Interval, n: int) -> Interval:
    if n >= 0:
        return x[0] * n, x[1] * n
    return x[1] * n, x[0] * n


def interval_mul_fraction(x: Interval, q: Fraction) -> Interval:
    n, d = q.numerator, q.denominator
    if n >= 0:
        return floor_div(x[0] * n, d), ceil_div(x[1] * n, d)
    return floor_div(x[1] * n, d), ceil_div(x[0] * n, d)


def positive_interval_mul_rational(x: Interval, n: int, d: int) -> Interval:
    """Multiply a nonnegative interval by the nonnegative rational n/d."""
    assert 0 <= x[0] <= x[1]
    assert n >= 0 and d > 0
    return floor_div(x[0] * n, d), ceil_div(x[1] * n, d)


def sqrt_fraction_bounds(q: Fraction) -> Interval:
    """Rigorous SCALE-fixed-point enclosure of sqrt(q), q >= 0."""
    assert q >= 0
    scaled_floor = (q.numerator * SCALE * SCALE) // q.denominator
    root = isqrt(scaled_floor)
    # The isqrt of floor(q*SCALE^2) already gives the lower endpoint; these
    # checks make the intended exact inequalities explicit.
    while (root + 1) ** 2 * q.denominator <= q.numerator * SCALE * SCALE:
        root += 1
    while root**2 * q.denominator > q.numerator * SCALE * SCALE:
        root -= 1
    if root**2 * q.denominator == q.numerator * SCALE * SCALE:
        return root, root
    return root, root + 1


def atan_rational_bounds(num: int, den: int) -> Interval:
    """Rigorous enclosure of atan(num/den) for 0 <= num/den < 1."""
    assert 0 <= num < den
    if num == 0:
        return 0, 0

    term = (
        floor_div(SCALE * num, den),
        ceil_div(SCALE * num, den),
    )
    partial: Interval = (0, 0)
    k = 0
    plus = True

    while True:
        partial = interval_add(partial, term) if plus else interval_sub(partial, term)

        # t_{k+1} = t_k * z^2 * (2k+1)/(2k+3), where
        # t_k = z^(2k+1)/(2k+1).
        factor_num = num * num * (2 * k + 1)
        factor_den = den * den * (2 * k + 3)
        next_term = positive_interval_mul_rational(term, factor_num, factor_den)

        if next_term[1] <= SERIES_STOP_UNITS:
            # Alternating-series remainder has the sign of the next term and
            # absolute value at most that term.
            if plus:  # the next term is negative
                return partial[0] - next_term[1], partial[1]
            return partial[0], partial[1] + next_term[1]

        term = next_term
        plus = not plus
        k += 1
        if k > 1_000_000:
            raise RuntimeError("atan series failed to converge")


def log_integer_bounds(p: int) -> Interval:
    """Rigorous enclosure of log(p), p >= 2."""
    assert p >= 2
    num, den = p - 1, p + 1
    term = (
        floor_div(SCALE * num, den),
        ceil_div(SCALE * num, den),
    )
    partial: Interval = (0, 0)
    k = 0

    while True:
        partial = interval_add(partial, term)
        factor_num = num * num * (2 * k + 1)
        factor_den = den * den * (2 * k + 3)
        next_term = positive_interval_mul_rational(term, factor_num, factor_den)

        # The remaining positive terms have successive ratio < z^2.
        tail_upper = ceil_div(
            next_term[1] * den * den,
            den * den - num * num,
        )
        if tail_upper <= SERIES_STOP_UNITS:
            return 2 * partial[0], 2 * (partial[1] + tail_upper)

        term = next_term
        k += 1
        if k > 1_000_000:
            raise RuntimeError("atanh series failed to converge")


def pi_bounds() -> Interval:
    a = atan_rational_bounds(1, 5)
    b = atan_rational_bounds(1, 239)
    return interval_sub(interval_mul_int(a, 16), interval_mul_int(b, 4))


def acos_fraction_bounds(x: Fraction) -> Interval:
    assert 0 < x < 1
    t = sqrt_fraction_bounds((1 - x) / (1 + x))
    lower_atan = atan_rational_bounds(t[0], SCALE)
    upper_atan = atan_rational_bounds(t[1], SCALE)
    return 2 * lower_atan[0], 2 * upper_atan[1]


def euler_phi(n: int) -> int:
    result = n
    x = n
    q = 2
    while q * q <= x:
        if x % q == 0:
            while x % q == 0:
                x //= q
            result -= result // q
        q += 1
    if x > 1:
        result -= result // x
    return result


def harmonic_exact(k: int) -> Fraction:
    return sum((Fraction(1, j) for j in range(1, k + 1)), Fraction(0))


def xi_witness(p: int, s: int) -> Fraction:
    return Fraction(12 * (s * (s + 1) - 1), 12 * s * s + (s - 1) * (p + 1))


def prime_window_exact(s: int, xi: Fraction) -> Fraction:
    m = s + 1
    k = Fraction(s, 1) // xi
    delta = Fraction(m, 1) - Fraction(k + 1, 1) * xi
    positive = max(Fraction(0), delta)
    return (
        Fraction(2 * m - 1, 2) * harmonic_exact(k)
        - k * xi
        + positive * positive / Fraction(2 * (k + 1), 1)
    )


def small_prime_cost_exact(p: int, s: int, xi: Fraction) -> Fraction:
    return s * s * xi - (s - 1) * (xi - Fraction(p + 1, 24) * xi * xi)


def tau_exact(p: int, s: int) -> Tuple[Fraction, Fraction, Fraction, Fraction]:
    xi = xi_witness(p, s)
    window = prime_window_exact(s, xi)
    small = small_prime_cost_exact(p, s, xi)
    tau = Fraction(2, (s + 1) ** 2) * (small + s * window)
    return tau, xi, window, small


WITNESSES: Dict[Tuple[int, int], Fraction] = {
    (2, 3): Fraction(16, 79),
    (2, 5): Fraction(6, 49),
    (2, 7): Fraction(1, 11),
    (2, 9): Fraction(3, 43),
    (2, 11): Fraction(3, 52),
    (2, 13): Fraction(2, 41),
    (2, 15): Fraction(3, 71),
    (2, 17): Fraction(3, 80),
    (2, 19): Fraction(1, 30),
    (2, 21): Fraction(1, 33),
    (2, 23): Fraction(1, 36),
    (2, 25): Fraction(1, 39),
    (2, 27): Fraction(1, 42),
    (2, 29): Fraction(1, 46),
    (3, 3): Fraction(4, 27),
    (3, 5): Fraction(7, 73),
    (3, 7): Fraction(1, 15),
    (3, 9): Fraction(1, 19),
    (3, 11): Fraction(3, 70),
    (5, 3): Fraction(7, 71),
    (5, 5): Fraction(1, 16),
    (7, 3): Fraction(1, 14),
}

CASES = (
    [(2, s) for s in range(3, 30, 2)]
    + [(3, s) for s in range(3, 12, 2)]
    + [(5, 3), (5, 5), (7, 3)]
)


def decimal_lower(value: int, digits: int = PRINT_DIGITS) -> str:
    """Decimal rounded toward -infinity from a scaled integer endpoint."""
    assert 0 <= digits <= PRECISION_DIGITS
    factor = 10 ** (PRECISION_DIGITS - digits)
    q = floor_div(value, factor)
    sign = "-" if q < 0 else ""
    q = abs(q)
    return f"{sign}{q // 10**digits}.{q % 10**digits:0{digits}d}"


def decimal_upper(value: int, digits: int = PRINT_DIGITS) -> str:
    """Decimal rounded toward +infinity from a scaled integer endpoint."""
    assert 0 <= digits <= PRECISION_DIGITS
    factor = 10 ** (PRECISION_DIGITS - digits)
    q = ceil_div(value, factor)
    sign = "-" if q < 0 else ""
    q = abs(q)
    return f"{sign}{q // 10**digits}.{q % 10**digits:0{digits}d}"


def margin_interval(
    p: int,
    s: int,
    y: Fraction,
    pi_interval: Interval,
    log_intervals: Dict[int, Interval],
) -> Tuple[Interval, Fraction, Fraction, Iterable[int]]:
    tau, xi, _window, _small = tau_exact(p, s)

    assert 1 < xi < s + 1
    assert Fraction(p + 1, 12) * xi < 1
    assert 0 < y < Fraction(1, p)

    lambda_interval = interval_sub(
        interval_mul_fraction(log_intervals[p], Fraction(12, p - 1)),
        interval_mul_fraction(pi_interval, 2 * y),
    )

    collision: Interval = (0, 0)
    layers = []
    c = p
    while c * y < 1:
        layers.append(c)
        x = c * y
        acos_interval = acos_fraction_bounds(x)
        sqrt_interval = sqrt_fraction_bounds(1 - x * x)
        bracket = interval_sub(acos_interval, interval_mul_fraction(sqrt_interval, x))
        term = interval_mul_fraction(
            interval_mul(pi_interval, bracket),
            Fraction(4 * euler_phi(c), c * c),
        )
        collision = interval_add(collision, term)
        c += p

    margin = interval_sub(
        interval_sub(
            interval_mul_int(lambda_interval, s),
            interval_fraction(Fraction(s + 1, 1) * tau),
        ),
        collision,
    )
    return margin, tau, xi, layers


def main() -> None:
    assert len(CASES) == 22
    assert set(CASES) == set(WITNESSES)

    pi_interval = pi_bounds()
    log_intervals = {p: log_integer_bounds(p) for p in (2, 3, 5, 7)}

    print("RIGOROUS RATIONAL FIXED-POINT CERTIFICATE FOR 22 HYBRID MARGINS")
    print(f"scale = 10^{PRECISION_DIGITS}")
    print("all interval operations use outward-rounded arbitrary-precision integers")
    print("all algebraic quantities are exact Fractions")
    print()
    print("p  s   xi                 tau_hyb                       Y      layers       margin interval")
    print("-" * 186)

    records = []
    for p, s in CASES:
        y = WITNESSES[(p, s)]
        margin, tau, xi, layers = margin_interval(p, s, y, pi_interval, log_intervals)
        assert margin[0] > 0, (p, s, margin)
        lo = decimal_lower(margin[0])
        hi = decimal_upper(margin[1])
        records.append((margin[0], p, s, y, lo, hi))
        print(
            f"{p:1d} {s:2d}  {str(xi):18s} {str(tau):29s} {str(y):7s} "
            f"{str(list(layers)):14s} [{lo}, {hi}]"
        )

    smallest = min(records)
    print()
    print("CERTIFIED: True")
    print("Smallest certified lower endpoint occurs at (p,s)=", (smallest[1], smallest[2]))
    print("Y =", smallest[3])
    print("lower endpoint =", smallest[4])
    print("upper endpoint =", smallest[5])


if __name__ == "__main__":
    main()
