"""
verify_dissymmetry_all_m.py

Independent verification of the unrooted series for all four values of
m treated in the paper (5, 6, 7, 8) by assembling the dissymmetry theorem
(Section 5.3) entirely from scratch in SymPy: G(x) = T_Cm(x) + T_S(x) -
T_{S-Cm}(x), including an explicit symbolic implementation of the dihedral
cycle index Z_Dm (equation 8, both parity branches), using only the
already-independently-verified rooted series as input.

Generalizes verify_dissymmetry_m6.py (m=6 only) to handle both parities.
A genuinely different computational path from every other script in this
repository: no Euler transform, no direct graph construction, just direct
symbolic assembly of the three dissymmetry-theorem components.

Author: Frederic G. Speyser
Run with: python3 verify_dissymmetry_all_m.py   (requires: pip install sympy)
"""
import sympy as sp

x = sp.symbols('x')


def trunc(expr, N):
    return sp.series(expr, x, 0, N + 1).removeO()


def verify_unrooted(m, a_rooted, N):
    s_expr = sum(a_rooted[k] * x**(1 + (m - 1) * k)
                 for k in range(len(a_rooted)) if 1 + (m - 1) * k <= N)
    s2 = s_expr.subs(x, x**2)

    if m % 2 == 1:
        KC = trunc(sp.Rational(1, 2) * (s_expr**(m - 1) + s2**((m - 1) // 2)), N)
    else:
        KC = trunc(sp.Rational(1, 2) * (s_expr**(m - 1) + s_expr * s2**((m - 2) // 2)), N)

    Gexp = 0
    i = 1
    while i * (m - 1) <= N:
        Gexp += KC.subs(x, x**i) / i
        i += 1
    Gexp = trunc(Gexp, N)

    E = trunc(sp.exp(Gexp), N)
    SX = trunc(x * (E - 1), N)
    SC = trunc(E - 1 - KC, N)
    TS = trunc(x * SC, N)
    TSCm = trunc(KC * SX, N)

    Sd = [s_expr.subs(x, x**d) for d in range(1, m + 1)]
    TCm = 0
    for d in sp.divisors(m):
        TCm += sp.totient(d) * Sd[d - 1]**(m // d)
    TCm = TCm / (2 * m)
    if m % 2 == 1:
        TCm += Sd[0] * Sd[1]**((m - 1) // 2) / 2
    else:
        TCm += (Sd[0]**2 * Sd[1]**((m - 2) // 2) + Sd[1]**(m // 2)) / 4
    TCm = trunc(TCm, N)

    G = sp.expand(trunc(TCm + TS - TSCm, N))
    poly = sp.Poly(G, x)
    return sorted([(int(d[0]), int(c)) for d, c in zip(poly.monoms(), poly.coeffs())])


# Rooted series a(k), already independently verified (PARI Euler-transform
# recurrence, PARI native series, Python exact-Fraction recursion).
ROOTED = {
    5: [1, 1, 3, 13, 62, 333, 1894],
    6: [1, 1, 4, 22, 140, 985, 7374],
    7: [1, 1, 4, 25, 176, 1397, 11757],
    8: [1, 1, 5, 37, 319, 3059, 31195],
}

# Known unrooted targets, for comparison.
KNOWN_UNROOTED = {
    5: [1, 1, 3, 8, 31],
    6: [1, 1, 4, 13, 67],
    7: [1, 1, 4, 14, 80],
    8: [1, 1, 5, 20, 143],
}

N_BY_M = {5: 22, 6: 27, 7: 32, 8: 37}

if __name__ == "__main__":
    for m in [5, 6, 7, 8]:
        result = verify_unrooted(m, ROOTED[m], N_BY_M[m])
        computed = [c for _, c in result]
        known = KNOWN_UNROOTED[m]
        ok = computed[:len(known)] == known
        print(f"m={m}: {'OK' if ok else 'MISMATCH'}, computed={computed[:len(known)]}, known={known}")
