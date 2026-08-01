"""
verify_dissymmetry_m6.py

Independent verification of the unrooted series for m=6 (A398035) by
assembling the dissymmetry theorem (Section 5.3 of the paper) entirely
from scratch in SymPy: G(x) = T_Cm(x) + T_S(x) - T_{S-Cm}(x), computed
via explicit symbolic Taylor-coefficient balancing, using only the
already-independently-verified rooted series as input.

This is a genuinely different computational path from every other script
in this repository: no Euler transform (unlike mgonal_cactus_series.py
and verify_pari.gp), no direct graph construction (unlike exhaustive_iso.py),
just direct symbolic assembly of the three dissymmetry-theorem components,
including an explicit implementation of the dihedral cycle index Z_D6.

Author: Frederic G. Speyser
Run with: python3 verify_dissymmetry_m6.py   (requires: pip install sympy)
"""
import sympy as sp

x = sp.symbols('x')

# Rooted series a(k) for m=6, already independently verified by four
# other methods (mgonal_cactus_series.py, verify_pari.gp, the Howroyd-style
# PARI recurrence, and separately by SymPy coefficient-balancing).
A_ROOTED_M6 = [1, 1, 4, 22, 140, 985, 7374, 57577, 463670, 3822418,
               32097451, 273570649, 2360512647, 20579056156, 180993484480,
               1603956849975, 14308385056256]

# Known target: A398035 (unrooted, m=6), for comparison.
A398035_KNOWN = [1, 1, 4, 13, 67, 372, 2419, 16551, 119995, 898848]


def trunc(expr, N):
    return sp.series(expr, x, 0, N + 1).removeO()


def verify(m, a_rooted, N):
    s_expr = sum(a_rooted[k] * x**(1 + (m - 1) * k)
                 for k in range(len(a_rooted)) if 1 + (m - 1) * k <= N)

    s2 = s_expr.subs(x, x**2)
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

    # T_Cm(x) via the dihedral cycle index Z_Dm (equation 8 of the paper),
    # m even case.
    Sd = [s_expr.subs(x, x**d) for d in range(1, m + 1)]
    TCm = 0
    for d in sp.divisors(m):
        TCm += sp.totient(d) * Sd[d - 1]**(m // d)
    TCm = TCm / (2 * m)
    TCm += (Sd[0]**2 * Sd[1]**((m - 2) // 2) + Sd[1]**(m // 2)) / 4
    TCm = trunc(TCm, N)

    G = sp.expand(trunc(TCm + TS - TSCm, N))
    poly = sp.Poly(G, x)
    return sorted([(int(d[0]), int(c)) for d, c in zip(poly.monoms(), poly.coeffs())])


if __name__ == "__main__":
    N = 26  # keeps runtime under ~1 minute; see the paper's repo notes
             # on SymPy's practical depth limit for this construction
    result = verify(6, A_ROOTED_M6, N)
    print("G(x) unrooted, m=6, via full dissymmetry-theorem assembly in SymPy:")
    print(result)
    print()
    computed = [c for _, c in result]
    print("n     computed    A398035        match")
    for i, (known, comp) in enumerate(zip(A398035_KNOWN, computed), start=1):
        print(f"{i:3d}   {comp:10d}   {known:10d}      {'OK' if comp == known else 'MISMATCH'}")
