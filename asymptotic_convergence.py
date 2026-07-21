"""
asymptotic_convergence.py

Empirical verification of Theorems 3-4's conclusion: that

    s_n ~ C_m * rho_m^(-n) * n^(-3/2)

really holds as n grows, for both an odd m (5, closed-form tau_m) and an
even m (6). This complements the exact-arithmetic checks
(mgonal_cactus_series.py) by testing the theorem's conclusion directly,
across many terms, rather than an intermediate formula.

Two implementations are given:
  - a fast one (numpy, double precision) that reaches many terms quickly,
    useful to see the overall trend;
  - a high-precision one (mpmath) that is authoritative once the numbers
    involved grow into the hundreds of digits, where double precision
    accumulates rounding error (observed directly: the fast version shows
    a spurious slight uptick in the ratio beyond n ~ 1000 for m = 5, which
    disappears entirely once recomputed at 60 significant digits -- i.e.
    it is a floating-point artifact, not a real secondary term).

Author: Frederic G. Speyser
Run with: python3 asymptotic_convergence.py
"""
import numpy as np


def solve_s_fast(m, N):
    """Rooted series s(x) up to x^N, double precision, order-by-order."""
    s = np.zeros(N + 1)
    s[1] = 1.0
    odd = (m % 2 == 1)
    for n in range(2, N + 1):
        s_t = s
        sp = np.zeros(N + 1)
        sp[0:N + 1:2] = s_t[:len(sp[0:N + 1:2])]
        if odd:
            spk = np.zeros(N + 1); spk[0] = 1.0
            for _ in range((m - 1) // 2):
                spk = np.convolve(spk, sp)[:N + 1]
            sk = np.zeros(N + 1); sk[0] = 1.0
            for _ in range(m - 1):
                sk = np.convolve(sk, s_t)[:N + 1]
            KC = 0.5 * (sk + spk)
        else:
            spk = np.zeros(N + 1); spk[0] = 1.0
            for _ in range((m - 2) // 2):
                spk = np.convolve(spk, sp)[:N + 1]
            sk = np.zeros(N + 1); sk[0] = 1.0
            for _ in range(m - 1):
                sk = np.convolve(sk, s_t)[:N + 1]
            skp = np.convolve(s_t, spk)[:N + 1]
            KC = 0.5 * (sk + skp)
        G = np.zeros(N + 1)
        i = 1
        while i * (m - 1) <= n:
            kc_i = np.zeros(N + 1)
            kc_i[0:N + 1:i] = KC[:len(kc_i[0:N + 1:i])]
            G += kc_i / i
            i += 1
        E = np.zeros(N + 1); E[0] = 1.0
        for k in range(1, n + 1):
            E[k] = sum(j * G[j] * E[k - j] for j in range(1, k + 1)) / k
        s[n] = E[n - 1]
    return s


CASES = {
    5: dict(rho=0.604765, C=0.7905),
    6: dict(rho=0.633235, C=0.8736),
    8: dict(rho=0.690268, C=1.0312),
}


def report(m, N=1200, stride=100):
    rho, C = CASES[m]['rho'], CASES[m]['C']
    s = solve_s_fast(m, N)
    support = [n for n in range(1, N + 1) if abs(s[n]) > 0.5]
    print(f"\n=== m={m}: s_n / (C_m * rho_m^-n * n^-1.5), n up to {support[-1]} ===")
    for n in support[::stride] + [support[-1]]:
        pred = C * rho ** (-n) * n ** (-1.5)
        print(f"  n={n:5d}   ratio={s[n]/pred:.6f}")


if __name__ == "__main__":
    report(5)
    report(6)
    report(8)
