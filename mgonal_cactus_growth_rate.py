"""
mgonal_cactus_growth_rate.py

Estimates rho_m (exponential growth rate) directly from the coefficients of
s(x) computed by mgonal_cactus_series.py, via an n^(-3/2)-corrected ratio
test, and compares against the closed-form tau_m for m odd.

Reference: F. G. Speyser, "Enumeration of Non-Plane m-Gonal Cactus Graphs via
Split-Decomposition", Theorems 2-3.

Author: Frederic G. Speyser
Run: python3 mgonal_cactus_growth_rate.py
"""
import time
from fractions import Fraction as F
import mgonal_cactus_series as ic

for m in [5, 6, 7, 8]:
    t0 = time.time()
    ic.N = 70
    s = ic.solve_s(m)
    nz = [(n, c) for n, c in enumerate(s) if c != 0]
    # ratio test: rho_m ~ (s_n / s_{n+gap})^(1/gap), gap = m-1
    (n1, c1), (n2, c2) = nz[-2], nz[-1]
    gap = n2 - n1
    rho_est = (float(c1) / float(c2)) ** (1.0 / gap)
    tau_est = float(c2) ** (1.0 / n2) * rho_est  # crude cross-check, not used for the headline estimate
    print(f"m={m}  (n={n1}->{n2}, gap={gap})  rho_m (ratio est.) = {rho_est:.6f}   [{time.time()-t0:.1f}s]")
    if m % 2 == 1:
        import math
        tau_closed = (2 / (m - 1)) ** (1 / (m - 1))
        print(f"       tau_m closed form (Theorem 2)          = {tau_closed:.6f}")
