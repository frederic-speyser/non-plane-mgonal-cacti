"""
verify_lemma5.py

Independent numerical verification of the three transition values quoted
in the paper's discussion of Conjecture 1 / Lemma 5 (Phi_y^(m+1) evaluated
at (rho_m, tau_m), the *numerical* criterion actually used in the text --
not the exact-but-circular Lemma 5 criterion itself, which needs s_{m+1}
(rho_m) rather than tau_m = s_m(rho_m)).

Important methodological note (kept here because it was learned the hard
way while writing this script): substituting the TRUE value s_{m+1}(rho_m)
in this formula, instead of tau_m, gives a substantially different number.
The two are genuinely different quantities -- see the paper's own discussion
of why tau_m is used as an approximation, with no established error bound.
This script reproduces the tau_m-based numbers quoted in the text, not the
(currently unusable) exact Lemma 5 criterion.

Author: Frederic G. Speyser
Run with: python3 verify_lemma5.py
"""
import numpy as np


def solve_s_fast(m, N):
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


def eval_series(coeffs, x):
    total, xp = 0.0, 1.0
    for c in coeffs:
        total += c * xp
        xp *= x
    return total


def phi_y_transition(m_next, rho_prev, tau_prev, N=400):
    """Phi_y^(m_next) at (rho_prev, tau_prev), with tau_prev substituted for
    the i=1 term's free variable, and the true s_{m_next} series used for
    the i>=2 tail (which is evaluated at safely small points rho_prev^i)."""
    s_next = solve_s_fast(m_next, N)
    odd = (m_next % 2 == 1)

    def KC_true_at(x):
        sx = eval_series(s_next, x)
        sx2 = eval_series(s_next, x ** 2)
        if odd:
            return 0.5 * (sx ** (m_next - 1) + sx2 ** ((m_next - 1) // 2))
        return 0.5 * (sx ** (m_next - 1) + sx * sx2 ** ((m_next - 2) // 2))

    tail = 0.0
    i = 2
    while rho_prev ** i > 1e-16:
        tail += KC_true_at(rho_prev ** i) / i
        i += 1

    sq = eval_series(s_next, rho_prev ** 2)
    if odd:
        K1 = 0.5 * (tau_prev ** (m_next - 1) + sq ** ((m_next - 1) // 2))
        dK1 = 0.5 * (m_next - 1) * tau_prev ** (m_next - 2)
    else:
        K1 = 0.5 * (tau_prev ** (m_next - 1) + tau_prev * sq ** ((m_next - 2) // 2))
        dK1 = 0.5 * ((m_next - 1) * tau_prev ** (m_next - 2) + sq ** ((m_next - 2) // 2))

    Phi = rho_prev * np.exp(K1 + tail)
    return Phi * dK1


CASES = [
    (6, 0.604765, 0.840896, 1.045),   # transition m=5 -> 6
    (7, 0.633235, 0.821008, 0.855),   # transition m=6 -> 7
    (8, 0.669930, 0.832683, 0.971),   # transition m=7 -> 8
]

if __name__ == "__main__":
    print("Transition   computed    reported (text)   |diff|")
    for m_next, rho, tau, reported in CASES:
        val = phi_y_transition(m_next, rho, tau)
        print(f"  -> m={m_next}     {val:.4f}      {reported}             {abs(val-reported):.4f}")
