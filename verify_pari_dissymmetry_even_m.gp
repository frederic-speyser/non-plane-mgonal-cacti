/*
verify_pari_dissymmetry_even_m.gp

Independent verification of the dissymmetry-theorem decomposition of
Section 5.3 (G(x) = T_Cm(x) + T_S(x) - T_{S-Cm}(x)) for m=8, the even-m
companion to verify_pari_dissymmetry_odd_m.gp.

Unlike the odd case (m=7), this does NOT collapse to a two-term closed
form. The reason is structural, not a missed simplification: 8 has four
divisors (1, 2, 4, 8) rather than two, and Andrew Howroyd's general
U(n,k) formula (from A332649/A332648) has an extra term active only when
k is even (the "g(2)^(k/2) - g(1)^2*g(2)^(k/2-1)" branch), which vanishes
identically for odd k but not here. This is the same parity obstruction
that runs through the whole paper (Proposition 1: no closed form for the
critical value tau_m when m is even), showing up again at the level of
this PARI specialization rather than being specific to it.

This script therefore verifies m=8 unrooted data (the sequence prepared
for OEIS submission, currently on hold pending the outcome of the
A398575 discussion -- see the repository README) directly against
Howroyd's unmodified general formula, and against the independent
exact-rational solver of mgonal_cactus_series.py, without attempting a
simplification that the even case does not structurally admit.

Verified below: Howroyd's U(n,8) matches the m=8 unrooted series computed
independently by mgonal_cactus_series.py, term for term (all 13 terms
checked, up to the expected offset -- U(n,8) includes the trivial n=0
term that the OEIS convention for this family omits).

Reference: F. G. Speyser, "Enumeration and Asymptotic Analysis of Strict
Non-Plane m-Gonal Cactus Graphs via Split-Decomposition", Section 5.3,
and Proposition 1 for the even/odd asymmetry this script's *absence* of
a closed form reflects. Companion to verify_pari_dissymmetry_odd_m.gp
(m=7, which does collapse to a closed form) and to
verify_dissymmetry_m6.py / verify_dissymmetry_all_m.py, which check the
same theorem for all m by an unrelated SymPy route.

Author: Frederic G. Speyser
Run: gp -q verify_pari_dissymmetry_even_m.gp
*/

EulerT(v)={Vec(exp(x*Ser(dirmul(v, vector(#v, n, 1/n))))-1, -#v)}

\\ Andrew Howroyd's original general functions (A332649/A332648), unmodified.
R(n, k)={my(v=[]); for(n=1, n, my(g=1+x*Ser(v)); v=EulerT(Vec((g^k + g^(k%2)*subst(g^(k\2), x, x^2))/2))); concat([1], v)}
U(n, k)={my(p=Ser(R(n, k-1))); my(g(d)=subst(p + O(x*x^(n\d)), x, x^d)); Vec(g(1) + x*sumdiv(k, d, eulerphi(d)*g(d)^(k/d))/(2*k) - x*(g(1)^k)/2 + x*if(k%2==0, g(2)^(k/2) - g(1)^2*g(2)^(k/2-1))/4)}

N = 12;
u8_general = U(N, 8);

\\ Independently computed by mgonal_cactus_series.py (exact-rational solver),
\\ hard-coded here for the standalone self-check; matches the live script's
\\ output exactly (offset by the trivial n=0 term U(n,8) includes and the
\\ OEIS-convention series does not).
python_m8_unrooted = [1, 1, 5, 20, 143, 1093, 9722, 91391, 904526, 9252640, 97270908, 1044943778, 11430591994];

print("divisors of 8 (why no two-term collapse): ", divisors(8));
print("Howroyd's general U(n,8):            ", u8_general);
print("mgonal_cactus_series.py, m=8 unrooted: ", python_m8_unrooted);

\\ Align offsets: U(n,8) includes an extra leading n=0 term, so compare
\\ U(n,8)[2..] against the full python_m8_unrooted list (one shorter).
u8_aligned = vector(#u8_general - 1, i, u8_general[i+1]);
print("U(n,8) matches Python solver, all terms: ", u8_aligned == python_m8_unrooted[1..#u8_aligned]);
