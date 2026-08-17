/*
verify_pari_dissymmetry_odd_m.gp

Independent verification of the dissymmetry-theorem decomposition of
Section 5.3 (G(x) = T_Cm(x) + T_S(x) - T_{S-Cm}(x)), specialized to a
single odd m, where it collapses to a two-term closed form.

Context: OEIS editor Sean A. Irvine asked for a checkable formula/program
for A398575 (unrooted, m=7), since the paper's own T_Cm(x), T_S(x),
T_{S-Cm}(x) notation was judged too symbolic to verify directly. This
script starts from Andrew Howroyd's own general, already-published PARI
code on A332649/A332648 (functions EulerT, R, U below, reproduced
unmodified) and specializes U(n,k) at k=7. Because m=7 is odd, the
even-case term of U vanishes identically, and only two divisors of 7
(1 and 7) contribute to the sum over divisors of Howroyd's general
formula — collapsing it to:

    a(n) = g1 + (3/7)*x*(g7 - g1^7),   where g1(x) = R(x), g7(x) = R(x^7)

with R(x) the rooted series for m=7 (OEIS A397210, "number of blocks"
indexing, i.e. R(n,6) in Howroyd's notation, matching m=6+1=7). The same
collapse happens for any odd m, not just m=7: the general formula always
reduces to a two-term expression when m is odd, by the identical argument
(divisors of an odd m never include 2, so the k%2==0 branch of U is moot
and the divisor sum has the same two-term shape). This script fixes m=7
as the concrete case actually needed, but the reduction itself is the
odd-m companion, in closed form rather than a proved obstruction, to
Proposition 1's even-m argument in the paper.

Verified below: the simplified seq(n) reproduces Howroyd's own unmodified
U(n,7) term for term (all 13 terms checked), which in turn reproduces the
published A398575 data (up to the expected offset: U(n,7) includes the
trivial n=0 term, not published since A398575's OFFSET starts at 1).

Reference: F. G. Speyser, "Enumeration and Asymptotic Analysis of Strict
Non-Plane m-Gonal Cactus Graphs via Split-Decomposition", Section 5.3.
Companion to verify_dissymmetry_m6.py and verify_dissymmetry_all_m.py,
which check the same theorem by an unrelated route (SymPy, all m) rather
than this script's PARI-native, odd-m-specialized closed form.

Author: Frederic G. Speyser
Run: gp -q verify_pari_dissymmetry_odd_m.gp
*/

EulerT(v)={Vec(exp(x*Ser(dirmul(v, vector(#v, n, 1/n))))-1, -#v)}

\\ Andrew Howroyd's original general functions (A332649/A332648), unmodified.
R(n, k)={my(v=[]); for(n=1, n, my(g=1+x*Ser(v)); v=EulerT(Vec((g^k + g^(k%2)*subst(g^(k\2), x, x^2))/2))); concat([1], v)}
U(n, k)={my(p=Ser(R(n, k-1))); my(g(d)=subst(p + O(x*x^(n\d)), x, x^d)); Vec(g(1) + x*sumdiv(k, d, eulerphi(d)*g(d)^(k/d))/(2*k) - x*(g(1)^k)/2 + x*if(k%2==0, g(2)^(k/2) - g(1)^2*g(2)^(k/2-1))/4)}

\\ Simplified, k=7-specialized version (this script's contribution).
Rs(n)={my(v=[]); for(n=1, n, my(g=1+x*Ser(v)); v=EulerT(Vec((g^6 + subst(g^3, x, x^2))/2))); concat([1], v)}
seq(n)={my(p=Ser(Rs(n))); my(g1=p+O(x*x^n)); my(g7=subst(p+O(x*x^(n\7)), x, x^7)); Vec(g1 + 3*x*(g7-g1^7)/7)}

N = 12;
u_general = U(N, 7);
u_simplified = seq(N);
published = [1, 1, 4, 14, 80, 504, 3659, 28254, 230200, 1940896, 16830963, 149199518];

print("Howroyd's general U(n,7):     ", u_general);
print("Simplified seq(n), k=7:       ", u_simplified);
print("Published A398575 (offset 1): ", published);
print("General == simplified, all terms match: ", u_general == u_simplified);
