\\ verify_pari.gp
\\
\\ Independent verification, in PARI/GP, of the enumerative series (rooted
\\ case, m = 5, 6, 7, 8) and of Theorem 2's closed form for tau_m (m odd).
\\ This is a fourth independent implementation, alongside:
\\   - mgonal_cactus_series.py (Python, exact Fraction arithmetic)
\\   - a SymPy coefficient-balancing recomputation (see the paper's working
\\     notes)
\\   - Andrew Howroyd's own PARI code on OEIS A332648/A332649 (Euler
\\     transform / Dirichlet convolution approach)
\\ Here, PARI's native truncated power series arithmetic (type Ser) is used
\\ directly -- a different code path from all of the above.
\\
\\ Author: Frederic G. Speyser
\\ Run with: gp -q verify_pari.gp

N = 30;

\\ ---------- Part 1: rooted series s(x), all four values of m ----------
compute_series(m) = {
  my(s, KC, Gsum, s_new);
  s = x + O(x^(N+1));
  for(iter=1, 40,
    if(Mod(m,2)==1,
      KC = (1/2)*(s^(m-1) + subst(s,x,x^2)^((m-1)/2));
    ,
      KC = (1/2)*(s^(m-1) + s*subst(s,x,x^2)^((m-2)/2));
    );
    Gsum = 0;
    for(i=1, N, Gsum = Gsum + subst(KC,x,x^i)/i);
    s_new = x + x*(exp(Gsum) - 1);
    s = s_new + O(x^(N+1));
  );
  return(s);
}

print("=== Rooted series s(x), m = 5, 6, 7, 8 (PARI native Ser arithmetic) ===");
for(mi=0,3, m=5+mi; s=compute_series(m); v=Vec(s); print("m=",m,"  ",v));

\\ ---------- Part 2: Theorem 2, m odd -- tau_m closed form ----------
print();
print("=== Theorem 2 (m odd): tau_m, numerical root-finding vs closed form ===");
for(k=0,3,m=5+2*k; tn=solve(t=0.5,1.0,(1/2)*(m-1)*t^(m-1)-1); tc=(2/(m-1))^(1/(m-1)); print("m=",m,"  tau(PARI numeric root)=",tn,"  tau(closed form)=",tc,"  diff=",abs(tn-tc)))
