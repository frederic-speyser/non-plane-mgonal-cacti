\\ verify_pari_euler.gp
\\
\\ Independent PARI/GP verification of the rooted series for m=5,6,7,8,
\\ via the Euler-transform recurrence (the same principle as Andrew
\\ Howroyd's code on A398033, generalized here to all four values of m
\\ treated in the paper, both parities).
\\
\\ For m odd:  K_C(x) = (1/2)*(s(x)^(m-1) + s(x^2)^((m-1)/2))
\\ For m even: K_C(x) = (1/2)*(s(x)^(m-1) + s(x)*s(x^2)^((m-2)/2))
\\
\\ This is a different code path from verify_pari.gp (which uses PARI's
\\ native truncated power series arithmetic directly on the functional
\\ equation) -- here the Euler transform explicitly builds the MSET
\\ construction of Section 5.2 of the paper, term by term.
\\
\\ Author: Frederic G. Speyser
\\ Run with: gp -q verify_pari_euler.gp

EulerT(v)={Vec(exp(x*Ser(dirmul(v, vector(#v, n, 1/n))))-1, -#v)}

seq_odd(m, n)={my(v=[]); for(k=1, n, my(g=1+x*Ser(v)); v=EulerT(Vec((g^(m-1) + subst(g^((m-1)/2), x, x^2))/2))); concat([1], v)}

seq_even(m, n)={my(v=[]); for(k=1, n, my(g=1+x*Ser(v)); v=EulerT(Vec((g^(m-1) + g*subst(g^((m-2)/2), x, x^2))/2))); concat([1], v)}

seq(m, n) = if(m%2==1, seq_odd(m,n), seq_even(m,n))

known5 = [1,1,3,13,62,333,1894,11258,68990,432964];
known6 = [1,1,4,22,140,985,7374,57577,463670,3822418];
known7 = [1,1,4,25,176,1397,11757,103376,937179,8699140];
known8 = [1,1,5,37,319,3059,31195,331991,3643790,40943462];

check(m, known) = {
  my(computed = seq(m, #known-1));
  if(computed == known,
     print("m=", m, ": OK, exact match on ", #known, " terms"),
     print("m=", m, ": MISMATCH -- computed=", computed, "  known=", known)
  );
}

check(5, known5);
check(6, known6);
check(7, known7);
check(8, known8);
