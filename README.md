# Enumeration and Asymptotic Analysis of Strict Non-Plane m-Gonal Cactus Graphs - verification code

## About this research

A *cactus graph* is a connected graph in which every edge lies on at most
one cycle — equivalently, every block (maximal 2-connected piece) is either
a bridge or a cycle. Under the name *Husimi trees*, cacti were first studied
in statistical mechanics, and their enumeration has been a recurring problem
in analytic combinatorics ever since. This repository accompanies a paper
that enumerates *strict m-gonal cacti* — cacti in which every block is a
cycle of the same fixed length *m*, with no bridges at all — in the *free*
(non-plane) setting: the graph is counted as an abstract object, with no
cyclic order imposed on the blocks meeting at a shared vertex (as opposed
to the *plane*, or embedded, case, where such an order is fixed and
distinguishes otherwise-identical graphs).

Existing literature (Bahrani & Lumbroso, 2018; Bóna, Bousquet, Labelle &
Leroux, 2000) provides a general split-decomposition framework that covers
this free/non-plane case in principle, but never instantiates it
numerically for any fixed *m* ≥ 5, and never carries out an asymptotic
analysis for it. The paper closes that gap: it derives the functional
equations explicitly for fixed *m*, obtains a **closed-form expression**
for the critical growth-rate value when *m* is odd, **proves** that the
same closed-form approach is structurally obstructed when *m* is even
(not merely that it hasn't been found), and carries out a full singularity
analysis in both cases. It also conjectures, with numerical support, that
the exponential growth rate decreases strictly as *m* increases. The
resulting enumerative data has also been submitted to the OEIS (see
"Relation to existing OEIS arrays" and "Data availability" below).

The scripts below independently generate and cross-check every numerical
claim in the paper, by several unrelated methods, so that a reader does not
have to take any of the enumerative or asymptotic data on faith.

Independent verification scripts accompanying:

> Fr. G. Speyser, *Enumeration and Asymptotic Analysis of Strict Non-Plane m-Gonal Cactus Graphs via Split-Decomposition*, 2026.
> Submitted to *The Electronic Journal of Combinatorics*.
> Preprint: https://doi.org/10.5281/zenodo.21513752

## Core scripts

- **`mgonal_cactus_series.py`** — computes the rooted and unrooted enumeration
  series for strict *m*-gonal cactus graphs (free / non-plane case), for
  *m* = 5, 6, 7, 8, directly from the defining functional equations of the
  split-decomposition grammar (Sections 5.1–5.3 of the paper), using exact
  rational formal power series arithmetic (Python `Fraction`). This is the
  script that generated the enumerative data tabulated in the paper and
  submitted to the OEIS. All other scripts below verify this one's output
  by independent means.
- **`mgonal_cactus_growth_rate.py`** — estimates the exponential growth rate
  1/ρ*m* directly from the coefficients computed above, via an
  *n*−3/2-corrected ratio test, independent of the analytic
  criticality system (Theorems 2–4). Used as supporting numerical evidence
  for Conjecture 1.

## Supplementary verification scripts

These go beyond what the paper itself requires; they were added as further,
independent cross-checks and are not needed to follow the paper. They were
developed progressively, one question at a time, as the paper and its OEIS
submissions were prepared — which is why several of them target the same
underlying data by genuinely different computational routes: that
redundancy is the point, not an oversight, and the CHANGELOG below records
the order in which each check was added.

- **`verify_pari.gp`** (PARI/GP) — recomputes the rooted series for all four
  *m* using PARI's **native truncated power series arithmetic**, solving
  the functional equation of §5.1 directly (a different code path from
  `mgonal_cactus_series.py`'s Fraction-based recursion). Also verifies
  Theorem 2's closed form for τ*m* (*m* odd) by independent numerical
  root-finding.
- **`verify_pari_euler.gp`** (PARI/GP) — a *second*, differently-structured
  PARI verification of the same four rooted series, this time via an
  **explicit Euler-transform recurrence**: the same principle Andrew
  Howroyd used in his own PARI code on A398033 (m=5), generalized here to
  all four values of *m* and both parities, with a built-in self-check
  against the known first ten terms of each. Kept alongside `verify_pari.gp`
  deliberately: the two use unrelated internal machinery (native series
  solving vs. explicit combinatorial recurrence) even though both are
  PARI/GP, so agreement between them is a genuine cross-check, not a
  restatement of the same computation in the same language.
- **`verify_pari_dissymmetry_odd_m.gp`** (PARI/GP) — verifies the
  dissymmetry-theorem decomposition of §5.3 (G(x) = T_Cm(x) + T_S(x) -
  T_{S-Cm}(x)) by yet another, unrelated route from
  `verify_dissymmetry_m6.py` / `verify_dissymmetry_all_m.py`: starting
  from Andrew Howroyd's own general PARI code on A332649/A332648, it
  specializes his U(n,k) formula to m=7 and shows it collapses, for any
  odd *m*, to a two-term closed form a(n) = g1 + (3/m)·x·(g_m − g1^m).
  Written in response to an OEIS editor's request (Sean A. Irvine, on
  A398575) for a formula checkable without unpacking the paper's
  T_Cm/T_S/T_{S-Cm} notation directly. Confirms, term for term, that the
  specialized closed form and Howroyd's unmodified general formula agree
  exactly (all 13 terms checked), and that both reproduce the published
  A398575 data.
- **`verify_pari_dissymmetry_even_m.gp`** (PARI/GP) — the even-*m*
  companion to the script above, for m=8. Unlike the odd case, **no
  two-term closed form exists here**: 8 has four divisors rather than
  two, and Howroyd's general formula carries an extra term active only
  for even *k*, neither of which cancels the way they do for odd *m* —
  the same parity obstruction that runs through the whole paper
  (Proposition 1), showing up again at this level. The script verifies
  m=8 unrooted data directly against Howroyd's unmodified general
  formula and against `mgonal_cactus_series.py`'s independent exact-
  rational solver (12 terms, exact agreement), without claiming a
  simplification the even case does not structurally admit.
- **`exhaustive_iso_m8.py`** — a *different kind* of check for m=8
  unrooted, since no algebraic simplification was available: builds
  strict 8-gonal cacti with 1, 2, and 3 blocks directly as graphs (no
  generating function, no PARI, no functional equation at all) and
  deduplicates by graph isomorphism (via `networkx`), recovering the
  counts 1, 1, 5 independently — the same method `exhaustive_iso.py`
  already uses for m=5, applied here to the one case (m=8, unrooted)
  that has no closed-form verification available.
- **`split_tree_v2.py`** — a from-scratch, brute-force split-decomposition
  search (Definition 1), used to test Theorem 1 directly on small graphs:
  positive cases (genuine strict *m*-gonal cacti) and negative cases (a
  chord added inside a block, a bridge between two blocks, a cycle of the
  wrong length) — checking that the characterization's condition (a) holds
  or fails exactly as it should.
- **`exhaustive_iso.py`** — builds strict 5-gonal cacti with 1, 2, and 3
  blocks directly as graphs (no functional equation involved at all) and
  deduplicates by graph isomorphism (via `networkx`), recovering the counts
  1, 1, 3 independently of any of the series computations above.
- **`asymptotic_convergence.py`** — checks empirically that s_n really
  converges to C_m·ρ_m^(-n)·n^(-3/2) as n grows (Theorems 3-4). Covers
  *m* = 5, 6, 7, 8 — one odd/even pair was verified first (5 and 6, then
  8), and *m* = 7 was added subsequently to complete the coverage of every
  value of *m* the paper actually treats numerically, for both parities.
  Includes a documented finding: a naive double-precision run shows a
  spurious uptick in the ratio beyond n ≈ 1000 for m=5, which disappears
  entirely at 60-digit precision — a floating-point artifact, not a real
  secondary term.
- **`verify_lemma5.py`** — independently reproduces the three numerical
  transition values (1.045, 0.855, 0.971) quoted in the paper's discussion
  of Conjecture 1, computed from the τ_m-substitution described in the
  text (see the script's docstring for why this differs from the exact,
  currently-unusable Lemma 5 criterion itself).
- **`verify_boltzmann.py`** — independently re-implements the Boltzmann
  sampler of Section 5.4 (the Burnside-style identity/reflection stabilizer
  choice, built recursively rather than sampled from the known coefficient
  distribution) and validates it against the exact enumerative data for
  m=5. A stated simplification (the dominant i=1 term of the MSET
  construction only) is documented in the script.
- **`verify_lemmas234.py`** — independent reconstruction of the graph G
  from a graph-labeled tree satisfying Theorem 1's conditions (a)-(d),
  checking Lemmas 2 and 4 directly (and Lemma 3 implicitly) on a 2-block
  and a 3-block example. The script's docstring records a genuine
  construction bug found and fixed during this verification (in the
  script itself, not the paper) — a useful cross-check of how precisely
  the paper's stated construction must be followed.
- **`verify_dissymmetry_m6.py`** — the *first* version of the dissymmetry-
  theorem verification (Section 5.3), written specifically for m=6
  (A398035) while that pair of sequences was being prepared: assembles
  G(x) = T_Cm(x) + T_S(x) - T_{S-Cm}(x) from scratch in SymPy, including
  an explicit symbolic implementation of the dihedral cycle index Z_D6
  (equation 8, even-*m* branch), using only the already-verified rooted
  series as input.
- **`verify_dissymmetry_all_m.py`** — the *generalized* version of the
  script above: the same dissymmetry-theorem assembly, but written to
  handle both parities and all four values of *m* the paper treats (5, 6,
  7, 8) in a single script, with a self-check against the known unrooted
  data for each. Kept alongside `verify_dissymmetry_m6.py` as the record
  of how the general version was arrived at, one case at a time.

## Usage

```bash
python3 mgonal_cactus_series.py
python3 mgonal_cactus_growth_rate.py
gp -q verify_pari.gp
gp -q verify_pari_euler.gp
gp -q verify_pari_dissymmetry_odd_m.gp
gp -q verify_pari_dissymmetry_even_m.gp
python3 exhaustive_iso_m8.py  # requires: pip install networkx
python3 split_tree_v2.py
python3 exhaustive_iso.py    # requires: pip install networkx
python3 asymptotic_convergence.py
python3 verify_lemma5.py
python3 verify_boltzmann.py
python3 verify_lemmas234.py
python3 verify_dissymmetry_m6.py      # requires: pip install sympy
python3 verify_dissymmetry_all_m.py   # requires: pip install sympy
```

No dependencies beyond the Python standard library, `numpy`, `sympy`, and
`networkx` (for `exhaustive_iso.py` and `exhaustive_iso_m8.py`). The four
`.gp` scripts require PARI/GP.

## Relation to existing OEIS arrays

The eight sequences computed here coincide, term for term, with columns
*k* = 5, 6, 7, 8 of two general arrays already on the OEIS:
[A332648](https://oeis.org/A332648) (rooted case) and
[A332649](https://oeis.org/A332649) (unrooted case), both by Andrew Howroyd.
This was identified during the OEIS submission process and independently
re-verified here by two methods: direct comparison against the published
data of A332648/A332649, and a symbolic recomputation using a different
algorithm from the one in `mgonal_cactus_series.py`. Both confirm exact
agreement.

An OEIS editorial review (Sean A. Irvine, August 2026) further noted that,
absent a formulation specific to a single column, standard practice for
arrays such as A332649 is not to create individual column entries beyond
*k* = 4 — a point the array's own CROSSREFS confirms. For the two rooted
entries (A397210, A397546), the paper's Theorem 2 / Proposition 1 (a
closed-form critical value for odd *m*, provably unreachable by the same
method for even *m*) was judged to supply exactly this kind of
column-specific content, and both entries are defended on that basis.

For the unrooted case, the same test was applied honestly to both parities,
with two different outcomes:

- **m=7 (odd):** the dissymmetry-theorem decomposition collapses to an
  explicit two-term closed form (`verify_pari_dissymmetry_odd_m.gp`),
  giving A398575 exactly the kind of column-specific, independently
  checkable content Irvine's request was looking for.
- **m=8 (even):** no such collapse exists — verified directly
  (`verify_pari_dissymmetry_even_m.gp`), not merely absent from a search.
  The data is independently confirmed correct by two unrelated methods
  (Howroyd's general formula, and brute-force graph construction in
  `exhaustive_iso_m8.py`), but correctness is not the same question as
  whether a standalone OEIS entry is warranted: with no column-specific
  formulation, the honest answer is that it is not, on the same standard
  applied to m=7. **This sequence is therefore not submitted to the OEIS**
  as a standalone entry, and none is currently planned. The data and both
  verification scripts live in this repository so the computation remains
  available, reproducible, and citable regardless of its OEIS status —
  useful, for instance, should A332649 ever be extended with a dedicated
  b-file for this column instead of a standalone sequence.

## Data availability

The integer sequences computed by `mgonal_cactus_series.py` have been
submitted to the OEIS (On-Line Encyclopedia of Integer Sequences), or are
in the process of being submitted, with one deliberate exception (see
below). Status so far:

- *m* = 5, rooted: [A398033](https://oeis.org/A398033) (approved)
- *m* = 5, unrooted: [A397250](https://oeis.org/A397250) (approved)
- *m* = 6, rooted: [A398034](https://oeis.org/A398034) (approved)
- *m* = 6, unrooted: [A398035](https://oeis.org/A398035) (approved)
- *m* = 7, rooted: [A397210](https://oeis.org/A397210) (approved;
  defended as a standalone entry on the strength of Theorem 2 / Proposition 1)
- *m* = 7, unrooted: [A398575](https://oeis.org/A398575) (proposed, pending
  review — the checkable closed-form PARI program in
  `verify_pari_dissymmetry_odd_m.gp` was prepared in response to an
  editor's request and is pending submission to the PROG field)
- *m* = 8, rooted: [A397546](https://oeis.org/A397546) (approved;
  same defense as A397210)
- *m* = 8, unrooted: **not submitted to the OEIS, and not planned.**
  Unlike its odd-*m* counterpart, this column admits no closed-form
  simplification (see `verify_pari_dissymmetry_even_m.gp`), so it does not
  meet the same column-specific-content standard applied to A398575. The
  data is independently computed and verified by two unrelated methods
  (`verify_pari_dissymmetry_even_m.gp` and `exhaustive_iso_m8.py`), and is
  kept in this repository for transparency and reuse — e.g. as a
  candidate b-file addition to A332649 — rather than as a standalone
  sequence.

## Acknowledgments

Thanks to Andrew Howroyd, a long-standing OEIS editor whose own work centers
on Pólya enumeration, for identifying the relation to A332648/A332649 (both
his own OEIS arrays) and for an independent PARI verification of the data;
to Sean A. Irvine for a close editorial reading that prompted a sharper,
parity-by-parity account of why some of these entries are defensible as
standalone sequences and others are not, and for the resulting closed-form
PARI specialization in `verify_pari_dissymmetry_odd_m.gp`; and to Stefano
Spezia, Robert C. Lyons, and Michael De Vlieger for their review of the
OEIS submissions above. Further contributors may be added here as the
review of the remaining submissions proceeds.

## Citation

If you use this code, please cite the paper above. A citable archive of this
repository is available via Zenodo: https://doi.org/10.5281/zenodo.21461100

## Author

Frédéric G. Speyser — Independent Researcher, France - https://orcid.org/0000-0002-1767-5325

## License

MIT (see `LICENSE`).
