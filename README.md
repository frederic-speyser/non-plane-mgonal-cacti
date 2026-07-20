

Readme · MD
# Enumeration of Non-Plane m-Gonal Cactus Graphs — verification code
 
Independent verification scripts accompanying:
 
> Fr. G. Speyser, *Enumeration of Non-Plane m-Gonal Cactus Graphs via Split-Decomposition*.
> Submitted to *The Electronic Journal of Combinatorics*.
> Preprint: arXiv:XXXX.XXXXX (Ref. TBA)
 
## Core scripts
 
- **`mgonal_cactus_series.py`** — computes the rooted and unrooted enumeration
  series for strict *m*-gonal cactus graphs (free / non-plane case), for
  *m* = 5, 6, 7, 8, directly from the defining functional equations of the
  split-decomposition grammar (Sections 5.1–5.3 of the paper), using exact
  rational formal power series arithmetic. Running it reproduces the
  coefficients tabulated in Sections 6.1–6.2 of the paper and prints the
  first 25 terms of each of the eight sequences submitted to the OEIS.
- **`mgonal_cactus_growth_rate.py`** — estimates the exponential growth rate
  1/ρ*m* directly from the coefficients computed above, via an
  *n*<sup>−3/2</sup>-corrected ratio test, independent of the analytic
  criticality system (Theorems 2–4). Used as supporting numerical evidence
  for Conjecture 1.
## Supplementary verification scripts
 
These go beyond what the paper itself requires; they were added as further,
independent cross-checks and are not needed to follow the paper.
 
- **`verify_pari.gp`** (PARI/GP) — recomputes the rooted series for all four
  *m* using PARI's native truncated power series arithmetic (a different
  code path from the two Python scripts above), and verifies Theorem 2's
  closed form for τ*m* (*m* odd) by independent numerical root-finding.
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
  converges to C_m·ρ_m^(-n)·n^(-3/2) as n grows (Theorems 3–4), for both an
  odd *m* (5) and an even *m* (6, the case whose proof required a
  correction — see the paper's working notes). Includes a documented
  finding: a naive double-precision run shows a spurious uptick in the
  ratio beyond n ≈ 1000, which disappears at 60-digit precision — i.e. a
  floating-point artifact, not a real effect.
- **`verify_lemma5.py`** — independently reproduces the three numerical
  transition values (1.045, 0.855, 0.971) quoted in the paper's discussion
  of Conjecture 1, computed from the τ_m-substitution described in the
  text (see the script's docstring for why this differs from the exact,
  currently-unusable Lemma 5 criterion itself).
## Usage
 
```bash
python3 mgonal_cactus_series.py
python3 mgonal_cactus_growth_rate.py
gp -q verify_pari.gp
python3 split_tree_v2.py
python3 exhaustive_iso.py    # requires: pip install networkx
python3 asymptotic_convergence.py
python3 verify_lemma5.py
```
 
No dependencies beyond the Python standard library, `numpy`, and
`networkx` (for `exhaustive_iso.py` only). `verify_pari.gp` requires PARI/GP.
 
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
 
## Data availability
 
The integer sequences computed by `mgonal_cactus_series.py` are submitted to
the OEIS (On-Line Encyclopedia of Integer Sequences). Accession numbers so
far:
 
- *m* = 5, rooted: [A398033](https://oeis.org/A398033)
- *m* = 5, unrooted: [A397250](https://oeis.org/A397250)
- *m* = 6, rooted: [A398034](https://oeis.org/A398034)
- *m* = 6, unrooted: [A398035](https://oeis.org/A398035)
- *m* = 7, 8 (rooted and unrooted): submissions pending (Ref. TBA)
## Acknowledgments
 
Thanks to Andrew Howroyd — a long-standing OEIS editor whose own work centers
on Pólya enumeration — for identifying the relation to A332648/A332649 (both
his own OEIS arrays) and for an independent PARI verification of the data.
Thanks also to Stefano Spezia, Robert C. Lyons, and Sean A. Irvine for their
review of the OEIS submissions above.
 
## Citation
 
If you use this code, please cite the paper above. A citable archive of this
repository is available via Zenodo: XXXXXXXXXXXXX (Ref. TBA).
 
## Author
 
Frédéric G. Speyser — ORCID: https://orcid.org/0000-0002-1767-5325
 
## License
 
MIT (see `LICENSE`).
