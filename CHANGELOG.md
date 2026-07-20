# Changelog

All notable changes to this repository are documented in this file, with a
date on each entry.

## [Unreleased]

(nothing yet)

## [1.1] - 2026-07-20

### Added
- `verify_boltzmann.py` — independent re-implementation of the Boltzmann
  sampler described in Section 5.4 of the paper (Burnside-style
  identity/reflection stabilizer choice, built recursively), validated
  against the exact enumerative data for m=5.

### Changed
- `README.md` — added description and usage line for `verify_boltzmann.py`.

## [1.0] - 2026-07-20

### Added
- `mgonal_cactus_series.py` — computes rooted and unrooted enumeration
  series for strict *m*-gonal cactus graphs, *m* = 5, 6, 7, 8, directly
  from the split-decomposition functional equations.
- `mgonal_cactus_growth_rate.py` — independent estimate of the exponential
  growth rate 1/ρ*m* from the coefficients above.
- `verify_pari.gp` — independent PARI/GP recomputation of the rooted series
  (native power series arithmetic) and of Theorem 2's closed form for
  tau_m, m odd.
- `split_tree_v2.py` — brute-force split-decomposition search testing
  Theorem 1 directly on small graphs, including negative cases (chord,
  bridge, wrong cycle length).
- `exhaustive_iso.py` — direct graph construction and isomorphism-based
  counting (k=1,2,3 blocks), independent of the functional equation
  entirely.
- `asymptotic_convergence.py` — empirical verification that s_n converges
  to C_m*rho_m^(-n)*n^(-3/2) as n grows, for m=5 and m=6.
- `verify_lemma5.py` — independent recomputation of the three numerical
  transition values quoted in the Conjecture 1 discussion.
- README rewritten: introduction to the research, OEIS accession numbers,
  relation to A332648/A332649, acknowledgments.

### Notes
- Data submitted to the OEIS: A398033, A397250, A398034, A398035 (under
  review, none formally accepted yet); m = 7, 8 sequences prepared,
  submission pending.
- A398034/A398035 (m=6) held in "editing" status pending Andrew Howroyd's
  go-ahead, per his explicit request.
- Preprint and Zenodo archive: not yet available (see README).
