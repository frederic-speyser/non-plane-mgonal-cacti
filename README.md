# Enumeration of Non-Plane m-Gonal Cactus Graphs — verification code

Independent verification scripts accompanying:

> Fr. G. Speyser, *Enumeration of Non-Plane m-Gonal Cactus Graphs via Split-Decomposition*.
> Submitted to *The Electronic Journal of Combinatorics*.
> Preprint: arXiv:XXXX.XXXXX (Ref. TBA)

## What this repository contains

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

## Usage

```bash
python3 mgonal_cactus_series.py
python3 mgonal_cactus_growth_rate.py
```

No dependencies beyond the Python standard library (`fractions`).

## Data availability

The integer sequences computed by `mgonal_cactus_series.py` are submitted to
the OEIS (On-Line Encyclopedia of Integer Sequences); accession numbers will
be added here once assigned.

## Citation

If you use this code, please cite the paper above. A citable archive of this
repository is available via Zenodo: XXXXXXXXXXXXX (Ref. TBA).

## Author

Frédéric G. Speyser — ORCID: 0009-0006-2989-656X

## License

MIT (see `LICENSE`).
