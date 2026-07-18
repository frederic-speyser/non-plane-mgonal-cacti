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

Thanks to Andrew Howroyd, a long-standing OEIS editor whose own work centers
on Pólya enumeration, for identifying the relation to A332648/A332649 (both
his own OEIS arrays) and for an independent PARI verification of the data.
Thanks also to Stefano Spezia, Robert C. Lyons, and Sean A. Irvine for their
review of the OEIS submissions above.

## Citation

If you use this code, please cite the paper above. A citable archive of this
repository is available via Zenodo: XXXXXXXXXXXXX (Ref. TBA).

## Author

Frédéric G. Speyser — ORCID: 0009-0006-2989-656X

## License

MIT (see `LICENSE`).
