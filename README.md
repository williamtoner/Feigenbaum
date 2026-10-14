# Feigenbaum constants: record expansions

The Feigenbaum constants **δ** and **α** for unimodal maps x ↦ 1 − μ|x|^z with a critical point of degree z = 2, 4, 6, 8.

`index.html` is the project page, with interactive bifurcation diagrams, the renormalisation operator, records and a digit explorer.

| Constant | Value | Digits | Previous best published | File |
|---|---|---:|---|---|
| δ, degree 2 | 4.66920160910299067185… | **10,000** | 1,019 (Broadhurst 1999)* | `digits/delta_degree2.txt` |
| α, degree 2 | 2.50290787509589282228… | **10,000** | 1,018 (Broadhurst 1999)* | `digits/alpha_degree2.txt` |
| δ, degree 4 | 7.28468621707334336430… | **5,000** | 325 proven (Burbanks, Osbaldestin & Thurlby 2021) | `digits/delta_degree4.txt` |
| α, degree 4 | 1.69030297140524485334… | **5,000** | 331 proven (Burbanks, Osbaldestin & Thurlby 2021) | `digits/alpha_degree4.txt` |
| δ, degree 6 | 9.29624683277137008283… | **5,000** | 46 (Briggs 1991) | `digits/delta_degree6.txt` |
| α, degree 6 | 1.46774245031990094445… | **5,000** | 59 (Mathar 2010) | `digits/alpha_degree6.txt` |
| δ, degree 8 | 10.9486242659415904255… | **5,000** | 38 (Briggs 1991) | `digits/delta_degree8.txt` |
| α, degree 8 | 1.35801727913805034548… | **5,000** | 42 (Briggs 1991) | `digits/alpha_degree8.txt` |

\* Molteni (2016) reported about 10,000 digits of the classical constants but did not publish them.

The digit files use OEIS b-file format: one `n digit` per line, with n starting at the number of digits before the decimal point.

## Method

The fixed point g of the period-doubling renormalisation operator satisfies g(1)·g(x) = g(g(g(1)·x)), and α = −1/g(1).
- **Representation.** Write g(x) = f(|x|^z) and expand f in Chebyshev polynomials.
- **Solve.** Enforce the equation at n collocation nodes, then solve the n×n system by Newton iteration in high-precision fixed-point arithmetic.
- **δ.** δ is the leading eigenvalue of the operator linearised at g, and is solved the same way.

This extends Molteni's method (arXiv:1602.02357) to general even degree.

## Verification

These are computed digits, not a proof.
1. **Two discretisations.** Each constant was computed twice, with n and about 1.1n nodes and different precisions. The runs agreed to 10,066+ digits (degree 2) and 5,067+ digits (degrees 4, 6, 8). The files stop 10 digits short of that.
2. **Independent implementation.** `src/alt.py` uses a Legendre basis, different nodes and normalisation, Arb ball arithmetic and full Newton steps, and shares no code with `src/fb.py`. It agrees on every digit it computed: 3,015 for degree 2 and 1,514 or 1,515 for degrees 4, 6, 8.
3. **Literature.** The digits agree with every previously published digit, including the 325 and 331 rigorously proven digits for degree 4.

## Code

```
src/fb.py          main solver (Chebyshev collocation)
src/alt.py         independent check (Legendre basis, Arb ball arithmetic)
src/proto.py       small mpmath prototype (used to bootstrap fb.py)
src/compare.py     count agreeing digits between result files / digit files
src/site_data.py   rebuild assets/data.js for index.html from digits/
```

```bash
uv sync
# compute degree-4 alpha and delta to 600 digits (about a minute on 8 cores) and check against the files
uv run python src/fb.py --z 4 --digits 600 --workers 8 --out z4.json
uv run python src/compare.py z4.json digits/alpha_degree4.txt alpha
uv run python src/compare.py z4.json digits/delta_degree4.txt delta
# independent implementation
uv run python src/alt.py --z 4 --digits 200 --threads 4 --out z4_alt.json
uv run python src/compare.py z4_alt.json digits/delta_degree4.txt delta
```

The full runs (`--digits 10050` for degree 2, `--digits 5050` for the others, plus a check run with `--nfactor 1.1`) took a few hours each on about 30 cores.

## References

- M. J. Feigenbaum, J. Stat. Phys. 19 (1978) 25–52.
- K. Briggs, *A precise calculation of the Feigenbaum constants*, Math. Comp. 57 (1991) 435–439.
- R. J. Mathar, arXiv:1008.4608 (2010).
- A. Molteni, *An efficient method for the computation of the Feigenbaum constants to high precision*, arXiv:1602.02357 (2016).
- A. D. Burbanks, A. H. Osbaldestin, J. A. Thurlby, J. Math. Phys. 62 (2021) 112701.
- M. Breden, J. Gonzalez, J. D. Mireles James, arXiv:2409.20457 (2024).
