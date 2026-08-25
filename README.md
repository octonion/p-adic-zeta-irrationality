# Hybrid arithmetic holonomy and twenty-two individual $p$-adic zeta values

**Christopher D. Long**

This repository contains the LaTeX source, compiled manuscript, and exact rational fixed-point interval certificate for the arithmetic-holonomy argument developed in the paper.

## Abstract

We develop a hybrid finite-place refinement of arithmetic holonomy for the Eisenstein--Eichler connection on the genus-zero modular curves $X_0(p)$, $p\in\{2,3,5,7\}$. The denominator staircase of the Eichler primitive cannot be realized by globally defined scalar solution coordinates: the nonzero modular weights create a genuine torus obstruction. We instead use two different local mechanisms. At auxiliary primes below a cutoff, the Hasse invariant algebraizes the unipotent comparison coordinate modulo the prime and a genuine-source rank--nullity argument produces one determinant-level saving on a large source kernel. Above the cutoff, we keep the calculation on the algebraic de Rham frame torsor; exact divided-Frobenius relations for the raw Eichler rows and a coefficient-ring Taylor no-backflow lemma recover the Calegari--Dimitrov--Tang prime window without scalarizing the packet. The resulting arithmetic cost is

$$
\tau_{p,s}(\xi)=\frac{2}{(s+1)^2}\left(S_{p,s}(\xi)+s I^{s+1}_{\xi}(\xi)\right),
$$

where $S_{p,s}$ is the small-prime Hasse-codimension cost and $I^{s+1}_{\xi}$ is the published prime-window function. Combining this with Calegari's overconvergent Eisenstein family, Buzzard's analytic continuation, an exact modular Jensen formula, and Bost's slopes inequality proves the irrationality of at least twenty-two individual Kubota--Leopoldt values:

$$
\zeta_2(s)\quad(3\le s\le29,\ s\text{ odd}),\qquad
\zeta_3(s)\quad(3\le s\le11,\ s\text{ odd}),
$$

$$
\zeta_5(3),\quad \zeta_5(5),\quad \zeta_7(3).
$$

Every numerical inequality is supplied with an independent rational fixed-point interval certificate. No opposite Hodge filtration, rational Tate scalarization, or rowwise transport of literal $n^e$ denominators is used.

## Repository files

| File | Description |
| --- | --- |
| [`hybrid_arithmetic_holonomy.tex`](hybrid_arithmetic_holonomy.tex) | Complete LaTeX source of the manuscript. |
| [`hybrid_arithmetic_holonomy.pdf`](hybrid_arithmetic_holonomy.pdf) | Compiled PDF of the manuscript. |
| [`hybrid_rational_interval_certificate.py`](hybrid_rational_interval_certificate.py) | Standard-library Python certificate using exact `Fraction` arithmetic and outward-rounded integer fixed-point intervals. It reconstructs all 22 numerical margins independently from the formulas in the paper. |
| [`hybrid_rational_interval_certificate_output.txt`](hybrid_rational_interval_certificate_output.txt) | Stored output of the rational interval certificate for all 22 cases. |
| [`hybrid_rational_interval_certificate_sha256.txt`](hybrid_rational_interval_certificate_sha256.txt) | SHA-256 manifest for the certificate script and stored certificate output. |

## Main theorem in the current draft

The manuscript proves, subject to the hypotheses and argument developed there, the irrationality of

$$
\zeta_2(s)\quad\text{for }s=3,5,\ldots,29,
$$

$$
\zeta_3(s)\quad\text{for }s=3,5,\ldots,11,
$$

and

$$
\zeta_5(3),\qquad \zeta_5(5),\qquad \zeta_7(3).
$$

Because these conclusions include individual irrationality statements beyond the classical cases, independent specialist review remains appropriate before publication.

## Reproducing the manuscript

A standard LaTeX installation with the packages listed in the source is sufficient. For example:

```bash
pdflatex hybrid_arithmetic_holonomy.tex
pdflatex hybrid_arithmetic_holonomy.tex
pdflatex hybrid_arithmetic_holonomy.tex
```

The resulting PDF should be `hybrid_arithmetic_holonomy.pdf`.

## Reproducing the rational interval certificate

The certificate requires only Python 3 and the Python standard library:

```bash
python3 hybrid_rational_interval_certificate.py \
  > hybrid_rational_interval_certificate_output.txt
```

The script performs the algebraic calculations with exact rational arithmetic. Transcendental quantities are enclosed using arbitrary-precision integer fixed-point intervals with explicit remainder bounds; no floating-point or external interval package is used in the certification step.

To verify the stored script and output against the repository manifest:

```bash
sha256sum -c hybrid_rational_interval_certificate_sha256.txt
```

The smallest certified margin is the $(p,s)=(5,5)$ case:

$$
\mathfrak M_{5,5}\left(\frac{29}{27},\frac{1}{16}\right)
\in
\left[
0.131799356827016832557664457131479890735826671587851068262655,
0.131799356827016832557664457131479890735826671587851068262656
\right].
$$
