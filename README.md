# SDDF — Spectral Curvature G[u]

Extension of **LOGOS**, within [`vega-vault`](https://github.com/Rylow999/vega-vault).

Study of spectral curvature

$$G[u] = \int \left(rac{d\ln E(k)}{d\ln k}ight)^2 d\ln k$$

as a diagnostic of the Kolmogorov inertial range, with exact closed form, a bias-corrected spectral exponent estimator, and a test of the log-periodic signature proposed by Migdal (arXiv:2604.12207, arXiv:2511.02165) for decaying turbulence.

**Key result:** $G^*[u] = rac{25}{12}\ln Re + b(\delta) + \mathcal{O}(Re^{-1})$, with prefactor $25/12=(25/9)\cdot(3/4)$ — the square of the Kolmogorov slope by the scaling exponent of $\eta\propto Re^{-3/4}$. Not a fit.

## Overview

The SDDF (Spectral Density Distribution Function) framework extends **LOGOS** within the `vega-vault` repository (see [Rylow999/vega-vault](https://github.com/Rylow999/vega-vault)). It provides tools for spectral curvature analysis, bias-corrected spectral exponent estimation, and log-periodic signature detection in decaying turbulence. The central observable $G[u]$ serves as a diagnostic of the Kolmogorov inertial range with exact closed-form theory and practical numerical estimation pipelines.

## Repository structure

```.
├── README.md              this file
├── AUDITORIA.md           complete history of review rounds (v1→v2→v3→A7)
├── CONTRASTE_MIGDAL_DNS.md complete guide to contrast the Migdal signature in real DNS
├── paper/
│   └── PAPER.md           consolidated document with all results
├── codigo/
│   ├── sddf_core.py       core numerics (reading, local slope, truncation)
│   ├── sddf_exact.py      closed-form analytic shapes
│   ├── estimador_principal.py  recommended estimator (q_corrected), reusable API
│   ├── modelos_espectrales.py  Pao / Pope / chimera / 2D parametrizations
│   ├── barrido_re_v3.py   generates CSV spectra from data/spectra
│   ├── sddf_completo.py   full pipeline (regenerates all CSV/PNG)
│   └── rust/              reference implementation in Rust (compiled)
├── datos/                 15 result tables + 27 spectral CSV files (includes spectrum_grueso.csv, recovered)
├── figuras/               10 figures (PNG)
└── tests/
    └── test_estimador_principal.py   regression against the paper
```

## Key results

- **Closed-form theory:** $G^*[u] = rac{25}{12}\ln Re + b(\delta) + \mathcal{O}(Re^{-1})$, with $25/12 = (25/9)\cdot(3/4)$, derived from $(5/3)^2 \cdot 3/4$. Numerical verification shows 0.06-0.07% error across $Re=100..2000$, with prefactor desviscosado $2.0832$ vs exacto $2.0833$ (bias $2.6\%$ measured and explained via finite-$\eta$ correction).
- **DNS validation (JHTDB box):** Evidence consistent with the theoretical regime: $q = 1.60 \pm 0.02$ for $k\eta \in [0.02, 0.18]$, $\sim 4\%$ below K41 ($5/3=1.667$). Including $k \leq 128$ (onset of dissipation) $q$ rises to $1.76	ext{--}1.92$ — the window matters.
- **Bias-corrected spectral estimator:** The recommended `estimador_principal` provides $q_{corrected}$ with documented pipeline; validation on synthetic noise + DNS real confirms intermitency with correct $eta$: error $2	imes10^{-7}$.
- **Spectral detector (smoothed):** The smoothed slope detector with two-sided inertial window is conservative by design: excludes $window/2$ border points, uses mobile regression (`smoothed_slope`), and validated against synthetic noise ($0	ext{--}10\%$) and DNS real. For DNS data with $\sim 2\%$ concha-to-concha SEM in 8-block averaging, $w=5$ is recommended; $w=9	ext{--}11$ for noisier spectra.
- **Migdal log-periodic test: negative honest result.** Applied to real JHTDB data (box $256^3$, iso $1024^3$) and own DNS 2D decay: no detectable signature. Bootstrap null (600 bootstraps, concha uncertainty $E_{sem}$) gives $p=0.43	ext{--}0.45$, with $\omega_{obs}$ at the grid edge $\Omega\in[0.3,8]$, identifying the cubic detrending artifact. Synthetic calibration (v3.2) is reinforced; the signature requires greater inertial range resolution → motivates $4096^3$.
- **Hypo-dissipative extension:** Closed-form observable with generalized dissipation $\gamma$: $G^*_\gamma = rac{25}{16}\gamma\ln Re + b(\gamma)$. $\gamma=2$ is the best fit. Explicitly not a claim of Navier-Stokes solution without forcing — declared.
- **DNS 2D own solver:** Pseudo-spectral standard solver (Boffetta-Ecke: 2/3 de-alias, RK2, exact dissipation+friction per integration step, constant physical forcing $F_0$ on shells $k\in[3,5]$). Three hard-learned physics lessons documented in code: (1) $\hat{\psi} = +\hat{\omega}/k^2$ — the negative sign inverts $u,v$ and reverses enstrophy cascade; (2) forcing amplitude is PHYSICAL (Parseval: $|f|_{	ext{rms}} = f_0$ requires $\sum|f_k|^2 = (f_0\cdot N^2)^2$) — with unit-norm $f_k$ the forcing was $\sim10^{-5}$ and the cascade never filled the spectrum; (3) $\chi pprox f_0^2$ (not 1): real dissipation scale $k_d \sim (f_0^2/
u^3)^{1/6}$ — $Re_{	ext{sat}} \sim 9400$ at $N=128$. Window declaration $k\in[6,20]$ (apriori, same for all $
u$) gives $G^* \in [24.4, 29.9]$ with $a = dG^*/d\ln(1/
u) = -6.74$, negative sign contrary to prediction $+4.5$ — the **window contamination artifact** is the finding: it explains the "impossible" $G^*pprox3634$ from the external 2D paper structurally, with our own data.
- **Key negative results honestly labeled:** (1) Migdal log-periodic not detectable in 2D own decay (14 times, same score) nor in real 3D DNS (box: score 22.5, null med 18.7, $p=0.45$; iso$1024$: score 47.7, null med 45.1, $p=0.43$). (2) Automatic detector fails on real data: spikes with concha-to-concha noise ($\sim10	ext{--}15\%$ at low $k$) and cuts 3-8 points; two-sided inertial window finds no valid window. Requires smoothed slope (window 5) for real data. (3) DNS 2D window artifact explains structurally the "impossible" $G^*pprox3634$ from the external 2D paper — the window declaration without a physically defined cutoff contaminates the integral differently at each $
u$.

## Reproduction

```bash
cd codigo
python3 barrido_re_v3.py ../datos/espectros
python3 sddf_completo.py ../datos/espectros ../datos ../figuras
./rust/bin/sddf ../datos/espectros/spectrum_Re_1000.csv 0.5
cd .. && python3 tests/test_estimador_principal.py
```

## Using the estimator on your own spectrum

```python
from estimador_principal import medir_exponente_espectral

# k, E: arrays of your spectrum. eta: Kolmogorov scale (nu^3/eps)^(1/4).
# beta: model constant for your cutoff (do NOT assume 5.2 from Pope without verifying — see AUDITORIA.md, findings A3/A6).
q, mu, k_corte = medir_exponente_espectral(k, E, eta, delta=0.10, beta=2.25)
```

## Citation

For citations, refer to `paper/PAPER.md` sección 11 for the full document with all results, including the closed-form derivation, DNS validation evidence, and honest labeling of negative results. The central formula $G^*[u] = rac{25}{12}\ln Re + b(\delta) + \mathcal{O}(Re^{-1})$ derives from $(5/3)^2 \cdot 3/4$ and is verified numerically across $Re=100..2000$.

## Related work

The observable $G[u]$ shares the "observer collapse" pattern with the $ho$-law of the repo [Rylow999/fhrr-rho-collapse](https://github.com/Rylow999/fhrr-rho-collapse): when the spectrum is measured on a discrete grid, the truncation error does not decay monotonically — it oscillates. The substrate (the flow) has the structure; the instrument (the grid) is the one that distorts it near its resolution limit. The multi-observer experiment (repo [Rylow999/rho-law](https://github.com/Rylow999/rho-law)) shows that there are invariant quantities between observers (ratios between windows, $CVpprox0.24$) and quantities that are artifacts ($G_{	ext{total}}$, $CVpprox0.7$). See [NOUS/RHO_LAW/](https://github.com/Rylow999/NOUS/RHO_LAW/) for the unifying cross-domain framework.

## Status & honest labeling

Refer to `AUDITORIA.md` for the complete history and `paper/PAPER.md` sección 11 for open limitations (recovery of `spectrum_grueso.csv` and recalculation of the 2D paper, both blocked by lack of input data, not by method). Key honest labelings:

- The central formula is **algebraically proven** (theorem); DNS evidence is consistent with the physical regime but does not prove the theorem universally.
- All three Migdal test results are formally labeled as **honest negative results** (type "medium"): they eliminate a family of predictions and delimit the validity regime of the program, not a failure.
- The DNS 2D window artifact is the structural explanation for the "impossible" $G^*pprox3634$ from the external 2D paper — the window declaration without a physically defined cutoff contaminates the integral differently at each $
u$.
- The $arepsilon$ "spectral" reported in earlier versions ($9.46	imes$ the documented value integrating the full range) is an artifact of a non-periodic $256^3$ spectrum with Hanning window; the gradient-validated $arepsilon$ ($0.96	imes$ documented) is the reliable quantity.

## Motto

*Per Aspera, Ad Astra.*
