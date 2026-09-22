#!/usr/bin/env python3
"""
exp_null_model_periodograma.py  (pendiente #1, README v3.1)

Null model para el test log-periodico de Migdal.

El score actual (P.max() / median(P)) NO es un test estadistico: fabrica
picos "significativos" (SNR ~ 19) incluso sin senal inyectada, porque el
residuo del periodograma es ruido correlado con estructura deterministica
quadratica (la curvatura del corte de disipacion) + artefacto de grilla.

Que hace este script:
  1. Construye el periodograma "observado" sobre espectros con amplitudes
     [0.0, 0.002, 0.005, 0.02, 0.05] (misma ventana fija que el bloque 14
     de sddf_completo.py).
  2. Genera N_NULL espectros nulos (amp=0) con ruido de grilla/debitaje
     numericamente equivalente (semblemos el mismo ruido irreducible que el
     pipeline: jitter log-normal sobre E y jitter en k? -- NO: el pipeline
     usa espectros analiticos EXACTOS. El unico "ruido" es de discretizacion
     numerica de la pendiente local y la misma grilla. Por eso para el null
     usamos amp=0 con el MISMO tratamiento y ademas una familia de nulos con
     beta y ventana jittereadas para capturar sensibilidad instrumental).
  3. Empirical p-value: fraccion de nulos cuyo max score >= score observado.
  4. Presion multiple: correccion Bonferroni sobre las 5 amplitudes, y
     reporte de lambda global = max score del null, que da un umbral de
     deteccion completo (familia de un solo test).

Salida: datos/16_null_model_periodograma.csv + resumen en consola.

Uso:  python3 exp_null_model_periodograma.py [N_NULL] [seed]
"""
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import modelos_espectrales as ME
import sddf_exact as EX
from sddf_core import local_slope

OUT = Path(__file__).parent.parent / "datos" / "16_null_model_periodograma.csv"

# --- mismos parametros que el bloque 14 de sddf_completo.py ---------------
NU_LP = 1e-10                      # nu -> Re = 1e10
K_LP = np.logspace(0, 7, 50000)
ETA_LP = EX.eta_of(NU_LP)
OMEGA_INY = 2.0
AMPS = [0.0, 0.002, 0.005, 0.02, 0.05]
OM_GRID = np.linspace(0.3, 8.0, 1500)
X_LO, X_HI = 1e-6, EX.x_c(0.25, beta=ME.BETA_PAO)


def periodograma(x, y, om):
    """Detrend cubico + |int y e^{-i w x} dx|^2. Identico al pipeline.

    Vectorizado via trapecio cumulativo: T(w) = int y e^{-i w x} dx se
    calcula para todos los w de una vez con la regla del trapecio en forma
    matricial aprovechando la grilla uniforme en x (ln k equiespaciado):
    T(w) = e^{-i w x0} * [ (y0+y1)/2 + sum_{j>=1} y_j e^{-i w j dx} * dip ]
    Usamos FFT seria overkill; basta un producto matriz-vector por bloques.
    """
    y = y - np.polyval(np.polyfit(x, y, 3), x)
    # T(w) = sum_j w_j y_j e^{-i om x_j},  w_j = pesos trapecio
    w = np.empty_like(y)
    w[0] = 0.5 * (x[1] - x[0])
    w[-1] = 0.5 * (x[-1] - x[-2])
    w[1:-1] = 0.5 * (x[2:] - x[:-2])
    z = w * y
    P = np.empty(len(om))
    BLOCK = 128
    for i in range(0, len(om), BLOCK):
        omb = om[i:i + BLOCK]
        ph = np.exp(-1j * np.outer(omb, x))      # (B, n)
        P[i:i + BLOCK] = np.abs(ph @ z) ** 2
    return P


def score_de(e_amp0, k, m):
    """Peak-to-median del periodograma de la pendiente residuo."""
    lk = np.log(k[m])
    s_loc = local_slope(k[m], e_amp0[m]) + 5 / 3
    P = periodograma(lk, s_loc, OM_GRID)
    return float(P.max() / np.median(P)), float(OM_GRID[int(np.argmax(P))])


def espectro_observado(amp):
    return ME.E_log_periodico(K_LP, NU_LP, beta=ME.BETA_PAO,
                              amp=amp, omega=OMEGA_INY)


def espectro_nulo(rng):
    """
    Null: amp=0, mismos parametros Pao pero jitter instrumental leve:
    - beta ~ N(BETA_PAO, 0.5%)      (incertidumbre en constante del modelo)
    - nu   ~ lognormal(NU_LP, 1%)   (corre la ventana k*eta)
    - ruido multiplicativo lognormal sigma=1e-4 (ruido de medicion minimo:
      con un singleton determinista, cualquier pico es 100% artefacto
      de window+detrend; el jitter captura sensibilidad a esos errores).
    """
    beta = ME.BETA_PAO * (1.0 + 0.005 * rng.standard_normal())
    nu = NU_LP * np.exp(0.01 * rng.standard_normal())
    e = ME.E_pao(K_LP, nu, beta=beta)
    e = e * np.exp(1e-4 * rng.standard_normal(len(e)))
    return e


def main():
    n_null = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 20260922
    rng = np.random.default_rng(seed)
    m = (K_LP * ETA_LP > X_LO) & (K_LP * ETA_LP < X_HI)

    # 1. scores observados -------------------------------------------------
    obs = {}
    print("[obs] amplitudes observadas:")
    for amp in AMPS:
        s, w = score_de(espectro_observado(amp), K_LP, m)
        obs[amp] = (s, w)
        print(f"  amp={amp:6.3f}  SNR={s:10.2f}  omega_peak={w:.4f}")

    # 2. nulos -------------------------------------------------------------
    print(f"[null] generando {n_null} espectros nulos (seed={seed})...")
    null_scores = np.empty(n_null)
    null_omegas = np.empty(n_null)
    for i in range(n_null):
        null_scores[i], null_omegas[i] = score_de(espectro_nulo(rng), K_LP, m)
        if (i + 1) % 100 == 0:
            print(f"  {i+1}/{n_null}")

    # 3. p-values empiricos + umbrales ------------------------------------
    lam_global = float(null_scores.max())
    q95 = float(np.quantile(null_scores, 0.95))
    q99 = float(np.quantile(null_scores, 0.99))
    print(f"[null] score nulo: media={null_scores.mean():.2f} "
          f"max={lam_global:.2f} q95={q95:.2f} q99={q99:.2f}")

    rows = []
    print(f"\n{'amp':>6s} {'SNR_obs':>10s} {'p_empirico':>11s} "
          f"{'p_Bonf(5)':>10s} {'>q95?':>6s} {'>max?':>6s}")
    for amp in AMPS:
        s, w = obs[amp]
        p_emp = float((1 + np.sum(null_scores >= s)) / (n_null + 1))
        p_bonf = min(1.0, p_emp * len(AMPS))
        rows.append([amp, OMEGA_INY, f"{s:.4f}", f"{w:.4f}",
                     f"{p_emp:.5f}", f"{p_bonf:.5f}",
                     int(s > q95), int(s > lam_global)])
        print(f"{amp:6.3f} {s:10.2f} {p_emp:11.5f} {p_bonf:10.5f} "
              f"{str(s > q95):>6s} {str(s > lam_global):>6s}")

    with OUT.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["amp", "omega_iny", "snr_obs", "omega_peak",
                    "p_empirico", "p_bonferroni_5", "excede_q95",
                    "excede_null_max"])
        w.writerows(rows)
    print(f"\n[ok] {OUT}")


if __name__ == "__main__":
    main()
