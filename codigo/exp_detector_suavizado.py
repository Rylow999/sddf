#!/usr/bin/env python3
"""
exp_detector_suavizado.py  (v3.4, leccion metodologica del DNS real)

Validacion cruzada del detector de rango inercial sobre pendiente
SUAVIZADA (inertial_window_smoothed) contra el detector puntual
(inertial_window):

  A. SINTETICO (donde el puntual funciona): espectros Pao/quimera con
     ruido multiplicativo creciente (0%, 1%, 5%, 3%...) para ver cuando
     cada detector empieza a fallar.
  B. DNS REAL (donde el puntual falla): el espectro promediado de 8
     bloques del box 1024^3 (spectrum_jhtdb_box_hanning.csv).

Salida: datos/19_detector_suavizado.csv + resumen.
"""
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import modelos_espectrales as ME
from sddf_core import (compute_spectral_curvature, g_normalizado,
                       inertial_window, inertial_window_smoothed,
                       q_efectivo, smoothed_slope)

ROOT = Path(__file__).parent.parent
DNS_CSV = ROOT / "datos" / "espectros" / "spectrum_jhtdb_box_hanning.csv"
OUT = ROOT / "datos" / "19_detector_suavizado.csv"

ETA_DOC = 2.87e-3
NU = 1.85e-4
EPS_DOC = 0.0928


def evaluar(k, E, etiqueta, semilla=None):
    filas = []
    for w in [3, 5, 9, 11]:
        lo_s, hi_s, ok_s = inertial_window_smoothed(k, E, delta=0.5,
                                                    window=w)
        m = (k >= lo_s) & (k <= hi_s) if ok_s else None
        if m is not None and m.sum() >= 4:
            G = compute_spectral_curvature(k[m], E[m])
            q = q_efectivo(k[m], E[m])
            span = np.log(hi_s / lo_s)
        else:
            G = q = span = np.nan
        filas.append([etiqueta, w, lo_s, hi_s, int(ok_s), span, G, q])
    # puntual para comparar
    lo_p, hi_p, ok_p = inertial_window(k, E, delta=0.5)
    m = (k >= lo_p) & (k <= hi_p) if ok_p else None
    if m is not None and m.sum() >= 4:
        G = compute_spectral_curvature(k[m], E[m])
        q = q_efectivo(k[m], E[m])
        span = np.log(hi_p / lo_p)
    else:
        G = q = span = np.nan
    filas.append([etiqueta, "puntual", lo_p, hi_p, int(ok_p), span, G, q])
    return filas


def parte_sintetica():
    print("=" * 74)
    print("A. SINTETICO: espectro quimera + ruido multiplicativo creciente")
    print("=" * 74)
    rng = np.random.default_rng(42)
    k = np.logspace(0, 3, 400)          # 3 decadas, como un DNS mediano
    eta = (NU**3 / EPS_DOC) ** 0.25
    E0 = ME.E_pao(k, NU, beta=ME.BETA_PAO)

    filas = []
    print(f"{'ruido':>7s} {'metodo':>9s} {'k_low':>8s} {'k_high':>9s} "
          f"{'ok':>4s} {'span':>6s} {'G*':>9s} {'q':>7s}")
    for noise in [0.0, 0.01, 0.03, 0.05, 0.10]:
        E = E0 * np.exp(noise * rng.standard_normal(len(k)))
        filas += evaluar(k, E, f"sint_{noise:.2f}")
        for f in filas[-5:]:
            print(f"{f[0][5:]:>7s} {str(f[1]):>9s} {f[2]:8.2f} {f[3]:9.2f} "
                  f"{f[4]:4d} {f[5]:6.2f} {f[6]:9.3f} {f[7]:7.4f}")
    return filas


def parte_dns():
    print("\n" + "=" * 74)
    print("B. DNS REAL: espectro promedio 8 bloques (box 1024^3, hanning)")
    print("=" * 74)
    rows = list(csv.reader(DNS_CSV.open()))
    k = np.array([float(r[0]) for r in rows[1:]])
    E = np.array([float(r[1]) for r in rows[1:]])
    filas = evaluar(k, E, "dns_real")
    print(f"{'metodo':>9s} {'k_low':>8s} {'k_high':>9s} {'ok':>4s} "
          f"{'span':>6s} {'G*':>9s} {'q':>7s}")
    for f in filas:
        print(f"{str(f[1]):>9s} {f[2]:8.2f} {f[3]:9.2f} {f[4]:4d} "
              f"{f[5]:6.2f} {f[6]:9.3f} {f[7]:7.4f}")
    return filas


def main():
    filas = parte_sintetica() + parte_dns()
    with OUT.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["caso", "metodo", "k_low", "k_high", "ok",
                    "span_nats", "G_estrella", "q_efectivo"])
        w.writerows(filas)
    print(f"\n[ok] {OUT}")


if __name__ == "__main__":
    main()
