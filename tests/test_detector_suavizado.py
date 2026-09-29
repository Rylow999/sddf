#!/usr/bin/env python3
"""Tests de regresion para el detector suavizado (v3.4).

Leccion del contraste con DNS real: el detector puntual (inertial_window)
se dispara con el ruido concha-a-concha del espectro real; el suavizado
(inertial_window_smoothed) encuentra el rango correcto.
"""
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent / "codigo"))
from sddf_core import (compute_spectral_curvature, g_normalizado,
                       inertial_window, inertial_window_smoothed,
                       q_efectivo, smoothed_slope)

try:
    import modelos_espectrales as ME
except ImportError:
    ME = None


def test_smoothed_slope_k41():
    """K41 puro: la pendiente suavizada tiene que ser -5/3 en el interior."""
    k = np.logspace(0, 3, 200)
    E = k ** (-5 / 3)
    s = smoothed_slope(k, E, window=5)
    interior = s[5:-5]
    assert np.allclose(interior, -5 / 3, atol=1e-9), \
        f"K41: s suavizada = {interior[:3]}... debe ser -5/3"
    print("[OK] smoothed_slope: K41 puro reproduce -5/3 exacto en interior")


def test_smoothed_vs_puntual_sintetico_limpio():
    """Sin ruido, suavizado y puntual encuentran la misma ventana."""
    if ME is None:
        print("[SKIP] modelos_espectrales no disponible")
        return
    k = np.logspace(0, 3, 400)
    E = ME.E_pao(k, 1e-4, beta=ME.BETA_PAO)
    lo_s, hi_s, ok_s = inertial_window_smoothed(k, E, delta=0.5, window=5)
    lo_p, hi_p, ok_p = inertial_window(k, E, delta=0.5)
    assert ok_s and ok_p, "sin ruido ambos detectores deben encontrar ventana"
    assert abs(hi_s - hi_p) / hi_p < 0.05, \
        f"k_high difiere: suav={hi_s:.2f} vs puntual={hi_p:.2f}"
    print(f"[OK] sintetico limpio: suavizado k_high={hi_s:.2f} ~ "
          f"puntual {hi_p:.2f}")


def test_puntual_falla_con_ruido_suavizado_no():
    """Con ruido 1%: el puntual pierde la ventana, el suavizado w=9 no."""
    if ME is None:
        print("[SKIP]")
        return
    rng = np.random.default_rng(7)
    k = np.logspace(0, 3, 400)
    E = ME.E_pao(k, 1e-4, beta=ME.BETA_PAO) * np.exp(
        0.01 * rng.standard_normal(len(k)))
    _, _, ok_p = inertial_window(k, E, delta=0.5)
    lo_s, hi_s, ok_s = inertial_window_smoothed(k, E, delta=0.5, window=9)
    assert not ok_p, "el puntual deberia fallar con ruido 1% (sem/slope~0.4)"
    assert ok_s, "el suavizado w=9 deberia encontrar la ventana con ruido 1%"
    span = np.log(hi_s / lo_s)
    assert span > 3, f"ventana suavizada demasiado corta: {span:.2f} nats"
    print(f"[OK] ruido 1%: puntual falla (ok=0), suavizado w=9 encuentra "
          f"k in [{lo_s:.2f},{hi_s:.2f}] ({span:.2f} nats)")


def test_dns_real():
    """DNS real (8 bloques, hanning): puntual falla, suavizado w=5 encuentra.

    El detector es conservador en el borde bajo por diseño (excluye los
    primeros window//2 puntos): k_low >= 12 para el grid del box.
    """
    csv_path = (Path(__file__).parent.parent / "datos" / "espectros" /
                "spectrum_jhtdb_box_hanning.csv")
    if not csv_path.exists():
        print("[SKIP] spectrum_jhtdb_box_hanning.csv no existe (corre "
              "exp_dns_real_fullbox.py primero)")
        return
    rows = list(csv.reader(csv_path.open()))
    k = np.array([float(r[0]) for r in rows[1:]])
    E = np.array([float(r[1]) for r in rows[1:]])
    _, _, ok_p = inertial_window(k, E, delta=0.5)
    lo_s, hi_s, ok_s = inertial_window_smoothed(k, E, delta=0.5, window=5)
    assert not ok_p, "el puntual deberia fallar en el DNS real (v3.3)"
    assert ok_s, "el suavizado w=5 deberia encontrar ventana en DNS real"
    m = (k >= lo_s) & (k <= hi_s)
    q = q_efectivo(k[m], E[m])
    assert 1.3 < q < 1.8, f"q fuera del rango fisico esperado: {q:.4f}"
    print(f"[OK] DNS real: puntual falla, suavizado w=5 encuentra "
          f"k in [{lo_s:.1f},{hi_s:.1f}], q={q:.4f} (K41 1.667)")


def main():
    test_smoothed_slope_k41()
    test_smoothed_vs_puntual_sintetico_limpio()
    test_puntual_falla_con_ruido_suavizado_no()
    test_dns_real()
    print("\nTodos los tests pasaron.")


if __name__ == "__main__":
    main()
