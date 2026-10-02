"""Espectro del box entero desde cache local: promedio de los 8 sub-bloques.

Es el mismo metodo validado de exp_dns_real_fullbox.py (8 bloques octantes,
Hanning, w=5), pero con la RAM nueva puedo guardar los bloques y pasar el
promedio una sola vez, ademas de hacer el test de A/B con ventana y sin.
"""
import csv
import sys
import time
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, "/home/delorien/sddf/codigo")
from sddf_core import (compute_spectral_curvature, g_normalizado,
                       inertial_window_smoothed, q_efectivo,
                       smoothed_slope)

ROOT = Path("/home/delorien/sddf")
H5 = ROOT / "datos" / "jhtdb_cache" / "coarse_t420.hdf5"
OUT_ESP = ROOT / "datos" / "espectros"
OUT_TXT = ROOT / "datos" / "20_dns_real_box_local_resumen.txt"
CACHE_OUT = ROOT / "datos" / "jhtdb_cache" / "blocks"

NU = 1.85e-4
EPS_DOC = 0.0928
ETA_DOC = 2.87e-3
DX = 2 * np.pi / 1024
BS = 256
A_K = 2 * np.pi / (BS * DX)

fx = np.fft.fftfreq(BS) * BS
KX, KY, KZ = np.meshgrid(fx, fx, fx, indexing="ij")
kr = np.rint(np.sqrt(KX**2 + KY**2 + KZ**2)).astype(np.int16).ravel()
nb = int(kr.max()) + 1

def espectro(u, kr, a_k, window=True):
    n = u.shape[0]
    if window:
        w1 = np.hanning(n)
        u = u * w1[:, None, None]
        norm = float(np.mean(w1**2))
    else:
        norm = 1.0
    s = np.zeros(nb)
    for c in range(3):
        uh = np.fft.fftn(u[..., c]) / n**3
        ek = 0.5 * np.abs(uh) ** 2
        s += np.bincount(kr, weights=ek.ravel(), minlength=nb)
        del uh, ek
    return s[1:] / norm

def eps_gradientes(u, nu=NU, dx=DX):
    g = {}
    for c in range(3):
        for ax in range(3):
            g[(c, ax)] = np.gradient(u[..., c], axis=ax) / dx
    s2 = 0.0
    for i in range(3):
        for j in range(3):
            d = 0.5 * (g[(i, j)] + g[(j, i)])
            s2 += float((d**2).mean())
    return 2 * nu * s2

origenes = [(x, y, z) for x in (0, 768) for y in (0, 768)
            for z in (0, 768)]

print("cargando directamente desde el H5 local (mas rapido que vitro los "
      "bloques cacheados)...", flush=True)
t0 = time.time()
Es_h, Es_c, eps_g, rms_comp = [], [], [], []
with h5py.File(H5, "r") as f:
    dset = f["sims/sim0/420"]
    for i, (x0, y0, z0) in enumerate(origenes):
        u = dset[x0:x0 + BS, y0:y0 + BS, z0:z0 + BS, :3].astype(np.float64)
        m = u.mean(axis=(0, 1, 2))
        rms_comp.append(np.sqrt((u**2).mean(axis=(0, 1, 2))))
        u = u - m
        eps_g.append(eps_gradientes(u))
        Es_h.append(espectro(u, kr, A_K, window=True))
        Es_c.append(espectro(u, kr, A_K, window=False))
        del u
        print(f"  [{i+1}/8] ({x0},{y0},{z0}) eps_g={eps_g[-1]:.4f} "
              f"{time.time()-t0:.0f}s", flush=True)

Es_h = np.array(Es_h)
Es_c = np.array(Es_c)
E_avg_h = Es_h.mean(axis=0)
E_sem_h = Es_h.std(axis=0, ddof=1) / np.sqrt(8)
E_avg_c = Es_c.mean(axis=0)
k = A_K * np.arange(1, nb)

# eps espectral acumulado (definicion: int k^2 E(k) dk, E = suma shell
# y dk = a_k es el paso en k del grid slab (1,4,1) -> nshell = a_k * n)
n_arr = np.arange(1, nb)
eps_sp = float(2 * NU * np.sum(k**2 * E_avg_h) * A_K)
print(f"\n[eps] gradientes: {np.mean(eps_g):.4f} +/- "
      f"{np.std(eps_g, ddof=1)/np.sqrt(8):.4f} vs doc {EPS_DOC} "
      f"({np.mean(eps_g)/EPS_DOC:.2f}x)")
print(f"[eps] espectral (sum k^2 E * dk, todo el rango): {eps_sp:.4f} "
      f"({eps_sp/EPS_DOC:.2f}x doc)")

# out CSV
OUT_ESP.mkdir(parents=True, exist_ok=True)
for tag, Ev, Es in (("hanning", E_avg_h, E_sem_h), ("crudo", E_avg_c, None)):
    if Es is None:
        Es = Es_c.std(axis=0, ddof=1) / np.sqrt(8)
    with (OUT_ESP / f"spectrum_jhtdb_box_{tag}_v35.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["k", "E", "E_sem", "k_eta"])
        for a, b, c in zip(k, Ev, Es):
            w.writerow([f"{a:.6e}", f"{b:.6e}", f"{c:.6e}", f"{a*ETA_DOC:.4f}"])

# SDDF
for tag, Ev in (("crudo", E_avg_c), ("hanning", E_avg_h)):
    s_s = smoothed_slope(k, Ev, window=5)
    lo, hi, ok = inertial_window_smoothed(k, Ev, delta=0.5, window=5)
    m = (k >= lo) & (k <= hi)
    if ok and m.sum() >= 4:
        G = compute_spectral_curvature(k[m], Ev[m])
        q = q_efectivo(k[m], Ev[m])
        print(f"[{tag}] k in [{lo:.1f},{hi:.1f}] G*={G:.4f} q={q:.4f}")
    else:
        print(f"[{tag}] detector suavizado: ok={ok}")

txt = []
txt.append("SDDF v3.5 — box periódico completo (cache local, 8 bloques)")
txt.append("=" * 60)
txt.append(f"eps gradientes = {np.mean(eps_g):.4f} ± "
           f"{np.std(eps_g, ddof=1)/np.sqrt(8):.4f} vs doc 0.0928 "
           f"({np.mean(eps_g)/EPS_DOC:.2f}x)")
txt.append(f"eps espectral = {eps_sp:.4f} ({eps_sp/EPS_DOC:.2f}x doc)")
txt.append("SDDF sobre ventana declarada (delta=0.5, w=5): ver tabla arriba")
OUT_TXT.write_text("\n".join(txt) + "\n")
print(f"\n[ok] {OUT_TXT}")
