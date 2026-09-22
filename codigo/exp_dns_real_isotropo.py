#!/usr/bin/env python3
"""
exp_dns_real_isotropo.py  (pendiente #3, README v3.1)

Validacion de G*[u] con DNS REAL: JHTDB isotropic1024coarse (256^3).

Descarga el snapshot de velocidad del mirror publico en HuggingFace
(ArielLubonja/johns-hopkins-turbulence-database, un solo timestep),
calcula el espectro de energia isotropo E(k) por FFT + bineado en
conchas de |k|, y corre el pipeline SDDF (corte por pendiente con
interpolacion y medidor del exponente corregido) sobre datos reales.

Salida:
  datos/espectros/spectrum_jhtdb_isotropic1024coarse_t{t}.csv
  datos/17_dns_real_resumen.txt

Uso: python3 exp_dns_real_isotropo.py [timestep 0-9]
"""
import sys
import urllib.request
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import sddf_exact as EX
from sddf_core import (compute_spectral_curvature, g_normalizado,
                       truncate_by_slope_interp, q_efectivo)

try:
    import h5py
except ImportError:
    sys.exit("Falta h5py: .venv/bin/pip install h5py")

HF_URL = ("https://huggingface.co/datasets/ArielLubonja/"
          "johns-hopkins-turbulence-database/resolve/main/"
          "isotropic1024-coarse-velocity.h5")
CACHE = Path(__file__).parent.parent / "datos" / "jhtdb_cache"
OUT_ESP = Path(__file__).parent.parent / "datos" / "espectros"
OUT_TXT = Path(__file__).parent.parent / "datos" / "17_dns_real_resumen.txt"

# Parametros oficiales del dataset isotropic1024coarse (FORCED, Re_l ~ 433)
NU = 1.85e-4            # viscosidad cinematica
DT_STORED = 0.02        # separacion entre snapshots almacenados
NX = 256                # lado de la grilla (coarse)


def descargar(dest: Path):
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        return
    print(f"[dl] bajando {HF_URL} (~2 GB) ...")
    urllib.request.urlretrieve(HF_URL, dest)
    print(f"[dl] ok -> {dest}")


def espectro_isotropo(u: np.ndarray):
    """E(k) 1D isotropo por FFT 3D + bineado en conchas |k|=entero."""
    n = u.shape[0]
    uh = np.fft.fftn(u, axes=(0, 1, 2)) / n**3       # hat(u) (parseval-safe)
    ek3d = 0.5 * np.sum(np.abs(uh)**2, axis=3)       # (n,n,n)
    fx = np.fft.fftfreq(n) * n
    kx, ky, kz = np.meshgrid(fx, fx, fx, indexing="ij")
    kr = np.sqrt(kx**2 + ky**2 + kz**2).astype(int)
    kmax = n // 2
    E = np.zeros(kmax)
    cnt = np.bincount(kr.ravel(), minlength=kmax + 1)
    suma = np.bincount(kr.ravel(), weights=ek3d.ravel(), minlength=kmax + 1)
    k_int = np.arange(0, kmax)
    # densidad de modos ~ 4pi k^2: E(k) = suma_concha / N_realizaciones  (ya
    # incluido en ek3d); para el espectro que DEFINEN en DNS: E(k) = suma
    # sobre la concha de 1/2|u_hat|^2 (sin dividir), lo que da la forma
    # clasica E(k) ~ k^2 en k->0 y -5/3 en inercial.
    Eshell = suma[:kmax]
    return k_int[1:], Eshell[1:]


def main():
    tstep = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    f5 = CACHE / "isotropic1024-coarse-velocity.h5"
    descargar(f5)
    with h5py.File(f5, "r") as f:
        print("[h5] keys:", list(f.keys()))
        dset = None
        for name in f.keys():
            d = f[name]
            if isinstance(d, h5py.Dataset) and d.ndim >= 4:
                dset = d
                print(f"[h5] usando '{name}' shape={d.shape}")
                break
        # shape esperado: (256,256,256,3,10)
        u = dset[..., tstep] if dset.ndim == 5 else dset[..., 0]
    print(f"[t] timestep={tstep}, u.shape={u.shape}")

    k, E = espectro_isotropo(u)
    OUT_ESP.mkdir(parents=True, exist_ok=True)
    csv_path = OUT_ESP / f"spectrum_jhtdb_isotropic1024coarse_t{tstep}.csv"
    with csv_path.open("w") as fcsv:
        fcsv.write("k,E\n")
        for ki, ei in zip(k, E):
            fcsv.write(f"{ki:.6e},{ei:.6e}\n")
    print(f"[csv] {csv_path}")

    # --- pipeline SDDF sobre datos reales --------------------------------
    eps_oficial = 0.0928  # disipacion media oficial isotropic1024coarse (Re_l=433)
    eta = (NU**3 / eps_oficial) ** 0.25
    kt, et, nel, ok, kc = truncate_by_slope_interp(k, E, delta=0.5)
    gn = g_normalizado(kt, et)
    q = q_efectivo(kt, et)
    G = compute_spectral_curvature(kt, et)
    span = np.log(kt[-1] / kt[0])
    g_pred = (25 / 12) * np.log(1 / (eta * kt[0]))  # solo referencia relativa

    txt = []
    txt.append(f"JHTDB isotropic1024coarse  snapshot t={tstep} (256^3, Re_l~433)")
    txt.append(f"nu = {NU}, eps = 0.0928 (oficial), eta = {eta:.4e} (grid units)")
    txt.append(f"corte (delta=0.5): k en [{kt[0]:.1f}, {kt[-1]:.1f}]  "
               f"({len(kt)} pts, interp ok={ok})")
    txt.append(f"span inercial = {span:.3f} nats")
    txt.append(f"G* = {G:.4f}")
    txt.append(f"<s^2> normalizado = {gn:.4f}  (K41 puro: 25/9 = 2.7778)")
    txt.append(f"q efectivo = sqrt(<s^2>) = {q:.4f}  (K41: 5/3 = 1.6667)")
    txt.append(f"sesgo vs K41: {(q - 5/3) / (5/3) * 100:+.2f} %")
    out = "\n".join(txt)
    print("\n" + out)
    OUT_TXT.write_text(out + "\n")
    print(f"\n[ok] {OUT_TXT}")


if __name__ == "__main__":
    main()
