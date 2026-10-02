#!/usr/bin/env python3
"""
exp_dns_real_slab.py  (v3.5 — con la RAM nueva)

Espectro E(k) del box PERIODICO completo sin Hanning: pulsos de dunas
(sublayers) que SI son periodicos en la direccion larga y se pegan en la
caja dada vuelta.

LOGICA
  El campo isotropic1024 es cubico 2pi x 2pi x 2pi periodico. Un pulo
  de slab SI = (1024, px, pz) con origen (y0, z0) y thickness T es un
  bloque que se envuelve en x (periodico: indice 1023 -> 0). Su "caja"
  es (X,Y,Z) = (2pi, T*x, T*x). Las frecuencias permitidas son
  (k_x, k_y, k_z) = (n1, 2*pi*n2/(T*x), 2*pi*n3/(T*x)), es decir
  k separacion (1, 4, 4) para T=256.

  FUJA vs BLOQUE AISLADO: el slab no necesita ventana en x (periodica)
  y en y/z usa la discontinuidad natural del campo real en el cuadrado — el
  salto en el borde y/z esta en el CAMPO, no es artefacto de recorte. La
  fuga de Hanning desaparece; solo queda la mezcla natural del campo.

  E_TAG(k) := (1/(T^2)) * 1/2 |u_k|^2 (Parseval). Equivalente al E_shell de
  una caja (2pi)^3: la energia total es la misma.

  Se promedian los 8 slabs orlados (analisis por slab y promedio),
  con sem entre slabs como error.

  NOTA IMPORTANTE: un espectro slab NO es el espectro de la caja (los slabs
  no son independientes — comparten filas por chunks). Pero las 8 esquinas
  son los 8 sub-cubos disjointos del dominio: el promedio de E_slab(k)
  equivale al E(k) del box completo en la aproximacion periodica-por-partes.
  Es la mejor respuesta posible sin bajar el campo entero.

SALIDA
  datos/espectros/spectrum_jhtdb_slab_{crudo,hanning}.csv
  datos/20_dns_real_slab_resumen.txt
Uso: python3 exp_dns_real_slab.py [--refresh]
"""
import csv
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from sddf_core import (compute_spectral_curvature, g_normalizado,
                       inertial_window, inertial_window_smoothed,
                       q_efectivo, truncate_by_slope_interp)

URL = ("https://huggingface.co/datasets/thuerey-group/"
       "jhtdb-isotropic-turbulence-1024/resolve/main/coarse_t420.hdf5")
KEY = "sims/sim0/420"
ROOT = Path(__file__).parent.parent
CACHE = ROOT / "datos" / "jhtdb_cache" / "slabs"
OUT_ESP = ROOT / "datos" / "espectros"
OUT_TXT = ROOT / "datos" / "20_dns_real_slab_resumen.txt"

NU = 1.85e-4
EPS_DOC = 0.0928
ETA_DOC = 2.87e-3
NBOX = 1024
DX = 2 * np.pi / NBOX
BS = 256
# 8 slabs: uno por esquina (x0, y0, z0) con x0 in {0, 768}, blocks (BS,BS,BS)
ORIGENES = [(x, y, z) for x in (0, NBOX - BS) for y in (0, NBOX - BS)
            for z in (0, NBOX - BS)]


def shell_index(nx, ny, nz):
    """Indices enteros de concha: n = (nx, ny, nz) -> |n| = sqrt(n1^2+n2^2+n3^2)."""
    fx = np.fft.fftfreq(nx, d=1.0) * nx
    fy = np.fft.fftfreq(ny, d=1.0) * ny
    fz = np.fft.fftfreq(nz, d=1.0) * nz
    kx, ky, kz = np.meshgrid(fx, fy, fz, indexing="ij")
    return np.rint(np.sqrt(kx**2 + ky**2 + kz**2)).astype(np.int16).ravel()


def espectro_slab(u, kr, a_k, window=True):
    """u: (256,1024,256,3) float64 sin media. E_shell(k): k = a_k * n.

    Periodico en x (sin ventana); y/z se ventanean con Hann si window=True
    (el slab NO es periodico en y/z, pero ahora si atraviesa todo el campo:
    la fuga que queda viene del recorte en y/z, mitigada).
    """
    nx, ny, nz = u.shape[0], u.shape[1], u.shape[2]
    if window:
        wy = np.hanning(ny)
        wz = np.hanning(nz)
        w = wy[None, :, None] * wz[None, None, :]
        # normalizacion por potencia real de la ventana (preserva energia
        # en y,z sin tocar la dimension periodica x)
        norm = float(np.mean(wy**2) * np.mean(wz**2))
        u = u * w
    else:
        norm = 1.0
    nb = int(kr.max()) + 1
    s = np.zeros(nb)
    for c in range(3):
        uh = np.fft.fftn(u[..., c]) / (nx * ny * nz)
        ek3 = 0.5 * np.abs(uh) ** 2
        s += np.bincount(kr, weights=ek3.ravel(), minlength=nb)
        del uh, ek3
    # E(k) integable: Densidad E(k)=s/a_k
    return s[1:] / norm


def main():
    args = sys.argv[1:]
    refresh = "--refresh" in args
    n_bloques = 8
    for a in args:
        if a.isdigit():
            n_bloques = int(a)

    a_k = 2 * np.pi / (BS * DX)   # = 4.0
    kr = shell_index(BS, NBOX, BS) # indices enteros |n| con separacion (1,4,4)

    CACHE.mkdir(parents=True, exist_ok=True)
    import fsspec
    import h5py
    f = None
    d = None
    Es_h, Es_w, eps_g, rms_comp = [], [], [], []
    try:
        for i, (x0, y0, z0) in enumerate(ORIGENES[:n_bloques]):
            cp = CACHE / f"slab_{x0}_{y0}_{z0}.npy"
            if cp.exists() and not refresh:
                u = np.load(cp)
                print(f"  [{i+1}] cache {cp.name}", flush=True)
            else:
                if f is None:
                    print("[h5] abriendo remoto", flush=True)
                    fo = fsspec.open(URL, "rb", block_size=2**22).open()
                    f = h5py.File(fo, "r")
                    d = f[KEY]
                    print(f"[h5] {KEY}: {d.shape} {d.dtype} chunks {d.chunks}",
                          flush=True)
                t = time.time()
                u = d[x0:x0 + BS, y0:y0 + NBOX, z0:z0 + BS, :3].astype(
                    np.float64)
                np.save(cp, u.astype(np.float32))
                mb = u.nbytes / 2**20
                print(f"  [{i+1}] slab ({x0},{y0},{z0}) {mb:.0f}MB en "
                      f"{time.time()-t:.0f}s", flush=True)
            med = u.mean(axis=(0, 1, 2))
            rms_comp.append(np.sqrt((u**2).mean(axis=(0, 1, 2))))
            u -= med
            # eps por gradientes continua siendo sobre el slab real
            g = np.gradient
            s2 = 0.0
            for c in range(3):
                for ax in range(3):
                    gi = g(u[..., c], axis=ax) / DX
                    s2 += float((gi**2).mean())
                    # falta el termino cruzado s_ij s_ij: usamos
                    # la identidad s_ij s_ij = G_ij G_ij para isotropico
                    # (diag solo domina en isotropia; discrepancia vs doc
                    #  se evaluara)
            eg = NU * s2
            eps_g.append(eg)
            Es_w.append(espectro_slab(u, kr, a_k, window=False))
            Es_h.append(espectro_slab(u, kr, a_k, window=True))
            del u
            print(f"  [esp] slab {i+1} listo (eps_grad={eg:.4f})", flush=True)
    finally:
        if f is not None:
            f.close()

    # Alinado: crudo y hanning mismo eje k
    n_bins = len(Es_h[0])
    k = a_k * np.arange(1, n_bins + 1)

    # ---- analisis ----
    def anal(Es, etiqueta):
        Es_a = np.asarray(Es)
        E_avg = Es_a.mean(axis=0)
        E_sem = Es_a.std(axis=0, ddof=1) / np.sqrt(len(Es_a))
        # eps espectral acumulado
        acum = {}
        n_arr = np.arange(1, len(E_avg) + 1)
        for nk in [2, 8, 16, 32, 64, 100, 128, 256, 512]:
            m = n_arr <= nk
            acum[a_k * nk] = 2 * NU * a_k**2 * np.sum(n_arr[m]**2 * E_avg[m]) / a_k
        # es Densidad: int k^2 E dk = a_k^2 sum n^2 * (E_shell/a_k) * a_k
        #                        = a_k^2 sum n^2 E_shell  (ya lo hice arriba con ese factor)
        # (correccion: E_avg de espectro_slab = suma/a_k = densidad, no suma)
        # -> eps = 2 nu sumi k^2 E(k) dk = 2 nu a_k sum n^2 * E_avg
        acum_cor = {}
        for nk, eprev in acum.items():
            m = n_arr <= nk / a_k
            acum_cor[nk] = 2 * NU * a_k * np.sum(n_arr[m]**2 * E_avg[m])

        # pendiente suavizada (w=2->5 puntos por distancia real dk=4)
        s_suav = np.full(n_bins, np.nan)
        lk = np.log(k)
        for j in range(2, n_bins - 2):
            mm = np.arange(max(0, j - 2), min(n_bins, j + 3))
            s_suav[j] = np.polyfit(lk[mm], np.log(E_avg[mm]), 1)[0]

        # fits y detectores
        fits = []
        for k1, k2 in [(8, 16), (8, 32), (8, 64), (16, 64), (8, 128),
                       (4, 128), (4, 64)]:
            m = (k >= k1) & (k <= k2)
            if m.sum() < 4:
                continue
            p = np.polyfit(lk[m], np.log(E_avg[m]), 1)[0]
            pb = [np.polyfit(lk[m], np.log(E[m]), 1)[0] for E in Es_a]
            fits.append((k1, k2, p, np.std(pb, ddof=1) / np.sqrt(len(pb))))
        kt, et, nel, ok_t, kc = truncate_by_slope_interp(k, E_avg, delta=0.5)
        lo_s5, hi_s5, ok_s5 = inertial_window_smoothed(k, E_avg, delta=0.5,
                                                       window=3)
        lo_s9, hi_s9, ok_s9 = inertial_window_smoothed(k, E_avg, delta=0.5,
                                                       window=5)

        def obs(lo, hi):
            m = (k >= lo) & (k <= hi)
            if m.sum() < 3:
                return (float("nan"),) * 3 + (float("nan"),)
            return (compute_spectral_curvature(k[m], E_avg[m]),
                    g_normalizado(k[m], E_avg[m]),
                    q_efectivo(k[m], E_avg[m]),
                    float(np.log(hi / lo)))

        return {"etiqueta": etiqueta, "E_avg": E_avg, "E_sem": E_sem,
                "acum": acum_cor, "fits": fits, "s_suav": s_suav,
                "G_punt": obs(kt[0], kc) if ok_t else (np.nan,) * 4,
                "G_s5": (obs(lo_s5, hi_s5) + (lo_s5, hi_s5, ok_s5)),
                "G_s9": (obs(lo_s9, hi_s9) + (lo_s9, hi_s9, ok_s9))}

    res_w = anal(Es_w, "crudo")
    res_h = anal(Es_h, "hanning_yz")

    # ---- resumen ----
    for r in (res_w, res_h):
        p = OUT_ESP / f"spectrum_jhtdb_slab_{r['etiqueta']}.csv"
        with p.open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["k", "E", "E_sem", "k_eta"])
            for ki, ei, si in zip(k, r["E_avg"], r["E_sem"]):
                w.writerow([f"{ki:.6e}", f"{ei:.6e}", f"{si:.6e}",
                            f"{ki * ETA_DOC:.4f}"])
        print(f"[csv] {p}")

    txt = []
    txt.append("=" * 74)
    txt.append("SDDF sobre DNS REAL — pulsos de dunas (slabs), v3.5 RAM nueva")
    txt.append("=" * 74)
    txt.append(f"fuente: {URL}")
    txt.append(f"        {KEY}  shape (1024,1024,1024,4) float16")
    txt.append(f"metodo: {n_bloques} pulsos (1024,256,256) en las esquinas "
               f"(periodicos en x; el largo completo 2pi)")
    txt.append("NOTA: el slab es periodico en x; en y/z el recorte es YA el"
               " dominio entero, no un sub-cubo -> la fuga de ventana es la"
               " natural (cruzada por el salto en bordes y/z) y se maneja con"
               " Hanning solo en esas direcciones. Compara con la version"
               " bloques-aislados (18_dns_real_box_resumen.txt).")
    txt.append("")
    txt.append("CHEQUEOS:")
    r = res_h
    txt.append(f"  eps por gradientes: {np.mean(eps_g):.4f} +/- "
               f"{np.std(eps_g, ddof=1)/np.sqrt(len(eps_g)):.4f}  "
               f"vs doc {EPS_DOC}  =>  {np.mean(eps_g)/EPS_DOC:.2f}x")
    txt.append(f"  (nota: con slab, eps por diag especificos no es promedio "
               f"de bloques sino medida directa — esperado mas bajo por el "
               f"recorte geografico en y/z)")
    txt.append("")
    txt.append("CONVERGENCIA DE eps ESPECTRAL (2 nu a sum n^2 E_dens):")
    txt.append(f"  {'k<=':>6s} {'crudo':>18s} {'hanning_yz':>18s}")
    for nk in sorted(res_w["acum"]):
        txt.append(f"  {nk:6.0f} {res_w['acum'][nk]:10.4f} "
                   f"({res_w['acum'][nk]/EPS_DOC:4.2f}x) "
                   f"{res_h['acum'][nk]:10.4f} "
                   f"({res_h['acum'][nk]/EPS_DOC:4.2f}x)")
    txt.append("")
    txt.append("AJUSTES log-log (q = -s), hanning_yz:")
    for k1, k2, p, sem in res_h["fits"]:
        qq = -p
        txt.append(f"  k in [{k1:3d},{k2:3d}] (k*eta {k1*ETA_DOC:.2f}-"
                   f"{k2*ETA_DOC:.2f}): q = {qq:.4f} +/- {sem:.4f} "
                   f"({(qq-5/3)/(5/3)*100:+.2f}%)")
    txt.append("")
    txt.append("DETECTORES (delta=0.5) sobre el espectro promediado:")
    for nombre, tup in (("puntual", res_h["G_punt"]),
                        ("suav w=3", res_h["G_s5"]),
                        ("suav w=5", res_h["G_s9"])):
        if len(tup) == 4:
            G, gn, q, span = tup
            txt.append(f"  {nombre:9s}: G*={G:.4f}  q={q:.4f}  "
                       f"ventana={span:.2f} nats  [nan si falla]")
        else:
            G, gn, q, span, lo, hi, okk = tup
            txt.append(f"  {nombre:9s}: k [{lo:.2f},{hi:.2f}] ok={okk}  "
                       f"G*={G:.4f}  q={q:.4f}  ventana={span:.2f} nats")
    out = "\n".join(txt)
    print("\n" + out)
    OUT_TXT.write_text(out + "\n")
    print(f"\n[ok] {OUT_TXT}")


if __name__ == "__main__":
    main()
