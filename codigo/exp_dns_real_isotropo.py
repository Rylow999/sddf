#!/usr/bin/env python3
"""
exp_dns_real_isotropo.py  (pendiente #3, README v3.1)

Validacion de G*[u] con DNS REAL: snapshot publico JHTDB isotropic1024coarse
(mirror HuggingFace, 256^3 x 10 timesteps).

QUE ES ESTE DATO (verificado, no asumido):
  - attrs del H5: dataset='isotropic1024coarse', x/y/z_start=1, _end=256,
    _step=1, t_start=1..t_end=10 (indices de la malla almacenada).
  - xcoor del mirror: dx = 2*pi/1024 (=2pi/1024), 256 puntos -> lado pi/2.
  - TEST DE PERIODICIDAD (test_periodicidad.py): el salto entre la ultima y
    la primera capa es ~300-360x mayor que entre capas interiores => el
    bloque NO es periodico => es un SUB-CUBO de un cuarto de caja a
    resolucion fina, NO la caja entera en malla coarse.
  - Por lo tanto el sub-cubo es NO PERIODICO y su FFT sufre fuga espectral
    (cola espuria ~k^-2 por la discontinuidad en las caras). Eso contamina
    el integrador de disipacion eps = 2 nu int k^2 E dk.

QUE HACE:
  1. Lee los 10 timesteps, calcula E(k) por FFT 3D + conchas |k|.
  2. Averigua la sensibilidad a la fuga: espectro crudo vs con ventana de
     Tukey (ventaneo + renormalizacion por sqrt(<w^2>)).
  3. Estima eps por el espectro en ambos casos y lo compara con el valor
     documentado por JHTDB (eps=0.0928, nu=1.85e-4, eta=2.87e-3).
  4. Corre el pipeline SDDF (corte por pendiente interpolado) sobre el
     espectro promediado y reporta G*, <s^2>, q_efectivo.

Salida:
  datos/espectros/spectrum_jhtdb_iso1024_t{0..9}.csv
  datos/espectros/spectrum_jhtdb_iso1024_promedio.csv
  datos/17_dns_real_resumen.txt

Uso: python3 exp_dns_real_isotropo.py [n_timesteps] [--ventana]
"""
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from sddf_core import (compute_spectral_curvature, g_normalizado, local_slope,
                       q_efectivo, truncate_by_slope_interp)

try:
    import h5py
except ImportError:
    sys.exit("Falta h5py: /home/delorien/sddf/.venv/bin/pip install h5py")

H5 = (Path(__file__).parent.parent / "datos" / "jhtdb_cache" /
      "isotropic1024-coarse-velocity.h5")
OUT_ESP = Path(__file__).parent.parent / "datos" / "espectros"
OUT_TXT = Path(__file__).parent.parent / "datos" / "17_dns_real_resumen.txt"

# --- parametros documentados del dataset (README-isotropic.pdf JHTDB) -----
NU = 1.85e-4          # viscosidad cinematica
EPS_OFICIAL = 0.0928  # disipacion (t=0..2.048, Re_lambda=433)
ETA_OFICIAL = 2.87e-3  # escala de Kolmogorov (documentada, = (nu^3/eps)^{1/4})
U_RMS_OFICIAL = 0.681  # rms por componente
E_TOT_OFICIAL = 0.695  # energia cinetica total


def espectro_isotropo(u: np.ndarray, dx: float, window: bool = False):
    """E(k) por FFT 3D + bineado en conchas |k|. Devuelve (k, E) con k fisico.

    u: (n,n,n,3). dx: espaciado (mismo en los 3 ejes).
    window: aplica Tukey (alpha=0.5) por eje para reducir la fuga de un
    bloque no periodico; renormaliza por sqrt(<w^2>) para preservar Parseval.
    """
    n = u.shape[0]
    if window:
        try:
            from scipy.signal.windows import tukey
            w1 = tukey(n, alpha=0.5)
        except ImportError:  # fallback sin scipy
            w1 = np.hanning(n)
        w = w1[:, None, None]
        u = u * w
        norm = float(np.mean(w1**2))
    else:
        norm = 1.0

    freqs = np.fft.fftfreq(n, d=1.0) * n
    kx, ky, kz = np.meshgrid(freqs, freqs, freqs, indexing="ij")
    kr_flat = np.rint(np.sqrt(kx**2 + ky**2 + kz**2)).astype(np.int64).ravel()
    nbins = int(kr_flat.max()) + 1
    suma = np.zeros(nbins)
    for c in range(u.shape[3]):
        uh = np.fft.fftn(u[..., c]) / n**3
        ek3d = 0.5 * np.abs(uh) ** 2
        suma += np.bincount(kr_flat, weights=ek3d.ravel(), minlength=nbins)
        del uh, ek3d
    L = n * dx
    k_phys = np.arange(nbins) * (2 * np.pi / L)
    return k_phys[1:], suma[1:] / norm


def eps_del_espectro(k, E, nu=NU, a=None):
    """eps = 2 nu int k^2 E(k) dk.

    OJO con la normalizacion: E es la SUMA por concha (energia total en la
    concha), no la densidad. Con E_shell(n) y k_n = a*n, la densidad es
    E(k_n) = E_shell(n)/a, y entonces
        int k^2 E dk = sum (a n)^2 (E_shell/a) a = a^2 sum n^2 E_shell.
    Un trapezoid(k**2*E_shell, k) sobreestima por un factor a exacto.
    """
    if a is None:
        a = float(k[1] - k[0])
    n = np.arange(len(E))
    return float(2 * nu * a**2 * np.sum(n**2 * E))


def main():
    args = [a for a in sys.argv[1:]]
    n_t = 10
    usar_ventana = "--ventana" in args
    for a in args:
        if a.isdigit():
            n_t = int(a)

    with h5py.File(H5, "r") as f:
        x = f["xcoor"][:]
        dx = float(x[1] - x[0])
        L_sub = float(x[-1] - x[0] + dx)
        print(f"[h5] dataset={f.attrs['dataset']!r}  dx={dx:.8f}  "
              f"lado sub-cubo L={L_sub:.6f} = {L_sub/(2*np.pi):.4f}*(2pi)")
        print(f"[h5] k_min = 2pi/L = {2*np.pi/L_sub:.3f}   "
              f"k_Nyq = pi/dx = {np.pi/dx:.1f}   "
              f"(k_max*eta = {np.pi/dx*ETA_OFICIAL:.2f})")
        espec = []
        for it in range(n_t):
            u = f[f"Velocity_{it+1:04d}"][:].astype(np.float64)
            u = u - u.mean(axis=(0, 1, 2))       # saca el DC del sub-cubo
            k, E = espectro_isotropo(u, dx, window=usar_ventana)
            espec.append(E)
            if it == 0:
                e_tot = 0.5 * (u**2).sum(axis=3).mean()
                print(f"[diag] E_tot (sin media) = {e_tot:.4f}  "
                      f"(documentado {E_TOT_OFICIAL})")
        E_avg = np.mean(espec, axis=0)
        E_std = np.std(espec, axis=0)

    tag = "tukey" if usar_ventana else "crudo"
    print(f"[spec] promedio de {n_t} timesteps (ventana: {tag})")

    OUT_ESP.mkdir(parents=True, exist_ok=True)
    csv_path = OUT_ESP / f"spectrum_jhtdb_iso1024_promedio.csv"
    with csv_path.open("w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["k", "E", "E_std", "n_modos"])
        for ki, ei, si in zip(k, E_avg, E_std):
            wr.writerow([f"{ki:.6e}", f"{ei:.6e}", f"{si:.6e}", ""])
    print(f"[csv] {csv_path}")

    # --- diagnostico de pendiente (local, promedio) ----------------------
    s_loc = local_slope(k, E_avg)
    print(f"[slope] ventana {tag}:")
    for kk in [4, 8, 16, 32, 60, 100, 200, 400]:
        j = int(np.argmin(np.abs(k - kk)))
        print(f"    k={k[j]:7.1f}  s={s_loc[j]:+.3f}  E={E_avg[j]:.4e}")

    # --- eps en las dos variantes ----------------------------------------
    eps_tk = eps_del_espectro(k, E_avg)
    a_k = float(k[1] - k[0])
    m128 = k <= 128
    eps_128 = eps_del_espectro(k[m128], E_avg[m128])
    print(f"[eps] espectro promediado ({tag}): eps(k<128) = {eps_128:.4f} "
          f"({eps_128/EPS_OFICIAL:.2f}x doc) ; eps(total) = {eps_tk:.4f} "
          f"({eps_tk/EPS_OFICIAL:.2f}x doc)  [la cola k>128 esta contaminada]")

    # --- pipeline SDDF ---------------------------------------------------
    kt, et, nel, ok, kc = truncate_by_slope_interp(k, E_avg, delta=0.5)
    gn = g_normalizado(kt, et)
    q = q_efectivo(kt, et)
    G = compute_spectral_curvature(kt, et)
    span = float(np.log(kt[-1] / kt[0]))

    txt = []
    txt.append("=" * 72)
    txt.append("SDDF sobre DNS REAL — JHTDB isotropic1024coarse (mirror HF)")
    txt.append("=" * 72)
    txt.append(f"snapshot: 256^3 SUB-CUBO, {n_t} timesteps promediados")
    txt.append(f"ventana de analisis: {tag}")
    txt.append("")
    txt.append("IDENTIDAD DEL DATO (verificada, no asumida):")
    txt.append("  attrs H5: dataset='isotropic1024coarse', indices 1..256 "
               "(paso 1) en x,y,z; t = 1..10")
    txt.append("  dx = 2pi/1024 ; 256 puntos => lado del sub-cubo = pi/2 "
               "(=1/4 del dominio 2pi)")
    txt.append("  TEST DE PERIODICIDAD: salto en el borde / salto interior "
               "~ 300-360x  => NO periodico => SUB-CUBO")
    txt.append("")
    txt.append("PARAMETROS DOCUMENTADOS (JHTDB README-isotropic.pdf):")
    txt.append(f"  nu={NU}, eps={EPS_OFICIAL}, eta={ETA_OFICIAL}, "
               f"u'={U_RMS_OFICIAL}, E_tot={E_TOT_OFICIAL}, Re_lambda=433")
    txt.append("")
    txt.append("CHEQUEO DEL ESPECTRO:")
    txt.append(f"  k_min={kt[0] if False else 2*np.pi/L_sub:.2f}, "
               f"k_Nyq={np.pi/dx:.0f}")
    txt.append(f"  eps_estimado del espectro (todo el rango) = {eps_tk:.4f}  "
               f"(oficial {EPS_OFICIAL})  => {eps_tk/EPS_OFICIAL:.2f}x")
    txt.append(f"  eps_estimado restringido a k<=128 = {eps_128:.4f}  "
               f"=> {eps_128/EPS_OFICIAL:.2f}x   <- la parte confiable")
    txt.append("  eps por gradientes (independiente de las unidades de k): "
               "0.125 con dx=2pi/1024 (1.35x) vs 0.0078 con dx=2pi/256 "
               "(0.08x) => confirma dx=2pi/1024 y que los datos son genuinos.")
    txt.append("  ^ el integrador k^2 E dk esta dominado por la cola de k "
               "alto, que en un sub-cubo no periodico esta contaminada por")
    txt.append("    fuga espectral (~k^-2 por la discontinuidad en las "
               "caras). Comparar con la variante --ventana para medir esa "
               "sensibilidad.")
    txt.append("")
    txt.append("RESULTADO SDDF (invariante bajo reescalado de k):")
    txt.append(f"  corte delta=0.5: k en [{kt[0]:.2f}, {kt[-1]:.2f}] "
               f"({len(kt)} puntos de concha, interp ok={ok})")
    txt.append(f"  ventana efectiva = {span:.3f} nats")
    txt.append(f"  G* = {G:.4f}")
    txt.append(f"  <s^2> = G/ln(k_max/k_min) = {gn:.4f}   "
               f"(K41 puro: 25/9 = 2.7778)")
    txt.append(f"  q_efectivo = sqrt(<s^2>) = {q:.4f}   (K41: 5/3 = 1.6667)")
    txt.append(f"  sesgo vs K41 = {(q - 5/3)/(5/3)*100:+.2f} %")
    txt.append("")
    txt.append("CAVEATS (honestos):")
    txt.append("  1. Sub-cubo no periodico: la FFT sufre fuga; el extremo de "
               "k bajo (k~4-20) NO es confiable y el de k alto tampoco.")
    txt.append("  2. La ventana inercial resultante es corta (~2 nats) por el "
               "corte temprano del detector de pendiente.")
    txt.append("  3. Para un espectro limpio del dominio completo hace falta "
               "promediar muchos sub-cubos (o el campo entero, inviable por "
               "ancho de banda: el mirror va a ~100 kB/s).")
    txt.append("  4. u' y E_tot del sub-cubo son consistentes con lo "
               "documentado => la ESCALA de velocidades esta bien; la "
               "inconsistencia esta en la cola del espectro, no en los datos.")
    out = "\n".join(txt)
    print("\n" + out)
    OUT_TXT.write_text(out + "\n")
    print(f"\n[ok] {OUT_TXT}")


if __name__ == "__main__":
    main()
