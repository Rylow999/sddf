#!/usr/bin/env python3
# exp_migdal_real.py - test log-periodico de Migdal sobre espectros JHTDB REALES.
#
# El null model (16_null_model_periodograma.csv) calibro el test sobre
# espectros SINTETICOS (amp=0.005 no pasa; amp=0.02 pasa, p_bonferroni=0.01).
# Aca aplicamos el MISMO observable a los espectros JHTDB reales cacheados,
# con el null honesto para datos reales: BOOTSTRAP de la incertidumbre por
# concha (E_sem del CSV): E_boot ~ N(E, E_sem), recomputar el score N_BOOT
# veces => p-valor empirico "ruido concha-a-concha solamente".
#
# Salida: datos/22_migdal_real.csv + datos/22_migdal_real_resumen.txt
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from sddf_core import local_slope, inertial_window_smoothed

ESPEC = Path(__file__).parent.parent / "datos" / "espectros"
OUT   = Path(__file__).parent.parent / "datos" / "22_migdal_real.csv"
RES   = Path(__file__).parent.parent / "datos" / "22_migdal_real_resumen.txt"

OM_GRID = np.linspace(0.3, 8.0, 1500)
N_BOOT  = 600
WINDOW  = 5
DELTA   = 0.5
S_REF   = -5.0/3.0

def periodograma(x, y, om):
    y = y - np.polyval(np.polyfit(x, y, 3), x)
    w = np.empty_like(y)
    w[0] = 0.5 * (x[1] - x[0])
    w[-1] = 0.5 * (x[-1] - x[-2])
    w[1:-1] = 0.5 * (x[2:] - x[:-2])
    z = w * y
    P = np.empty(len(om))
    BLOCK = 128
    for i in range(0, len(om), BLOCK):
        omb = om[i:i+BLOCK]
        ph = np.exp(-1j * np.outer(omb, x))
        P[i:i+BLOCK] = np.abs(ph @ z) ** 2
    return P

def score_res(k, E, m):
    lk = np.log(k[m])
    s_loc = local_slope(k[m], E[m]) + 5.0/3.0
    P = periodograma(lk, s_loc, OM_GRID)
    return float(P.max()/np.median(P)), float(OM_GRID[int(np.argmax(P))])

rows = []
lines = []
NOMBRES = ["spectrum_jhtdb_box_hanning.csv",
           "spectrum_jhtdb_box_promedio.csv",
           "spectrum_jhtdb_iso1024_promedio.csv"]
print("== Migdal bootstrap sobre espectros JHTDB reales ==", flush=True)
print(f"N_BOOT={N_BOOT} window={WINDOW} delta={DELTA}", flush=True)
for name in NOMBRES:
    d = np.genfromtxt(ESPEC / name, delimiter=",", names=True)
    k, E = np.asarray(d["k"], float), np.asarray(d["E"], float)
    if "E_sem" in d.dtype.names:
        Esem = np.asarray(d["E_sem"], float)
    else:
        Esem = 0.01 * np.abs(E)   # fallback: 1% del valor local
    Esem = np.where(np.isfinite(Esem) & (Esem > 0), Esem, 0.01*np.abs(E))
    kl, kh, ok = inertial_window_smoothed(k, E, delta=DELTA, s_ref=S_REF, window=WINDOW)
    m = (k >= kl) & (k <= kh)
    span = np.log(kh/kl) if ok else float("nan")
    if not ok or m.sum() < 8:
        lines.append(f"[{name}] ventana insuficiente (ok={ok}, n={int(m.sum())}) - omitido")
        print(lines[-1], flush=True)
        continue
    s_obs, w_obs = score_res(k, E, m)
    rng = np.random.default_rng(2244)
    boots = np.empty(N_BOOT)
    for b in range(N_BOOT):
        Eb = E + rng.standard_normal(len(E)) * Esem
        boots[b] = score_res(k, Eb, m)[0]
    p_emp = float(np.mean(boots >= s_obs))
    q95 = float(np.quantile(boots, 0.95))
    rows.append({"espectro": name, "k_low": round(float(kl),2), "k_high": round(float(kh),2),
                 "span_ln": round(span,3), "n_shells": int(m.sum()),
                 "score_obs": round(s_obs,3), "omega_obs": round(w_obs,3),
                 "null_med": round(float(np.median(boots)),3), "null_q95": round(q95,3),
                 "p_empirico": p_emp})
    lines.append(f"[{name}] window k=[{kl:.1f},{kh:.1f}] span={span:.2f} | score={s_obs:.3f} w={w_obs:.2f} | null med={np.median(boots):.3f} q95={q95:.3f} | p={p_emp:.4f}")
    print(lines[-1], flush=True)

with open(OUT, "w") as f:
    if rows:
        f.write(",".join(rows[0].keys()) + "\n")
        for r in rows:
            f.write(",".join(str(v) for v in r.values()) + "\n")
with open(RES, "w") as f:
    f.write("Test log-periodico de Migdal sobre espectros JHTDB reales (bootstrap E_sem)\n")
    f.write(f"N_BOOT={N_BOOT}, window={WINDOW}, delta={DELTA}, s_ref=-5/3\n\n")
    f.write("\n".join(lines) + "\n\n")
    f.write("Lectura: p_empirico < 0.05 => estructura log-periodica sobre el ruido\n")
    f.write("concha-a-concha. p alto => no detectable a esta resolucion (negativo honesto).\n")
print("MIGDAL_REAL_DONE", flush=True)
