#!/usr/bin/env python3
# exp_dns2d_analisis.py - analisis de los espectros DNS 2D propios (ventana DECLARADA).
#
# El detector automatico (inertial_window_smoothed, delta=0.5) falla en estos
# espectros: 43 conchas => ruido concha-a-concha ~10-15% => la ventana de
# delta se corta (la misma leccion de v3.4, peor resolucion que el box 1024^3).
# Solucion del propio paper: ventana DECLARADA a priori, igual para todos los
# nu (k in [6,20]: fuera del forzante k<=5, antes de la disipacion kd>=23.7).
# Misma ventana en todo nu => el coeficiente a = dG*/dlnRe es comparable.
#
# Salida: datos/24_dns2d_analisis.csv + datos/24_dns2d_analisis_resumen.txt
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from sddf_core import local_slope, compute_spectral_curvature

EDIR = Path(__file__).parent.parent / "datos" / "espectros_2d"
OUT  = Path(__file__).parent.parent / "datos" / "24_dns2d_analisis.csv"
RES  = Path(__file__).parent.parent / "datos" / "24_dns2d_analisis_resumen.txt"

K_LOW, K_HIGH = 6.0, 20.0   # declarada a priori, igual para todo nu
OM2 = np.linspace(0.3, 8.0, 800)

def periodograma2(x, y):
    om = OM2
    y = y - np.polyval(np.polyfit(x, y, 3), x)
    w = np.empty_like(y)
    w[0] = 0.5*(x[1]-x[0]); w[-1] = 0.5*(x[-1]-x[-2]); w[1:-1] = 0.5*(x[2:]-x[:-2])
    z = w*y
    P = np.empty(len(om))
    for i in range(0, len(om), 128):
        omb = om[i:i+128]
        ph = np.exp(-1j*np.outer(omb, x))
        P[i:i+128] = np.abs(ph @ z)**2
    return P

NUS = [7e-4, 5e-4, 3e-4]
rows = []
lines = []
print("== Analisis DNS 2D propia: ventana declarada k=[", K_LOW, ",", K_HIGH, "] ==", flush=True)
for nu in NUS:
    f = EDIR / f"espectro_2d_forzado_nu{nu:g}.csv"
    if not f.exists():
        lines.append(f"nu={nu:g}: espectro no encontrado - omitido")
        print(lines[-1], flush=True)
        continue
    d = np.genfromtxt(f, delimiter=",", names=True)
    k, E = np.asarray(d["k"], float), np.asarray(d["E"], float)
    m = (k >= K_LOW) & (k <= K_HIGH)
    ki, ei = k[m], E[m]
    G = compute_spectral_curvature(ki, ei)
    s = float(np.polyfit(np.log(ki), np.log(ei), 1)[0])
    rows.append({"nu": nu, "ln_inv_nu": round(float(np.log(1.0/nu)),4),
                 "G": round(G,4), "q_medido": round(s,4), "n_shells": int(m.sum())})
    lines.append(f"nu={nu:g}: G={G:.4f}  q={s:.4f}  (n={int(m.sum())} conchas)")
    print(lines[-1], flush=True)

if len(rows) >= 2:
    Gs = np.array([r["G"] for r in rows])
    x  = np.array([r["ln_inv_nu"] for r in rows])
    a_fit = float(np.polyfit(x, Gs, 1)[0])
    qs = np.array([r["q_medido"] for r in rows])
    pred_sc = float(np.median(qs))**2 * 0.5
    lines.append("")
    lines.append(f"a = dG*/dln(1/nu) = {a_fit:.4f}")
    lines.append(f"q medido (finite-Re): {np.median(qs):.4f}  -> q^2*gamma(=1/2) = {pred_sc:.4f}")
    lines.append("ideal asintotico: q=3 -> a = 9/2 = 4.5")
    print(lines[-3], flush=True)
    print(lines[-2], flush=True)

# Migdal sobre el DECAY propio (ventana declarada)
dec_dir = EDIR / "decay"
if dec_dir.exists():
    fs = sorted(dec_dir.glob("espectro_2d_decay_t*.csv"))
    lines.append("")
    lines.append(f"Migdal sobre el decay propio ({len(fs)} tiempos, ventana declarada):")
    print(lines[-1], flush=True)
    sc = []
    for f in fs:
        d = np.genfromtxt(f, delimiter=",", names=True)
        k, E = np.asarray(d["k"], float), np.asarray(d["E"], float)
        m = (k >= K_LOW) & (k <= K_HIGH)
        if m.sum() < 6: continue
        t = float(f.stem.split("_t")[1])
        lk = np.log(k[m])
        s_loc = local_slope(k[m], E[m]) + 3.0
        P = periodograma2(lk, s_loc)
        sc.append((t, float(P.max()/np.median(P))))
        lines.append(f"  t={t:.3f}: score={sc[-1][1]:.3f}")
    if sc:
        arr = np.array([s for _, s in sc])
        lines.append(f"  mediana={np.median(arr):.3f} max={arr.max():.3f} (rango de scores del null sintetico: med 18.7)")
        print(f"  decay: {len(sc)} scores, mediana={np.median(arr):.3f}", flush=True)

with open(OUT, "w") as f:
    if rows:
        f.write(",".join(rows[0].keys()) + "\n")
        for r in rows:
            f.write(",".join(str(v) for v in r.values()) + "\n")
with open(RES, "w") as f:
    f.write("Analisis DNS 2D propia — ventana DECLARADA k=[6,20] (a priori, igual en todo nu)\n")
    f.write("El detector automatico falla a esta resolucion (ruido concha, leccion v3.4)\n\n")
    f.write("\n".join(lines) + "\n")
print("ANALISIS_DONE", flush=True)
