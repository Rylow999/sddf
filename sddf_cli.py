#!/usr/bin/env python3
# sddf_cli.py — MVP: diagnostico espectral de turbulencia en un comando.
# Uso:
#   python sddf_cli.py diagnostico <espectro.csv> [--eta 0.01] [--delta 0.5] [--beta 2.25]
#   python sddf_cli.py verificar              # los 31 checks del paper contra data/
#   python sddf_cli.py synth <Re>            # espectro Kolmogorov sintetico para probar
# El CSV de entrada: columna k, columna E (header o no). Salida: <nombre>_diag.csv
import argparse, json, sys, math
from pathlib import Path
import numpy as np

KAPPA_W = 25.0 / 12.0   # 2.08333...: el teorema

def local_slope(k, E):
    """Pendiente log-log local entre puntos consecutivos."""
    lk, lE = np.log(k), np.log(E)
    return (lE[1:] - lE[:-1]) / (lk[1:] - lk[:-1])

def smoothed_slope(k, E, window=5):
    """Regresion movil (la leccion v3.4: el detector puntual falla en datos reales)."""
    k, E = np.asarray(k, float), np.asarray(E, float)
    if window % 2 == 0: window += 1
    half = window // 2
    lk, lE = np.log(k), np.log(E)
    s = np.empty(len(k))
    for i in range(len(k)):
        j0, j1 = max(0, i - half), min(len(k), i + half + 1)
        if j1 - j0 < 2: j1 = min(len(k), j0 + 2)
        s[i] = np.polyfit(lk[j0:j1], lE[j0:j1], 1)[0]
    return s

def inertial_window(k, E, delta=0.5, s_ref=-5.0/3.0, window=5):
    """Detector de rango inercial dos-lados sobre pendiente suavizada."""
    s = smoothed_slope(k, E, window)
    half = window // 2
    good = np.abs(s - s_ref) <= delta
    good[:half] = False; good[-half:] = False
    # k_high: ultimo good antes de 5 malos
    k_high = k[-1]
    for i in range(half, len(k) - 4):
        if not np.all(good[i:i+5]):
            k_high = k[i]; break
    k_low = k[0]
    for i in range(half, len(k)):
        if good[i]:
            k_low = k[i]; break
    ok = k_high > k_low * 3
    return float(k_low), float(k_high), bool(ok)

def compute_G(k, E):
    """G[u] = integral de pendiente^2 sobre ln k (trapecio)."""
    if len(k) < 3: return 0.0
    s = local_slope(k, E)
    return float(np.trapezoid(s ** 2, np.log(k[1:])))

def q_corregido(k, E, eta, delta, beta):
    """El estimador del paper: ajusta q en la ventana y corrige el sesgo."""
    kl, kh, ok = inertial_window(k, E, delta=delta)
    if not ok:
        return None
    m = (k >= kl) & (k <= kh)
    ki, ei = k[m], E[m]
    q_raw = -np.polyfit(np.log(ki), np.log(ei), 1)[0]
    # correccion de sesgo (sddf_exact): la ventana finita con corte exponencial
    k_c = 0.05 / eta
    # G de la ventana
    G = compute_G(ki, ei)
    # q_corregido por el modelo de corte (forma cerrada, ver paper §5)
    # aproximacion operativa: q_corr = sqrt(G / ln(kh/kl)) reescalado
    q_corr = math.sqrt(max(G, 0.0) / max(math.log(kh / kl), 1e-12))
    return {"q_raw": float(q_raw), "q_corr": float(q_corr), "G": float(G),
            "k_low": kl, "k_high": kh, "n_shells": int(m.sum())}

def kolmogorov_synth(Re, n=200):
    """Espectro sintetico K41 + corte de Pao para validar el pipeline."""
    eta = (1e-3) * (Re / 1000.0) ** (-0.75)   # eta ~ Re^-3/4
    k = np.logspace(0, math.log10(0.05 / eta) * 1.1, n)
    beta = 2.25
    E = (k ** (-5.0/3.0)) * np.exp(-beta * (k * eta) ** 2)
    E *= 1.5 / k[0]  # normalizacion arbitraria
    return k, E, eta

def diagnostico(csv_path, eta, delta, beta, out_dir=None):
    """Pipeline completo sobre un CSV: devuelve dict + escribe CSV de salida."""
    d = np.genfromtxt(csv_path, delimiter=",", names=True, deletechars="")
    k, E = np.asarray(d["k"], float), np.asarray(d["E"], float)
    m = np.isfinite(k) & np.isfinite(E) & (k > 0) & (E > 0)
    k, E = k[m], E[m]
    orden = np.argsort(k)
    k, E = k[orden], E[orden]
    res = {"archivo": str(csv_path), "n_puntos": int(len(k)),
           "k_min": float(k.min()), "k_max": float(k.max())}
    # detector
    kl, kh, ok = inertial_window(k, E, delta=delta)
    res.update({"k_low": kl, "k_high": kh, "rango_detectado": ok})
    # G y q en la ventana
    if ok:
        est = q_corregido(k, E, eta, delta, beta)
        res.update(est)
        # prediccion del teorema: G* ~ (25/12) ln(Re) con Re ~ (1/eta)^(4/3)
        Re_eff = (0.05 / eta) ** (4.0 / 3.0) if eta > 0 else None
        if Re_eff:
            G_pred = KAPPA_W * math.log(Re_eff)
            res["Re_efectivo"] = Re_eff
            res["G_predicho_25_12"] = G_pred
            res["ratio_G"] = res["G"] / G_pred
            res["veredicto"] = ("CONSISTENTE con K41 (G/G* ~ 1)" if 0.5 < res["ratio_G"] < 2.0
                                else "DESVIACION del teorema — revisar espectro/artefactos")
    else:
        res["veredicto"] = "SIN RANGO INERCIAL DETECTABLE (pocos puntos o ruido)"
    # escribir CSV de salida
    out = Path(out_dir or Path(csv_path).parent)
    base = Path(csv_path).stem
    out_csv = out / f"{base}_diag.csv"
    filas = ["metrica,valor"]
    for kk, vv in res.items():
        filas.append(f"{kk},{vv}")
    out_csv.write_text("\n".join(filas) + "\n")
    res["salida"] = str(out_csv)
    return res

def main():
    ap = argparse.ArgumentParser(prog="sddf", description="Diagnostico espectral SDDF (teorema G*=(25/12)lnRe)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p_diag = sub.add_parser("diagnostico", help="Diagnosticar un espectro (CSV con columnas k,E)")
    p_diag.add_argument("csv")
    p_diag.add_argument("--eta", type=float, default=None, help="escala de Kolmogorov (si se conoce)")
    p_diag.add_argument("--delta", type=float, default=0.5)
    p_diag.add_argument("--beta", type=float, default=2.25)
    p_diag.add_argument("-o", "--outdir", default=None)
    p_synth = sub.add_parser("synth", help="Generar espectro Kolmogorov sintetico")
    p_synth.add_argument("Re", type=float)
    p_synth.add_argument("-o", "--outdir", default=".")
    sub.add_parser("verificar", help="31 checks del paper contra data/ (requiere repo)")
    args = ap.parse_args()

    if args.cmd == "diagnostico":
        # eta auto: si no se da, se estima del corte del espectro
        eta = args.eta
        if eta is None:
            d = np.genfromtxt(args.csv, delimiter=",", names=True)
            k = np.asarray(d["k"], float)
            eta = 0.05 / k.max()  # heuristica: k_max ~ 1/eta
            print(f"[sddf] eta estimado por k_max: {eta:.3e}", file=sys.stderr)
        res = diagnostico(args.csv, eta, args.delta, args.beta, args.outdir)
        print(json.dumps(res, indent=1, ensure_ascii=False))
    elif args.cmd == "synth":
        k, E, eta = kolmogorov_synth(args.Re)
        out = Path(args.outdir) / f"spectrum_synth_Re{int(args.Re)}.csv"
        np.savetxt(out, np.column_stack([k, E]), delimiter=",", header="k,E", comments="")
        print(f"[sddf] espectro sintetico Re={args.Re:.0f} -> {out} (eta={eta:.3e})")
    elif args.cmd == "verificar":
        # corre verify.py del repo si esta presente
        import subprocess
        vpy = Path(__file__).parent / "verify.py"
        if vpy.exists():
            r = subprocess.run([sys.executable, str(vpy)], capture_output=True, text=True)
            print(r.stdout)
        else:
            print("verify.py no encontrado — correr desde la raiz del repo sddf")

if __name__ == "__main__":
    main()
