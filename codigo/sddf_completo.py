#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sddf_completo.py — pipeline completo del analisis de curvatura espectral.

Genera todos los CSV y PNG del paquete. Bloques:

   1  reproduccion del bug del lector (piso absoluto 1e-30)
   2  validacion del integrador contra la forma cerrada exacta
   3  sensibilidad del corte por amplitud (por que 165.99 Re^-0.175 es falso)
   4  corte por pendiente -> ley logaritmica (grilla vs interpolado)
   5  descomposicion del prefactor y AJUSTE DESVISCOSADO
   6  barrido exacto Re = 1e2 .. 1e12 y convergencia del prefactor
   7  convergencia de grilla
   8  forma cerrada termino a termino
   9  estimador normalizado <s^2> y exponente efectivo
  10  sesgo de Re finito vs senal de intermitencia
  11  comparacion de modelos espectrales (Pao / quimera / Pope completo)
  12  sensibilidad a beta
  13  prediccion para turbulencia 2D
  14  test log-periodico (firma tipo Migdal)

Uso: python3 sddf_completo.py [dir_espectros] [dir_salida_datos] [dir_figuras]
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import modelos_espectrales as ME
import sddf_exact as EX
from sddf_core import (read_spectrum_csv, compute_spectral_curvature,
                       local_slope, truncate_tail, truncate_by_slope,
                       truncate_by_slope_interp, truncate_by_x,
                       g_normalizado)

IN_DIR = sys.argv[1] if len(sys.argv) > 1 else "../datos/espectros"
OUT_DIR = sys.argv[2] if len(sys.argv) > 2 else "../datos"
FIG_DIR = sys.argv[3] if len(sys.argv) > 3 else "../figuras"
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

EPS, C_K, BETA = 1.0, 1.5, 5.2
NU = {100: 0.01, 200: 0.005, 500: 0.002, 1000: 0.001, 2000: 0.0005}
RE_LIST = sorted(NU)
DELTAS = [0.10, 0.25, 0.50, 1.00]
S_K41 = -5.0 / 3.0

C1, C2, C3, C4 = "#1f5673", "#c1666b", "#48a9a6", "#8a6fa8"
LOG = []


def say(s=""):
    print(s)
    LOG.append(str(s))


def wcsv(name, header, rows):
    p = os.path.join(OUT_DIR, name)
    with open(p, "w") as f:
        f.write(",".join(header) + "\n")
        for r in rows:
            f.write(",".join(str(x) for x in r) + "\n")
    return p


def path_re(re, pref="spectrum"):
    return os.path.join(IN_DIR, f"{pref}_Re_{re}.csv")


def fig(name):
    p = os.path.join(FIG_DIR, name)
    plt.tight_layout()
    plt.savefig(p, dpi=140)
    plt.close()
    return p


# ======================================================================
say("=" * 78)
say("BLOQUE 1 -- REPRODUCCION DEL BUG (piso absoluto 1e-30 en el lector)")
say("=" * 78)
rows1 = []
say(f"{'archivo':14s} {'n':>5s} {'aplast':>7s} {'G_bug':>13s} {'G_sin_piso':>14s} "
    f"{'G*_ampl':>10s} {'G*_pend':>9s}")
for re in RE_LIST:
    p = path_re(re)
    kf, ef, _ = read_spectrum_csv(p, floor=1e-30)
    g_bug = compute_spectral_curvature(kf, ef)
    n_fl = int((ef <= 1e-30).sum())
    k, e, _ = read_spectrum_csv(p, floor=None)
    g_nf = compute_spectral_curvature(k, e)
    kt, et, _ = truncate_tail(k, e, 1e-6)
    g_rel = compute_spectral_curvature(kt, et)
    ks, es, _, ok = truncate_by_slope(k, e, 0.5)
    g_sl = compute_spectral_curvature(ks, es)
    say(f"Re={re:<11d} {len(kf):5d} {n_fl:7d} {g_bug:13.4f} {g_nf:14.2f} "
        f"{g_rel:10.4f} {g_sl:9.4f}")
    rows1.append([f"spectrum_Re_{re}.csv", re, len(kf), n_fl, f"{g_bug:.6f}",
                  f"{g_nf:.6e}", f"{g_rel:.6f}", f"{g_sl:.6f}", int(ok)])
kf, ef, _ = read_spectrum_csv(os.path.join(IN_DIR, "spectrum_fino.csv"), floor=1e-30)
g_fino = compute_spectral_curvature(kf, ef)
say(f"{'fino (K41)':14s} {len(kf):5d} {0:7d} {g_fino:13.4f}   "
    f"exacto (25/9)ln100 = {25/9*np.log(100):.4f}")
rows1.append(["spectrum_fino.csv", "-", len(kf), 0, f"{g_fino:.6f}",
              f"{g_fino:.6e}", f"{g_fino:.6f}", f"{g_fino:.6f}", 1])
wcsv("01_resultados_G.csv",
     ["archivo", "Re", "n_puntos", "n_aplastados", "G_original_bug",
      "G_sin_piso", "G_corte_amplitud_1e-6", "G_corte_pendiente_d0.5",
      "corte_valido"], rows1)
say("Los dos Re mas altos NO tienen puntos aplastados y aun asi daban 1995 y 557:")
say("el piso era solo la mitad del problema; la otra mitad es que la integral")
say("sin ventana definida no converge (la cola aporta x_max^(8/3)).")

# ======================================================================
say()
say("=" * 78)
say("BLOQUE 2 -- VALIDACION DEL INTEGRADOR CONTRA LA FORMA CERRADA")
say("=" * 78)
rows2 = []
for re in RE_LIST:
    k, e, _ = read_spectrum_csv(path_re(re), floor=None)
    gn = compute_spectral_curvature(k, e)
    eta = EX.eta_of(NU[re])
    ga = EX.G_exact(k[0] * eta, k[-1] * eta, beta=BETA)
    err = 100 * abs(gn - ga) / abs(ga)
    say(f"Re={re:5d} eta={eta:.6e}  num={gn:14.4f}  exacto={ga:14.4f}  err={err:.6f} %")
    rows2.append([re, NU[re], f"{eta:.9e}", f"{gn:.9e}", f"{ga:.9e}", f"{err:.8f}"])
wcsv("02_validacion_analitica.csv",
     ["Re", "nu", "eta", "G_numerico", "G_exacto", "error_pct"], rows2)

# ======================================================================
say()
say("=" * 78)
say("BLOQUE 3 -- SENSIBILIDAD DEL CORTE POR AMPLITUD (el falso Re^-0.175)")
say("=" * 78)
umbrales = [1e-3, 1e-4, 1e-5, 1e-6, 1e-8, 1e-10, 1e-12]
rows3, sens = [], {}
say(f"{'umbral':>8s} " + " ".join(f"{'Re'+str(r):>9s}" for r in RE_LIST)
    + f" {'A':>9s} {'alpha':>8s} {'R2':>9s}")
for th in umbrales:
    gs = []
    for re in RE_LIST:
        k, e, _ = read_spectrum_csv(path_re(re), floor=None)
        kt, et, _ = truncate_tail(k, e, th)
        gs.append(compute_spectral_curvature(kt, et))
    A, al, r2 = EX.fit_power(RE_LIST, gs)
    sens[th] = (gs, A, al, r2)
    say(f"{th:8.0e} " + " ".join(f"{g:9.3f}" for g in gs)
        + f" {A:9.2f} {al:8.4f} {r2:9.5f}")
    rows3.append([f"{th:.0e}"] + [f"{g:.6f}" for g in gs]
                 + [f"{A:.4f}", f"{al:.6f}", f"{r2:.6f}"])
wcsv("03_sensibilidad_umbral.csv",
     ["umbral"] + [f"G_Re{r}" for r in RE_LIST] + ["A", "alpha", "R2"], rows3)
als = [sens[t][2] for t in umbrales]
say(f"alpha varia entre {min(als):.4f} y {max(als):.4f} (factor {max(als)/min(als):.2f})")
say("con R2 > 0.97 SIEMPRE. El R2 alto mide suavidad, no fisica.")

# ======================================================================
say()
say("=" * 78)
say("BLOQUE 4 -- CORTE POR PENDIENTE: LEY LOGARITMICA")
say("=" * 78)
rows4, G_by_delta, G_by_delta_int = [], {}, {}
say(f"{'delta':>6s} {'x_c teor':>9s} " + " ".join(f"{'Re'+str(r):>8s}" for r in RE_LIST)
    + f" {'a':>8s} {'b':>9s} {'R2':>8s}")
for d in DELTAS:
    xc_t = EX.x_c(d, beta=BETA)
    gs, gsi, kcs, xcs = [], [], [], []
    for re in RE_LIST:
        k, e, _ = read_spectrum_csv(path_re(re), floor=None)
        ks, es, _, ok = truncate_by_slope(k, e, d)
        gs.append(compute_spectral_curvature(ks, es))
        ki, ei, _, ok2, kc = truncate_by_slope_interp(k, e, d)
        gsi.append(compute_spectral_curvature(ki, ei))
        kcs.append(kc)
        xcs.append(kc * EX.eta_of(NU[re]))
    a, b, r2 = EX.fit_log(RE_LIST, gs)
    ai, bi, r2i = EX.fit_log(RE_LIST, gsi)
    G_by_delta[d], G_by_delta_int[d] = gs, gsi
    say(f"{d:6.2f} {xc_t:9.6f} " + " ".join(f"{g:8.4f}" for g in gs)
        + f" {a:8.5f} {b:9.5f} {r2:8.5f}")
    say(f"{'':6s} {'x_c med':>9s} " + " ".join(f"{x:8.5f}" for x in xcs)
        + f"   [interp: a={ai:.5f} b={bi:.5f} R2={r2i:.5f}]")
    rows4.append([d, f"{xc_t:.8f}"] + [f"{g:.6f}" for g in gs]
                 + [f"{g:.6f}" for g in gsi] + [f"{k:.4f}" for k in kcs]
                 + [f"{x:.6f}" for x in xcs]
                 + [f"{a:.6f}", f"{b:.6f}", f"{r2:.6f}",
                    f"{ai:.6f}", f"{bi:.6f}", f"{r2i:.6f}"])
wcsv("04_ley_logaritmica.csv",
     ["delta", "x_c_teorico"] + [f"G_grilla_Re{r}" for r in RE_LIST]
     + [f"G_interp_Re{r}" for r in RE_LIST] + [f"k_corte_Re{r}" for r in RE_LIST]
     + [f"x_c_medido_Re{r}" for r in RE_LIST]
     + ["a_grilla", "b_grilla", "R2_grilla", "a_interp", "b_interp", "R2_interp"], rows4)

# ======================================================================
say()
say("=" * 78)
say("BLOQUE 5 -- PREFACTOR: Re FINITO vs ASINTOTA, Y AJUSTE DESVISCOSADO")
say("=" * 78)
rows5 = []
say(f"{'delta':>6s} {'a_num':>9s} {'a_interp':>9s} {'a_exacto':>9s} {'a_DESVISC':>10s} "
    f"{'25/12':>8s} {'sesgo%':>8s} {'mu_eff':>8s}")
for d in DELTAS:
    gs = np.array(G_by_delta[d])
    gsi = np.array(G_by_delta_int[d])
    gex = EX.G_star(np.array(RE_LIST, float), d, beta=BETA)
    a, _, _ = EX.fit_log(RE_LIST, gs)
    ai, _, _ = EX.fit_log(RE_LIST, gsi)
    ae, _, _ = EX.fit_log(RE_LIST, gex)
    adv, _, _ = EX.fit_log(RE_LIST, EX.deviscosify(RE_LIST, gsi, beta=BETA))
    sesgo = 100 * (ae - EX.A_ASYMP_K41) / EX.A_ASYMP_K41
    mu = EX.mu_from_a(ae)
    say(f"{d:6.2f} {a:9.5f} {ai:9.5f} {ae:9.5f} {adv:10.5f} "
        f"{EX.A_ASYMP_K41:8.5f} {sesgo:8.3f} {mu:8.5f}")
    rows5.append([d, f"{EX.x_c(d, beta=BETA):.8f}", f"{a:.6f}", f"{ai:.6f}",
                  f"{ae:.6f}", f"{adv:.6f}", f"{EX.A_ASYMP_K41:.6f}",
                  f"{sesgo:.4f}", f"{mu:.6f}"])
wcsv("05_prefactor.csv",
     ["delta", "x_c", "a_numerico_grilla", "a_numerico_interp", "a_exacto_misma_ventana",
      "a_desviscosado", "a_asintotico_25_12", "sesgo_pct", "mu_efectivo_implicado"], rows5)
say("El ajuste desviscosado -- restarle a G* los terminos -2q*beta/Re y")
say("-(beta^2 p/2)/Re^2 ANTES de ajustar -- recupera 25/12 con 5 puntos.")

# ======================================================================
say()
say("=" * 78)
say("BLOQUE 6 -- BARRIDO EXACTO Re = 1e2 .. 1e12")
say("=" * 78)
RE_EXT = np.logspace(2, 12, 41)
rows6 = []
for d in DELTAS:
    g = EX.G_star(RE_EXT, d, beta=BETA)
    rows6.append([d] + [f"{v:.6f}" for v in g])
wcsv("06_barrido_extendido.csv", ["delta"] + [f"G_Re{r:.3e}" for r in RE_EXT], rows6)

REMAX = [2e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e10, 1e12]
rows6b = []
say(f"{'Re_max':>10s} {'a_cruda':>10s} {'a_desviscosada':>16s} {'desvio vs 25/12 %':>18s}")
for rm in REMAX:
    rr = np.logspace(2, np.log10(rm), 40)
    g = EX.G_star(rr, 0.5, beta=BETA)
    a, _, _ = EX.fit_log(rr, g)
    adv, _, _ = EX.fit_log(rr, EX.deviscosify(rr, g, beta=BETA))
    say(f"{rm:10.0e} {a:10.6f} {adv:16.6f} {100*(a-EX.A_ASYMP_K41)/EX.A_ASYMP_K41:18.4f}")
    rows6b.append([f"{rm:.0e}", f"{a:.6f}", f"{adv:.6f}",
                   f"{100*(a-EX.A_ASYMP_K41)/EX.A_ASYMP_K41:.6f}"])
wcsv("06b_pendiente_por_ventana.csv",
     ["Re_max", "a_cruda", "a_desviscosada", "desvio_vs_25_12_pct"], rows6b)

# ======================================================================
say()
say("=" * 78)
say("BLOQUE 7 -- CONVERGENCIA DE GRILLA (corte de grilla vs interpolado)")
say("=" * 78)
rows7 = []
g_ex_ref = EX.G_star(100.0, 0.5, beta=BETA)
say(f"{'N':>8s} {'G_grilla':>12s} {'err%':>9s} {'G_interp':>12s} {'err%':>9s}  "
    f"(exacto={g_ex_ref:.6f})")
for N in [100, 200, 500, 1000, 2000, 5000, 20000, 100000]:
    k = np.logspace(0, 3, N)
    e = ME.E_quimera(k, 0.01, beta=BETA)
    ks, es, _, _ = truncate_by_slope(k, e, 0.5)
    g1 = compute_spectral_curvature(ks, es)
    ki, ei, _, _, _ = truncate_by_slope_interp(k, e, 0.5)
    g2 = compute_spectral_curvature(ki, ei)
    e1 = 100 * abs(g1 - g_ex_ref) / g_ex_ref
    e2 = 100 * abs(g2 - g_ex_ref) / g_ex_ref
    say(f"{N:8d} {g1:12.6f} {e1:9.4f} {g2:12.6f} {e2:9.4f}")
    rows7.append([N, f"{g1:.8f}", f"{e1:.6f}", f"{g2:.8f}", f"{e2:.6f}",
                  f"{g_ex_ref:.8f}"])
wcsv("07_convergencia_grilla.csv",
     ["N_puntos", "G_corte_grilla", "error_grilla_pct", "G_corte_interp",
      "error_interp_pct", "G_exacto"], rows7)

# ======================================================================
say()
say("=" * 78)
say("BLOQUE 8 -- FORMA CERRADA TERMINO A TERMINO")
say("=" * 78)
rows8 = []
say(f"{'delta':>6s} {'Re':>6s} {'(3/4)q^2 lnRe':>14s} {'b(delta)':>10s} "
    f"{'t1 ~ -1/Re':>12s} {'t2 ~ -1/Re^2':>13s} {'G* total':>10s}")
for d in DELTAS:
    for re in RE_LIST:
        t1, t2 = EX.finite_Re_terms(re, beta=BETA)
        lead = EX.A_ASYMP_K41 * np.log(re)
        b = EX.b_delta(d, beta=BETA)
        tot = lead + b + t1 + t2
        say(f"{d:6.2f} {re:6d} {lead:14.6f} {b:10.6f} {t1:12.6f} {t2:13.6e} {tot:10.6f}")
        rows8.append([d, re, f"{lead:.8f}", f"{b:.8f}", f"{t1:.8e}", f"{t2:.8e}",
                      f"{tot:.8f}"])
wcsv("08_forma_cerrada_terminos.csv",
     ["delta", "Re", "termino_log_25_12", "b_delta", "t1_Re_menos1",
      "t2_Re_menos2", "G_star_exacto"], rows8)

# ======================================================================
say()
say("=" * 78)
say("BLOQUE 9 -- ESTIMADOR NORMALIZADO <s^2> = G*/ln(k_c/k_min)")
say("=" * 78)
rows9 = []
say(f"{'delta':>6s} {'Re':>6s} {'<s^2>':>9s} {'q_eff':>8s} {'mu_eff':>9s} "
    f"{'q_corregido':>12s} {'mu_corr':>9s}   (K41: 2.77778 / 1.66667)")


def q_corregido(G, k_min, k_c, eta, delta, beta, p=4/3, q0=5/3, iters=3):
    """
    Despeja q de la forma cerrada en vez de dividir a lo bruto:
      G = q^2 ln(x_c/x1) + 2 q beta (x_c^p - x1^p) + (beta^2 p/2)(x_c^2p - x1^2p)
    Los dos ultimos terminos son conocidos una vez fijados delta y beta.
    """
    x1, xc_ = k_min * eta, k_c * eta
    L = np.log(xc_ / x1)
    t3 = (beta ** 2 * p / 2.0) * (xc_ ** (2 * p) - x1 ** (2 * p))
    q = q0
    for _ in range(iters):
        t2 = 2.0 * q * beta * (xc_ ** p - x1 ** p)
        q = np.sqrt(max(G - t2 - t3, 0.0) / L)
    return q


for d in [0.10, 0.25]:
    for re in RE_LIST:
        k, e, _ = read_spectrum_csv(path_re(re), floor=None)
        ki, ei, _, _, kc = truncate_by_slope_interp(k, e, d)
        gn = g_normalizado(ki, ei)
        q = np.sqrt(gn)
        eta = EX.eta_of(NU[re])
        qc = q_corregido(compute_spectral_curvature(ki, ei), ki[0], kc, eta, d, BETA)
        rows9.append([d, re, f"{gn:.6f}", f"{q:.6f}", f"{q-5/3:.6f}",
                      f"{qc:.6f}", f"{qc-5/3:.6f}"])
        say(f"{d:6.2f} {re:6d} {gn:9.5f} {q:8.5f} {q-5/3:9.5f} {qc:12.6f} {qc-5/3:9.6f}")
wcsv("09_estimador_normalizado.csv",
     ["delta", "Re", "s2_medio", "q_efectivo_crudo", "mu_efectivo_crudo",
      "q_corregido", "mu_corregido"], rows9)
say("<s^2> crudo NO es invariante de Re: converge a 25/9 solo como 1/lnRe,")
say("que es MAS lento que el ajuste de pendiente. A Re=100 con delta=0.10 el")
say("sesgo es mu_eff = 0.084, tres veces la intermitencia fisica.")
say("Despejando q de la forma cerrada (columna q_corregido) el sesgo desaparece.")

# ======================================================================
say()
say("=" * 78)
say("BLOQUE 10 -- SESGO DE Re FINITO vs SENAL DE INTERMITENCIA")
say("=" * 78)
rows10 = []
say(f"{'mu':>7s} {'a=(3/4)(5/3+mu)^2':>19s}   interpretacion")
for mu, txt in [(0.000, "K41 puro"), (0.0212, "sesgo de Re finito del barrido"),
                (0.025, "intermitencia baja"), (0.033, "zeta_2 = 0.70 (DNS)"),
                (0.040, "intermitencia alta")]:
    a = EX.a_from_mu(mu)
    say(f"{mu:7.4f} {a:19.6f}   {txt}")
    rows10.append([f"{mu:.4f}", f"{a:.6f}", txt])
wcsv("10_sesgo_vs_intermitencia.csv",
     ["mu", "prefactor_a", "interpretacion"], rows10)
say("El sesgo de Re finito (a = 2.1366, mu_eff = 0.0212) es del MISMO ORDEN")
say("que la intermitencia fisica (mu ~ 0.033). Sin desviscosar, G[u] no puede")
say("distinguir 'rango inercial corto' de 'escalamiento anomalo'.")

# --- verificacion directa: espectro con mu=0.033 inyectado ------------
say()
say("Verificacion: barrido con mu=0.033 inyectado (modelo pao_beta2.25):")
gsm = []
for re in RE_LIST:
    k, e, _ = read_spectrum_csv(os.path.join(IN_DIR, f"intermitente_Re_{re}.csv"), floor=None)
    ki, ei, _, _, _ = truncate_by_slope_interp(k, e, 0.25, s_ref=-(5/3+0.033))
    gsm.append(compute_spectral_curvature(ki, ei))
a_m, _, _ = EX.fit_log(RE_LIST, gsm)
a_dv, _, _ = EX.fit_log(RE_LIST, EX.deviscosify(RE_LIST, gsm, beta=ME.BETA_PAO,
                                                q=5/3+0.033))
say(f"  a_cruda = {a_m:.5f}   a_desviscosada = {a_dv:.5f}   "
    f"teorico (3/4)(5/3+0.033)^2 = {EX.a_from_mu(0.033):.5f}")
say(f"  mu recuperado = {EX.mu_from_a(a_dv):.5f}  (inyectado 0.033)")

# ======================================================================
say()
say("=" * 78)
say("BLOQUE 11 -- COMPARACION DE MODELOS ESPECTRALES")
say("=" * 78)
rows11 = []
say(f"{'modelo':32s} {'beta':>6s} {'a_cruda':>9s} {'a_desvisc':>10s} {'b ajustada':>11s}")
# (beta_eff, p_eff) del comportamiento de f_eta a x pequeno, que es lo unico
# que entra en la correccion de Re finito. Para Pope, f_eta ~ exp(-beta x^4/(4 c^3)),
# o sea p_eff = 4 y la correccion decae como Re^-3: despreciable.
for pref, nombre, beta_m, p_m in [
        ("spectrum", "quimera Pao4/3 + beta Pope 5.2", ME.BETA_POPE, 4 / 3),
        ("pao", "Pao consistente beta=2.25", ME.BETA_PAO, 4 / 3),
        ("pope", "Pope completo f_L*f_eta", ME.BETA_POPE / (4 * ME.C_ETA_POPE ** 3), 4.0)]:
    gs = []
    for re in RE_LIST:
        p = os.path.join(IN_DIR, f"{pref}_Re_{re}.csv")
        k, e, _ = read_spectrum_csv(p, floor=None)
        if pref == "pope":
            pass   # con L=30 el aporte de f_L en k=1 ya es |ds|=0.027
        ki, ei, _, _, _ = truncate_by_slope_interp(k, e, 0.25)
        gs.append(compute_spectral_curvature(ki, ei))
    a, b, r2 = EX.fit_log(RE_LIST, gs)
    adv, _, _ = EX.fit_log(RE_LIST, EX.deviscosify(RE_LIST, gs, beta=beta_m, p=p_m))
    say(f"{nombre:32.32s} {beta_m:6.2f} {a:9.5f} {adv:10.5f} {b:11.5f}")
    rows11.append([nombre, f"{beta_m:.4f}", f"{p_m:.4f}"] + [f"{g:.6f}" for g in gs]
                  + [f"{a:.6f}", f"{adv:.6f}", f"{b:.6f}", f"{r2:.6f}"])
wcsv("11_modelos_espectrales.csv",
     ["modelo", "beta_eff", "p_eff"] + [f"G_Re{r}" for r in RE_LIST]
     + ["a_cruda", "a_desviscosada", "b", "R2"], rows11)
say("El prefactor sobrevive al cambio de modelo de corte; la ordenada no.")

# ======================================================================
say()
say("=" * 78)
say("BLOQUE 12 -- SENSIBILIDAD A beta (el parametro mal elegido)")
say("=" * 78)
rows12 = []
say(f"{'beta':>6s} {'x_c(0.25)':>10s} {'b(0.25)':>10s} {'a asintotico':>13s}")
for bt in [2.25, 2.40, 3.0, 4.0, 5.2, 7.2]:
    xc_ = EX.x_c(0.25, beta=bt)
    b_ = EX.b_delta(0.25, beta=bt)
    say(f"{bt:6.2f} {xc_:10.6f} {b_:10.6f} {EX.A_ASYMP_K41:13.6f}")
    rows12.append([bt, f"{xc_:.8f}", f"{b_:.8f}", f"{EX.A_ASYMP_K41:.8f}"])
wcsv("12_sensibilidad_beta.csv",
     ["beta", "x_c_delta0.25", "b_delta0.25", "a_asintotico"], rows12)
say(f"Pasar de beta=5.2 (Pope) a beta=2.25 (Pao consistente) mueve x_c un factor "
    f"{EX.x_c(0.25, beta=2.25)/EX.x_c(0.25, beta=5.2):.3f} y b en "
    f"{EX.b_delta(0.25, beta=2.25)-EX.b_delta(0.25, beta=5.2):+.3f}; a no se mueve.")

# ======================================================================
say()
say("=" * 78)
say("BLOQUE 13 -- PREDICCION PARA EL PAPER 2D")
say("=" * 78)
rows13 = [
    ["3D K41 (inercial)", "5/3", "3/4", f"{EX.a_general(5/3, 3/4):.6f}", "25/12"],
    ["3D con mu=0.033", "1.700", "3/4", f"{EX.a_general(1.7, 3/4):.6f}", "-"],
    ["2D cascada de enstrofia", "3", "1/2", f"{EX.a_general(3.0, 0.5):.6f}", "9/2"],
    ["2D Kraichnan-Batchelor con log", "3", "1/2", f"{EX.a_general(3.0, 0.5):.6f}", "9/2"],
    ["2D cascada inversa", "5/3", "0", f"{EX.a_general(5/3, 0.0):.6f}", "0"],
]
say(f"{'regimen':32s} {'q':>7s} {'d lnk_c/d lnRe':>15s} {'a = q^2 * exp':>14s}")
for r in rows13:
    say(f"{r[0]:32s} {r[1]:>7s} {r[2]:>15s} {r[3]:>14s}")
wcsv("13_prediccion_2d.csv",
     ["regimen", "exponente_espectral_q", "exponente_de_k_corte", "a_predicho", "fraccion"],
     rows13)
say("El paper 2D reporta G* ~ 3634 y G* ~ Re^0.70. Con q=3, G*=3634 exigiria")
say(f"ln(k_c/k_min) = 3634/9 = {3634/9:.1f} nats, o sea k_c/k_min = e^{3634/9:.0f}.")
say("Imposible: es el mismo artefacto de cola sin ventana que aca daba 182017.")

# ======================================================================
say()
say("=" * 78)
say("BLOQUE 14 -- TEST LOG-PERIODICO (firma tipo Migdal)")
say("=" * 78)
say("Ventana inercial larga a proposito: Re=1e10, k in [1,1e7], 200k puntos.")
say("Con la ventana corta del barrido original (2.7 nats de ln k) NO alcanza")
say("ni para un ciclo y la deteccion es imposible: ese es el primer resultado.")
NU_LP = 1e-10 ** 0.0 * 1e-10   # nu = 1e-10 -> Re = 1e10
k_lp = np.logspace(0, 7, 200000)
eta_lp = EX.eta_of(NU_LP)
OMEGA_INY = 2.0
rows14 = []


def periodograma(x, y, om):
    y = y - np.polyval(np.polyfit(x, y, 3), x)
    return np.array([abs(np.trapezoid(y * np.exp(-1j * w * x), x)) ** 2 for w in om])


om_grid = np.linspace(0.3, 8.0, 4000)
say(f"{'amp':>6s} {'ventana(nats)':>14s} {'Delta<s^2>':>12s} {'Delta q_eff':>12s} "
    f"{'omega detectada':>16s} {'SNR':>8s}")
for amp in [0.0, 0.002, 0.005, 0.02, 0.05]:
    e_lp = ME.E_log_periodico(k_lp, NU_LP, beta=ME.BETA_PAO, amp=amp, omega=OMEGA_INY)
    e_0 = ME.E_pao(k_lp, NU_LP, beta=ME.BETA_PAO)
    # ventana FIJA en k*eta para que la comparacion sea limpia
    x_lo, x_hi = 1e-6, EX.x_c(0.25, beta=ME.BETA_PAO)
    m = (k_lp * eta_lp > x_lo) & (k_lp * eta_lp < x_hi)
    gn1 = g_normalizado(k_lp[m], e_lp[m])
    gn0 = g_normalizado(k_lp[m], e_0[m])
    lk = np.log(k_lp[m])
    s_loc = local_slope(k_lp[m], e_lp[m]) + 5 / 3
    P = periodograma(lk, s_loc, om_grid)
    w_pk = om_grid[int(np.argmax(P))]
    snr = P.max() / np.median(P)
    span = lk[-1] - lk[0]
    say(f"{amp:6.3f} {span:14.2f} {gn1-gn0:12.3e} {np.sqrt(gn1)-np.sqrt(gn0):12.3e} "
        f"{w_pk:16.4f} {snr:8.1f}")
    rows14.append([amp, OMEGA_INY, f"{span:.4f}", f"{gn1-gn0:.6e}",
                   f"{np.sqrt(gn1)-np.sqrt(gn0):.6e}", f"{w_pk:.6f}", f"{snr:.3f}"])

# el mismo test en la ventana corta del barrido original
e_s = ME.E_log_periodico(K_SHORT := np.logspace(0, 3, 20000), 5e-4,
                         beta=ME.BETA_PAO, amp=0.02, omega=OMEGA_INY)
ki, ei, _, _, _ = truncate_by_slope_interp(K_SHORT, e_s, 0.25)
lk_s = np.log(ki)
P_s = periodograma(lk_s, local_slope(ki, ei) + 5 / 3, om_grid)
say(f"{'0.020':>6s} {lk_s[-1]-lk_s[0]:14.2f} {'-':>12s} {'-':>12s} "
    f"{om_grid[int(np.argmax(P_s))]:16.4f} {P_s.max()/np.median(P_s):8.1f}   <-- ventana corta: falla")
rows14.append([0.02, OMEGA_INY, f"{lk_s[-1]-lk_s[0]:.4f}", "", "",
               f"{om_grid[int(np.argmax(P_s))]:.6f}",
               f"{P_s.max()/np.median(P_s):.3f}"])

wcsv("14_test_log_periodico.csv",
     ["amplitud", "omega_inyectada", "ventana_ln_k_nats", "delta_s2_medio",
      "delta_q_efectivo", "omega_detectada", "SNR_periodograma"], rows14)
say("G[u] casi no se mueve (Delta q_eff ~ 1e-4 con amp=0.02), pero el")
say("periodograma de s(ln k) clava omega si la ventana tiene varios ciclos.")
say("Conclusion: para contrastar la escalera de Stokes de Migdal el observable")
say("NO es G[u]; es el espectro de Fourier de la pendiente local en ln k.")

# ======================================================================
# FIGURAS
# ======================================================================
say()
say("=" * 78)
say("FIGURAS")
say("=" * 78)

# fig1: diagnostico del bug
k100, e100, _ = read_spectrum_csv(path_re(100), floor=None)
kb, eb, _ = read_spectrum_csv(path_re(100), floor=1e-30)
f, ax = plt.subplots(1, 3, figsize=(15, 4.2))
ax[0].loglog(kb, eb, color=C2, lw=1.6, label="con piso 1e-30")
ax[0].loglog(k100, e100, color=C1, lw=1.2, ls="--", label="sin piso")
ax[0].set_xlabel("k"); ax[0].set_ylabel("E(k)")
ax[0].set_title("Re=100: la meseta artificial"); ax[0].legend(fontsize=8)
ax[1].semilogx(kb, local_slope(kb, eb), color=C2, lw=1.4, label="con piso")
ax[1].semilogx(k100, local_slope(k100, e100), color=C1, lw=1.2, ls="--", label="sin piso")
ax[1].axhline(-5 / 3, color="k", lw=0.8, ls=":")
ax[1].set_ylim(-60, 5); ax[1].set_xlabel("k"); ax[1].set_ylabel("d lnE / d lnk")
ax[1].set_title("salto espurio en la derivada"); ax[1].legend(fontsize=8)
ax[2].loglog(k100, local_slope(k100, e100) ** 2, color=C1, lw=1.2)
ax[2].set_xlabel("k"); ax[2].set_ylabel("s(k)$^2$")
ax[2].set_title("integrando: diverge como $x^{8/3}$")
say(fig("fig1_diagnostico_bug.png"))

# fig2: espectros y pendientes
f, ax = plt.subplots(1, 2, figsize=(11, 4.2))
cmap = plt.get_cmap("viridis")
for i, re in enumerate(RE_LIST):
    k, e, _ = read_spectrum_csv(path_re(re), floor=None)
    c = cmap(i / 4)
    ax[0].loglog(k, e, color=c, lw=1.3, label=f"Re={re}")
    ax[1].semilogx(k, local_slope(k, e), color=c, lw=1.3)
ax[0].loglog(k, 1.5 * k ** (-5 / 3), "k:", lw=0.9, label="$k^{-5/3}$")
ax[0].set_xlabel("k"); ax[0].set_ylabel("E(k)"); ax[0].legend(fontsize=8)
ax[0].set_title("barrido en Re")
ax[1].axhline(-5 / 3, color="k", ls=":", lw=0.9)
ax[1].set_ylim(-12, 0); ax[1].set_xlabel("k"); ax[1].set_ylabel("d lnE / d lnk")
ax[1].set_title("pendiente local: el corte es $k\\eta$ = cte")
say(fig("fig2_espectros.png"))

# fig3: los tres estimadores
f, ax = plt.subplots(1, 3, figsize=(15, 4.2))
ax[0].loglog(RE_LIST, [float(r[4]) for r in rows1[:5]], "o-", color=C2)
ax[0].set_title("G con el bug (sin sentido)"); ax[0].set_xlabel("Re"); ax[0].set_ylabel("G")
ax[1].loglog(RE_LIST, sens[1e-6][0], "o-", color=C4)
A, al, r2 = sens[1e-6][1], sens[1e-6][2], sens[1e-6][3]
rr = np.array(RE_LIST, float)
ax[1].loglog(rr, A * rr ** (-al), "k--", lw=0.9)
ax[1].set_title(f"corte por amplitud: $Re^{{-{al:.3f}}}$ (artefacto)")
ax[1].set_xlabel("Re")
for d in DELTAS:
    ax[2].semilogx(RE_LIST, G_by_delta_int[d], "o-", label=f"$\\delta$={d}")
ax[2].semilogx(rr, EX.A_ASYMP_K41 * np.log(rr) + EX.b_delta(0.5, beta=BETA),
               "k--", lw=0.9, label="(25/12) ln Re + b")
ax[2].set_title("corte por pendiente: ley logaritmica")
ax[2].set_xlabel("Re"); ax[2].legend(fontsize=8)
say(fig("fig3_G_vs_Re.png"))

# fig4: sensibilidad del umbral
f, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ax[0].semilogx(umbrales, [sens[t][2] for t in umbrales], "o-", color=C2)
ax[0].set_xlabel("umbral de amplitud"); ax[0].set_ylabel(r"$\alpha$ ajustado")
ax[0].set_title(r"$\alpha$ es una funcion del umbral, no del flujo")
ax[1].semilogx(umbrales, [sens[t][3] for t in umbrales], "s-", color=C1)
ax[1].set_ylim(0.95, 1.001); ax[1].set_xlabel("umbral de amplitud")
ax[1].set_ylabel("$R^2$"); ax[1].set_title("y el $R^2$ siempre es alto")
say(fig("fig4_sensibilidad.png"))

# fig5: convergencia
f, ax = plt.subplots(1, 2, figsize=(11, 4.2))
rr = np.array([float(r[0]) for r in rows6b])
ax[0].semilogx(rr, [float(r[1]) for r in rows6b], "o-", color=C1, label="ajuste crudo")
ax[0].semilogx(rr, [float(r[2]) for r in rows6b], "s-", color=C3, label="desviscosado")
ax[0].axhline(EX.A_ASYMP_K41, color="k", ls="--", lw=0.9, label="25/12")
ax[0].set_xlabel("$Re_{max}$ de la ventana"); ax[0].set_ylabel("prefactor a")
ax[0].legend(fontsize=8); ax[0].set_title("convergencia del prefactor")
NN = [int(r[0]) for r in rows7]
ax[1].loglog(NN, [float(r[2]) for r in rows7], "o-", color=C1, label="corte en grilla")
ax[1].loglog(NN, [max(float(r[4]), 1e-8) for r in rows7], "s-", color=C3, label="corte interpolado")
ax[1].set_xlabel("N puntos"); ax[1].set_ylabel("error [%]")
ax[1].legend(fontsize=8); ax[1].set_title("error de discretizacion")
say(fig("fig5_convergencia.png"))

# fig6: ley log extendida
f, ax = plt.subplots(figsize=(7, 4.6))
for i, d in enumerate(DELTAS):
    ax.semilogx(RE_EXT, EX.G_star(RE_EXT, d, beta=BETA), color=cmap(i / 3),
                lw=1.5, label=f"$\\delta$={d}")
    ax.semilogx(RE_LIST, G_by_delta_int[d], "o", color=cmap(i / 3), ms=5)
ax.set_xlabel("Re"); ax.set_ylabel("$G^*[u]$")
ax.set_title("$G^* = (25/12)\\,\\ln Re + b(\\delta) + O(Re^{-1})$  (puntos: numerico)")
ax.legend(fontsize=8)
say(fig("fig6_ley_log_extendida.png"))

# fig7: sesgo vs intermitencia
f, ax = plt.subplots(figsize=(7, 4.6))
mus = np.linspace(0, 0.06, 200)
ax.plot(mus, EX.a_from_mu(mus), color=C1, lw=1.8, label=r"$a=(3/4)(5/3+\mu)^2$")
ax.axhline(2.136648, color=C2, ls="--", lw=1.2,
           label="a medido sin desviscosar = 2.1366")
ax.axhline(EX.A_ASYMP_K41, color="k", ls=":", lw=1.0, label="25/12 (K41)")
ax.axvline(0.033, color=C3, ls="-.", lw=1.0, label=r"$\mu\approx0.033$ ($\zeta_2=0.70$)")
ax.plot([0.0212], [2.136648], "o", color=C2, ms=8)
ax.set_xlabel(r"$\mu$ (correccion de intermitencia)")
ax.set_ylabel("prefactor a"); ax.legend(fontsize=8)
ax.set_title("el sesgo de Re finito imita intermitencia")
say(fig("fig7_sesgo_intermitencia.png"))

# fig8: modelos espectrales
f, ax = plt.subplots(1, 2, figsize=(11, 4.2))
kk = np.logspace(0, 3, 3000)
for nm, fn, c in [("Pao beta=2.25", lambda k: ME.E_pao(k, 0.0005, beta=ME.BETA_PAO), C1),
                  ("quimera beta=5.2", lambda k: ME.E_quimera(k, 0.0005, beta=ME.BETA_POPE), C2),
                  ("Pope completo", lambda k: ME.E_pope(k, 0.0005, L=30.0), C3)]:
    ee = fn(kk)
    ax[0].loglog(kk, ee, color=c, lw=1.4, label=nm)
    ax[1].semilogx(kk, local_slope(kk, ee), color=c, lw=1.4, label=nm)
ax[0].set_xlabel("k"); ax[0].set_ylabel("E(k)"); ax[0].legend(fontsize=8)
ax[0].set_title("Re=2000: tres parametrizaciones")
ax[1].axhline(-5 / 3, color="k", ls=":", lw=0.9)
ax[1].set_ylim(-8, 1); ax[1].set_xlabel("k"); ax[1].set_ylabel("d lnE/d lnk")
ax[1].set_title("el corte llega 2.3x antes con beta=5.2")
say(fig("fig8_modelos.png"))

# fig9: test log-periodico
f, ax = plt.subplots(1, 2, figsize=(11, 4.2))
e_lp = ME.E_log_periodico(k_lp, NU_LP, beta=ME.BETA_PAO, amp=0.02, omega=OMEGA_INY)
m = (k_lp * eta_lp > 1e-6) & (k_lp * eta_lp < EX.x_c(0.25, beta=ME.BETA_PAO))
lk = np.log(k_lp[m])
yy = local_slope(k_lp[m], e_lp[m]) + 5 / 3
yy = yy - np.polyval(np.polyfit(lk, yy, 3), lk)
ax[0].plot(lk[::20], yy[::20], color=C1, lw=0.9)
ax[0].set_xlabel("ln k"); ax[0].set_ylabel("residuo de la pendiente")
ax[0].set_title("modulacion log-periodica (amp=0.02, $\\omega$=2)")
P = periodograma(lk, local_slope(k_lp[m], e_lp[m]) + 5 / 3, om_grid)
ax[1].plot(om_grid, P / P.max(), color=C2, lw=1.2, label="ventana larga (Re=1e10)")
ax[1].plot(om_grid, P_s / P_s.max(), color=C4, lw=1.0, alpha=0.7,
           label="ventana del barrido (Re=2000)")
ax[1].axvline(OMEGA_INY, color="k", ls="--", lw=0.9, label="$\\omega$ inyectado")
ax[1].set_xlabel("$\\omega$ (frecuencia en $\\ln k$)")
ax[1].set_ylabel("potencia normalizada")
ax[1].legend(fontsize=7); ax[1].set_title("se detecta solo si la ventana es larga")
say(fig("fig9_log_periodico.png"))

# fig10: prediccion 2D
f, ax = plt.subplots(figsize=(7, 4.6))
rr = np.logspace(2, 8, 200)
ax.semilogx(rr, EX.a_general(5 / 3, 3 / 4) * np.log(rr) - 4.13, color=C1, lw=1.6,
            label="3D K41: a = 25/12")
ax.semilogx(rr, EX.a_general(3.0, 0.5) * np.log(rr) - 8, color=C2, lw=1.6,
            label="2D enstrofia: a = 9/2")
ax.semilogx(rr, 0 * rr + 5, color=C3, lw=1.6, label="2D cascada inversa: a = 0")
ax.set_xlabel("Re"); ax.set_ylabel("$G^*[u]$"); ax.legend(fontsize=8)
ax.set_title("prediccion falsable para el paper 2D")
say(fig("fig10_prediccion_2d.png"))

with open(os.path.join(OUT_DIR, "00_salida_consola.txt"), "w") as f:
    f.write("\n".join(LOG) + "\n")
print("\nListo.")
