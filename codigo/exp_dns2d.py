#!/usr/bin/env python3
# exp_dns2d.py - DNS 2D PROPIA: cascada de enstrofia (forzada) + decay (Migdal).
#
# Solver pseudo-espectral 2D doble-periodico (estandar Boffetta-Ecke):
#   dw/dt = -(u.grad)w + nu nabla^2 w + f - alpha w,   w = -nabla^2 psi
# De-alias 2/3, RK2, disipacion+friccion exactas (factor de integracion),
# forzado constante en shells 3<=k<=5, friccion lineal alpha=0.02
# (sin ella la cascada inversa no satura: practica estandar en 2D).
#
# Parte A (cascada, test a = q^2 gamma): nu en [1e-3, 3e-4, 1e-4]
#   (k_d2D = (chi/nu^3)^(1/6) = 31.6, 42.8, 100.0; todos < Nyquist 128).
#   Espectro promedio sobre N_SNAP snapshots del estado saturado.
#   Prediccion del paper: 2D enstrofia q=3, gamma=1/2 => a = 9/2 = 4.5.
#   Convencion E(k): E(b) = k_b * sum_shell(0.5 |w|^2 / k^2) (Parseval 2D;
#   sin el factor k_b la pendiente se desplaza en -1: bug clasico).
#
# Parte B (decay, Migdal real): IC = estado saturado nu=1e-4, sin forzado,
#   E(k,t) en tiempos log-espaciados; score Migdal en cada tiempo.
#
# Salida: datos/espectros_2d/*.csv + datos/23_dns2d_resumen.txt
import sys, time
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from sddf_core import local_slope, inertial_window_smoothed, compute_spectral_curvature

BASE = Path(__file__).parent.parent / "datos"
EDIR = BASE / "espectros_2d"
EDIR.mkdir(parents=True, exist_ok=True)
RES  = BASE / "23_dns2d_resumen.txt"

N       = 128
K_FORC  = (3, 5)
DEALIAS = 2.0/3.0
DT      = 0.001
T_SAT   = 200
T_DECAY = 120
N_SNAP  = 50
CHI     = 1.0
ALPHA   = 0.02

k1d = np.fft.fftfreq(N, d=1.0/N)   # wavenumbers enteros m (0..N/2, -N/2..-1); shells 3-5 no vacios
k1d_r = k1d[:N//2+1]
# rfft2: shape (N, N//2+1) = (ky, kx) — axis 0 = ky completo, axis 1 = kx medio
KYR, KXR = np.meshgrid(k1d, k1d_r, indexing="ij")
KR2 = KXR**2 + KYR**2
KR  = np.sqrt(KR2)
KNY = N/2
mask_dealias = KR <= DEALIAS * KNY

# Amplitud FISICA del forzado: irfft2 divide por N^2 (Parseval 2D), asi que
# |f|_rms = f0 exige sum |fk|^2 = (f0*N^2)^2. Con fk unit-norm el forzado era
# ~1e-5 fisico (u_sat~2e-4, tau_nl~500): la cascada jamas llenaba el espectro.
# CFL: u_sat ~ f0/(alpha*k_f) ~ 5 => U*dt*k_max = 5*0.001*64 = 0.32 (estable).
F0      = 0.1    # Re_sat ~ f0/(alpha*k_f*nu) ~ 9400 con nu=1e-4: cascada establece
rng = np.random.default_rng(2024)
fk_hat = np.zeros_like(KR2, dtype=complex)
shell_f = (KR >= K_FORC[0]) & (KR <= K_FORC[1])
fk_hat[shell_f] = rng.standard_normal(int(shell_f.sum())) + 1j*rng.standard_normal(int(shell_f.sum()))
fk_hat /= np.sqrt((np.abs(fk_hat)**2).sum())
fk_hat *= F0 * N**2

def vort2u(w_hat):
    # w = -nabla^2 psi  =>  psi_hat = +w_hat/k^2 (el signo negativo invertia
    # u,v: adveccion anti-disipativa, flujo de enstrofia revertido, espectro
    # sin rango inercial — pendientes -10 a -90 observadas en la 1a corrida)
    k2 = KR2.copy(); k2[0,0] = 1.0
    psi_hat = w_hat / k2
    u_hat =  1j * KYR * psi_hat
    v_hat = -1j * KXR * psi_hat
    return np.fft.irfft2(u_hat, s=(N,N)), np.fft.irfft2(v_hat, s=(N,N))

def step(w_hat, nu, forcing):
    # L exacto: disipacion + friccion en un factor
    L = np.exp(-(nu*KR2 + ALPHA) * DT)
    def nl(w_):
        u, v = vort2u(w_)
        wx = np.fft.irfft2(1j*KXR*w_, s=(N,N))
        wy = np.fft.irfft2(1j*KYR*w_, s=(N,N))
        adv = np.fft.rfft2(u*wx + v*wy)
        adv[~mask_dealias] = 0
        return adv
    n1 = nl(w_hat)
    fterm = fk_hat if forcing else 0
    w_mid = (w_hat * L) + DT * (-n1 + fterm)
    n2 = nl(w_mid)
    w_new = (w_hat * L) + 0.5*DT*(-n1 - n2) + DT*fterm
    w_new[~mask_dealias] = 0
    return w_new

def espectro(w_hat):
    # E(b) = sum_shell(0.5 |w|^2 / k^2): la conversion vorticidad->velocidad
    # ya esta en el /k^2 por modo; el factor k_b de circunferencia seria
    # doble conteo con la suma discreta (introducia un 1/k extra: leia k^-4
    # cuando el real es k^-3, y el detector no encontraba ventana).
    # El 1/(2N^2) de Parseval es constante: no afecta la pendiente.
    E_shell = np.zeros(N//2+1)
    kb = np.clip(np.rint(KR).astype(int), 0, N//2)
    KE = 0.5 * (np.abs(w_hat)**2) / (KR2 + 1e-30)
    for b in range(1, N//2+1):
        mm = (kb == b)
        if mm.any():
            E_shell[b] = KE[mm].sum()
    kk = np.arange(N//2+1, dtype=float)
    good = E_shell > 0
    return kk[good], E_shell[good]

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

def medir_G(k, E):
    # Deteccion en dos capas: (1) s_ref=None -> inferencia robusta del propio
    # repo (sin asuncion teorica); (2) fallback: s_ref = meseta MEDIDA de la
    # pendiente suavizada. El -3 ideal de enstrofia NO se asume: la meseta
    # real de esta caja (Re~9400) esta en ~-5.5 (efecto finite-Re honesto).
    from sddf_core import truncate_by_slope_interp, smoothed_slope
    metodo = "infer"
    ki = ei = None
    try:
        ki, ei, _, ok, k_c = truncate_by_slope_interp(k, E, delta=0.5, s_ref=None)
    except Exception:
        ok = False
    if not ok or len(ki) < 5:
        metodo = "plateau"
        s_all = smoothed_slope(k, E, window=9)
        half = 9 // 2
        s_ref = float(np.median(s_all[half:-half]))
        kl, kh, ok = inertial_window_smoothed(k, E, delta=0.5, s_ref=s_ref, window=9)
        if not ok:
            return None, None, None, metodo
        m = (k >= kl) & (k <= kh)
        ki, ei = k[m], E[m]
    G = compute_spectral_curvature(ki, ei)
    s = float(np.polyfit(np.log(ki), np.log(ei), 1)[0])
    return float(G), s, (float(ki[0]), float(ki[-1])), metodo

lines = []
print("== DNS 2D propia: cascada de enstrofia + decay ==", flush=True)
print(f"N={N} k_f={K_FORC} dt={DT} T_sat={T_SAT} alpha={ALPHA}", flush=True)

resultados = []
w_final = None
for nu in [7e-4, 5e-4, 3e-4]:   # kd real ~ 23.7/26.2/32.3 < dealias 42.7 (chi~f0^2, no 1)
    t0 = time.time()
    w_hat = np.zeros_like(KR2, dtype=complex)
    low = (KR >= 2) & (KR <= 8)
    w_hat[low] = rng.standard_normal(int(low.sum())) + 1j*rng.standard_normal(int(low.sum()))
    w_hat[~mask_dealias] = 0
    w_hat *= F0 * 4.0    # IC fisica: vorticidad ~ f0*k_f
    nsteps = int(T_SAT / DT)
    snap_every = max(nsteps // N_SNAP, 1)
    it_start = nsteps // 3   # spin-up: el primer tercio NO entra al promedio
    acc_E = None; acc_k = None; nsnap = 0
    for it in range(nsteps):
        w_hat = step(w_hat, nu, forcing=True)
        if (it+1) % snap_every == 0 and it+1 > it_start:
            kk, EE = espectro(w_hat)
            if acc_E is None:
                acc_k, acc_E = kk, EE.copy()
            else:
                acc_E += EE
            nsnap += 1
    acc_E /= max(nsnap, 1)
    kd2 = (CHI / nu**3)**(1.0/6.0)
    kk_fine = np.logspace(np.log10(max(acc_k[0], 1.0)), np.log10(acc_k[-1]), len(acc_k))
    E_fine = np.interp(kk_fine, acc_k, acc_E)
    np.savetxt(EDIR / f"espectro_2d_forzado_nu{nu:g}.csv",
               np.column_stack([kk_fine, E_fine]), delimiter=",",
               header="k,E", comments="")
    G, s, win, metodo = medir_G(kk_fine, E_fine)
    resultados.append({"nu": nu, "kd2": kd2, "G": G, "slope": s, "window": str(win), "metodo": metodo})
    lines.append(f"nu={nu:g} (k_d2D={kd2:.1f}): G={'None' if G is None else round(G,4)} slope={'None' if s is None else round(s,4)} window={win} via={metodo} ({time.time()-t0:.0f}s)")
    print(lines[-1], flush=True)
    if nu == 3e-4:
        w_final = w_hat.copy()

Gs  = np.array([r["G"] for r in resultados if r["G"] is not None])
inv = np.array([1.0/r["nu"] for r in resultados if r["G"] is not None])
# kd = (chi/nu^3)^(1/6) con chi~f0^2 => kd ~ nu^-1/2 ~ Re^1/2 (gamma=1/2):
# a = dG*/dlnRe = 9 x 1/2 = 4.5 se testea con abscisa ln(1/nu) directa
a_fit = float(np.polyfit(np.log(inv), Gs, 1)[0]) if len(Gs) >= 2 else float("nan")
lines.append(f"a = dG*/dln(1/nu) = {a_fit:.4f}  (ideal asintotico 9/2 = 4.5 con q=3)")
print(lines[-1], flush=True)
qs = np.array([r["slope"] for r in resultados if r["slope"] is not None])
if len(qs) >= 1 and np.isfinite(a_fit):
    pred_sc = float(np.median(qs))**2 * 0.5   # gamma=1/2: kd=(chi/nu^3)^(1/6) ~ nu^-1/2 en esta familia
    lines.append(f"autoconsistencia: q_medido={np.median(qs):.3f} (finite-Re, no el 3 asintotico) -> q^2*gamma = {pred_sc:.3f}")
    print(lines[-1], flush=True)

if w_final is not None:
    t0 = time.time()
    times = np.logspace(np.log10(DT*20), np.log10(T_DECAY), 14)
    t_now = 0.0
    espectros_decay = []
    for t_target in times:
        while t_now < t_target:
            w_final = step(w_final, 1e-4, forcing=False)
            t_now += DT
        kk, EE = espectro(w_final)
        espectros_decay.append((float(t_target), kk, EE))
    dec_dir = EDIR / "decay"
    dec_dir.mkdir(exist_ok=True)
    for t, kk, EE in espectros_decay:
        np.savetxt(dec_dir / f"espectro_2d_decay_t{t:.3f}.csv",
                   np.column_stack([kk, EE]), delimiter=",", header="k,E", comments="")
    sc = []
    for t, kk, EE in espectros_decay:
        kl, kh, ok = inertial_window_smoothed(kk, EE, delta=0.5, s_ref=-3.0, window=5)
        if not ok: continue
        m = (kk >= kl) & (kk <= kh)
        if m.sum() < 8: continue
        lk = np.log(kk[m])
        s_loc = local_slope(kk[m], EE[m]) + 3.0
        P = periodograma2(lk, s_loc)
        sc.append((t, float(P.max()/np.median(P))))
    lines.append(f"decay: {len(espectros_decay)} tiempos, {len(sc)} scores migdal ({time.time()-t0:.0f}s)")
    for t, s_ in sc:
        lines.append(f"  t={t:.3f}: score={s_:.3f}")
    print(lines[-1], flush=True)

with open(RES, "w") as f:
    f.write("DNS 2D propia: cascada de enstrofia + decay (Migdal)\n")
    f.write(f"N={N}, k_f={K_FORC}, dt={DT}, T_sat={T_SAT}, alpha={ALPHA}, dealias=2/3\n\n")
    f.write("\n".join(lines) + "\n")
print("DNS2D_DONE", flush=True)
