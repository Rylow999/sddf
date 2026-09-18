#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sddf_core.py — nucleo numerico de la curvatura espectral

    G[u] = INT [ d(ln E)/d(ln k) ]^2  d(ln k)

Replica bit a bit el esquema de main.rs (derivadas centradas en el interior,
forward/backward en los bordes, trapecio sobre ln k) y agrega:

  * truncate_by_slope        corte por desviacion de pendiente (v1)
  * truncate_by_slope_interp corte INTERPOLADO en ln k (v2, recomendado):
                             elimina el error de cuantizacion del corte,
                             que es la fuente dominante de error en grillas
                             de pocos miles de puntos
  * truncate_by_x            corte directo en k*eta = x_c (si se conoce nu)
  * g_normalizado            G* / ln(k_c/k_min) = <s^2>, invariante de Re
  * q_efectivo               exponente espectral implicado por <s^2>

El piso absoluto 1e-30 del lector original se conserva SOLO como opcion
(floor=1e-30) para poder reproducir el bug historico.
"""
import numpy as np

S_K41 = -5.0 / 3.0


# ----------------------------------------------------------------------
# Lectura
# ----------------------------------------------------------------------
def read_spectrum_csv(path, floor=None):
    """
    Lee un CSV 'k,E(k)[,k_eta]'. Ignora comentarios '#' y la cabecera.
    floor=1e-30 reproduce el piso del lector original (bug historico).
    floor=None descarta los puntos con E<=0 o no finito, que es lo correcto.
    Devuelve (k, E, n_descartados).
    """
    k_vec, e_vec, dropped = [], [], 0
    with open(path, "r") as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split(",")
            if len(parts) < 2:
                continue
            try:
                k = float(parts[0].strip())
                e = float(parts[1].strip())
            except ValueError:
                continue  # cabecera
            if floor is not None:
                k_vec.append(k)
                e_vec.append(e if e > floor else floor)
            else:
                if np.isfinite(k) and k > 0 and np.isfinite(e) and e > 0:
                    k_vec.append(k)
                    e_vec.append(e)
                else:
                    dropped += 1
    return np.asarray(k_vec, float), np.asarray(e_vec, float), dropped


def read_meta(path):
    """Extrae los '# clave=valor' de la cabecera (nu, eta, beta, C_k, ...)."""
    meta = {}
    with open(path, "r") as fh:
        for line in fh:
            s = line.strip()
            if not s.startswith("#"):
                if s and not s.startswith("#"):
                    break
                continue
            if "=" in s:
                kk, vv = s[1:].split("=", 1)
                meta[kk.strip()] = vv.strip()
    return meta


# ----------------------------------------------------------------------
# Nucleo
# ----------------------------------------------------------------------
def local_slope(k, e):
    """Pendiente log-log local. Identico esquema que main.rs."""
    lk, le = np.log(np.asarray(k, float)), np.log(np.asarray(e, float))
    n = len(lk)
    d = np.zeros(n)
    if n < 2:
        return d
    d[1:-1] = (le[2:] - le[:-2]) / (lk[2:] - lk[:-2])
    d[0] = (le[1] - le[0]) / (lk[1] - lk[0])
    d[-1] = (le[-1] - le[-2]) / (lk[-1] - lk[-2])
    return d


def compute_spectral_curvature(k, e):
    """G[u] por trapecio sobre ln k. Identico a main.rs."""
    k = np.asarray(k, float)
    if len(k) < 3:
        return 0.0
    lk = np.log(k)
    d = local_slope(k, e)
    return float(np.trapezoid(d ** 2, lk))


# ----------------------------------------------------------------------
# Truncados
# ----------------------------------------------------------------------
def truncate_tail(k, e, rel_thresh=1e-6):
    """Corte por amplitud relativa. DEPRECADO: el exponente que se ajuste
    despues depende del umbral (ver 03_sensibilidad_umbral.csv)."""
    k, e = np.asarray(k, float), np.asarray(e, float)
    if len(k) < 3:
        return k, e, 0
    below = np.nonzero(e < rel_thresh * e[0])[0]
    if len(below) == 0:
        return k, e, 0
    cut = max(int(below[0]), 3)
    return k[:cut], e[:cut], len(k) - cut


def truncate_by_slope(k, e, delta=0.5, s_ref=S_K41):
    """Corte por desviacion de pendiente, en el primer punto de grilla."""
    k, e = np.asarray(k, float), np.asarray(e, float)
    if len(k) < 3:
        return k, e, 0, False
    s = local_slope(k, e)
    bad = np.nonzero(np.abs(s - s_ref) > delta)[0]
    if len(bad) == 0:
        return k, e, 0, True
    cut = int(bad[0])
    ok = cut >= 3          # bandera: si es False, el resultado NO significa nada
    cut = max(cut, 3)
    return k[:cut], e[:cut], len(k) - cut, ok


def truncate_by_slope_interp(k, e, delta=0.5, s_ref=S_K41):
    """
    Corte por desviación de pendiente, USANDO s_ref (no hardcode -5/3).

    El parámetro s_ref es el estado de referencia del observador.
    Si se pasa None, intenta inferir el rango inercial robusto y luego corta.
    Por defecto se interpreta K41 como en la publicación original.

    NOTA: Este código usa truncado por umbral + interpolación lineal en (ln k, s).
    En casos patológicos (ruido fuerte), no olvides considerar alternativas:
    - correlate smoother primitives
    - robust lowess
    - hacérselo previo el análisis (entrenar un clasificador simple) antes de cortar.
    """
    k, e = np.asarray(k, float), np.asarray(e, float)
    if len(k) < 3:
        return k, e, 0, False, np.nan
    s = local_slope(k, e)
    dev = np.abs(s - s_ref)
    bad = np.nonzero(dev > delta)[0]
    if len(bad) == 0:
        return k, e, 0, True, float(k[-1])
    j = int(bad[0])
    if j < 3:
        return k[:3], e[:3], len(k) - 3, False, float(k[2])
    lk = np.log(k)
    d0, d1 = dev[j - 1], dev[j]
    w = 0.0 if d1 == d0 else (delta - d0) / (d1 - d0)
    w = min(max(w, 0.0), 1.0)
    lkc = lk[j - 1] + w * (lk[j] - lk[j - 1])
    lec = np.log(e[j - 1]) + w * (np.log(e[j]) - np.log(e[j - 1]))
    k_new = np.concatenate([k[:j], [np.exp(lkc)]])
    e_new = np.concatenate([e[:j], [np.exp(lec)]])
    return k_new, e_new, len(k) - j, True, float(np.exp(lkc))


def truncate_by_x(k, e, eta, x_c):
    """Corte directo en k*eta = x_c (requiere conocer nu). Interpolado."""
    k, e = np.asarray(k, float), np.asarray(e, float)
    k_c = x_c / eta
    if k_c >= k[-1]:
        return k, e, 0, True, float(k[-1])
    j = int(np.searchsorted(k, k_c))
    if j < 3:
        return k[:3], e[:3], len(k) - 3, False, float(k[2])
    lk = np.log(k)
    w = (np.log(k_c) - lk[j - 1]) / (lk[j] - lk[j - 1])
    lec = np.log(e[j - 1]) + w * (np.log(e[j]) - np.log(e[j - 1]))
    k_new = np.concatenate([k[:j], [k_c]])
    e_new = np.concatenate([e[:j], [np.exp(lec)]])
    return k_new, e_new, len(k) - j, True, float(k_c)


# ----------------------------------------------------------------------
# Estimadores derivados
# ----------------------------------------------------------------------
def g_normalizado(k, e):
    """
    <s^2> = G / ln(k_max/k_min). Invariante de Re: no cuenta decadas,
    mide la forma del espectro. Para K41 puro da 25/9 = 2.7778.
    """
    k = np.asarray(k, float)
    if len(k) < 3:
        return np.nan
    span = np.log(k[-1] / k[0])
    if span <= 0:
        return np.nan
    return compute_spectral_curvature(k, e) / span


def q_efectivo(k, e):
    """Exponente espectral efectivo -q implicado por <s^2>: q = sqrt(<s^2>)."""
    gn = g_normalizado(k, e)
    return np.nan if not np.isfinite(gn) else np.sqrt(gn)


def g_original_con_bug(path):
    """G[u] tal como lo daba el lector original (piso absoluto 1e-30)."""
    k, e, _ = read_spectrum_csv(path, floor=1e-30)
    return compute_spectral_curvature(k, e), int((e <= 1e-30).sum())


# ---------------------------------------------------------------------------
# Estimador analítico directo de q (sin iteración de punto fijo)
# ---------------------------------------------------------------------------

def q_from_G_closed_form(G, k1, kc, beta, p):
    """
    Despeja q de la forma cerrada de G.

    La fórmula cerrada es:
        G = q^2 * L + q * A + B
    donde:
        L = ln(kc / k1)
        A = 2 * beta * (kc^p - k1^p)
        B = (beta^2 * p / 2) * (kc^(2p) - k1^(2p))

    Resolviendo la cuadrática en q:
        L*q^2 + A*q + (B - G) = 0
        q = (-A + sqrt(A^2 - 4*L*(B - G))) / (2*L)

    Preferir la solución positiva (q es una magnitud).

    Raises:
        ValueError: si el discriminante es negativo (los datos son
            inconsistentes con el modelo espectral) — esto es PREFERIBLE
            al comportamiento anterior, que devolvía q=0 en silencio.
    """
    L = np.log(kc / k1)
    if L <= 0:
        raise ValueError(f"L = ln(kc/k1) debe ser > 0, recibido L={L}")
    A = 2.0 * beta * (kc ** p - k1 ** p)
    B = (beta ** 2 * p / 2.0) * (kc ** (2 * p) - k1 ** (2 * p))
    discrim = A ** 2 - 4.0 * L * (B - G)
    if discrim < 0:
        raise ValueError(
            f"Datos inconsistentes con el modelo: discriminante={discrim:.4e} < 0. "
            f"Verifique la ventana [k1={k1:.3g}, kc={kc:.3g}] y el modelo (beta, p)."
        )
    return (-A + np.sqrt(discrim)) / (2.0 * L)


# ---------------------------------------------------------------------------
# Detector de rango inercial de DOS lados (P1 del informe de auditoría)
# ---------------------------------------------------------------------------

def inertial_window(k, e, delta=0.5, s_ref=S_K41):
    """
    Detecta AMBOS extremos del rango inercial:

    - lado bajo (k_low): primer k donde |s(k) - s_ref| < delta viniendo desde k[0]
      (el extremo de forzado/k grandes escalas dejan de distorsionar)
    - lado alto (k_high): último k donde |s(k) - s_ref| < delta antes de que la
      disipación rompa el inercial (equivale al truncate_by_slope_interp)

    Devuelve (k_low, k_high, ok). Si ok=False, la ventana no se pudo separar
    limpiamente (el espectro probablemente está dominado por no-inercial).
    """
    k = np.asarray(k, float)
    e = np.asarray(e, float)
    if len(k) < 6:
        return k[0], k[-1], False
    s = local_slope(k, e)
    good = np.abs(s - s_ref) <= delta
    # k_high: el primer punto después del cual TODO está mal (región disipación)
    # Estrategia: desde la izquierda, el último 'good' antes de que
    # empiece una bajada sostenida hacia la disipación (promedio de 5 malos
    # seguidos en adelante)
    k_high = k[-1]
    for i in range(len(k) - 4):
        if not np.all(good[i:i + 5]):
            k_high = k[i]
            break
    # k_low: primer 'good' desde la izquierda
    k_low = k[0]
    for i in range(len(k)):
        if good[i]:
            k_low = k[i]
            break
    ok = k_high > k_low * 3  # pedimos al menos 1 década de rango inercial
    return float(k_low), float(k_high), ok
