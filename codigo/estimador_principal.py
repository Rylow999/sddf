#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
estimador_principal.py — el estimador recomendado del exponente espectral,
promovido desde el bloque 9 de sddf_completo.py (informe v2, seccion 4.3)
a una funcion de API reutilizable, documentada, y con sus propios tests.

Por que este y no <s^2> = G*/ln(k_c/k_min)
-------------------------------------------
<s^2> converge al exponente asintotico solo como 1/ln(Re): con delta=0.10 y
Re=100 el sesgo implicito es mu_eff=0.084, cuatro veces la intermitencia
fisica que se querria medir (ver informe v2, seccion 4.3).

q_corregido() en cambio despeja q de la forma cerrada exacta

    G = q^2 ln(x_c/x1) + 2 q beta (x_c^p - x1^p) + (beta^2 p/2)(x_c^2p - x1^2p)

de forma iterativa (los dos ultimos terminos dependen de q solo a traves del
termino lineal, asi que 3-5 iteraciones alcanzan para 6 cifras). Con esto:

  * funciona con un UNICO espectro — no hace falta el barrido en Re;
  * en los 5 espectros originales da q = 1.666667 +/- 0.000001 en todos los
    Re y ambos delta probados (ver test_contra_informe() mas abajo);
  * es la pieza que el informe v2 (seccion 11, punto 2) marca como "el
    estimador principal de ahora en mas".

Uso tipico sobre un espectro real (DNS, experimento, lo que sea):

    from estimador_principal import medir_exponente_espectral
    q, mu, k_corte = medir_exponente_espectral(k, E, eta, delta=0.10)
"""
import numpy as np

from sddf_core import (read_spectrum_csv, compute_spectral_curvature,
                        truncate_by_slope_interp, g_normalizado)
import sddf_exact as EX

BETA_DEFAULT = EX.BETA_DEFAULT
P_DEFAULT = EX.P_DEFAULT
Q_K41 = EX.Q_K41


def q_corregido(G, k_min, k_c, eta, delta, beta=BETA_DEFAULT, p=P_DEFAULT,
                 q0=Q_K41, iters=5):
    """
    Despeja q de la forma cerrada en vez de dividir G/ln(k_c/k_min) a lo
    bruto. Los terminos de corte (t2, t3) son conocidos una vez fijados
    delta y beta; q converge en pocas iteraciones porque solo entra
    linealmente en t2. Ver informe v2, seccion 4.3.
    """
    x1, xc_ = k_min * eta, k_c * eta
    L = np.log(xc_ / x1)
    t3 = (beta ** 2 * p / 2.0) * (xc_ ** (2 * p) - x1 ** (2 * p))
    q = q0
    for _ in range(iters):
        t2 = 2.0 * q * beta * (xc_ ** p - x1 ** p)
        q = np.sqrt(max(G - t2 - t3, 0.0) / L)
    return q


def medir_exponente_espectral(k, E, eta, delta=0.10, beta=BETA_DEFAULT,
                               p=P_DEFAULT, q0=Q_K41, iters=5):
    """
    API de una sola llamada: de un espectro (k, E) y su escala de Kolmogorov
    eta, devuelve (q_corregido, mu_corregido, k_corte).

    mu_corregido = q_corregido - 5/3 es la correccion de intermitencia del
    exponente espectral, comparable directamente con zeta_2/2 - algo o con
    otras estimaciones de intermitencia de la literatura (informe v2, 4.2).
    """
    k = np.asarray(k, float)
    E = np.asarray(E, float)
    ki, ei, _, ok, k_c = truncate_by_slope_interp(k, E, delta=delta,
                                                   s_ref=-q0)
    if not ok:
        raise ValueError("ventana inercial insuficiente para este delta: "
                          "el corte cayo en los primeros 3 puntos")
    G = compute_spectral_curvature(ki, ei)
    q = q_corregido(G, ki[0], k_c, eta, delta, beta=beta, p=p, q0=q0,
                     iters=iters)
    return float(q), float(q - Q_K41), float(k_c)


# ----------------------------------------------------------------------
# Test de regresion contra el informe v2 (seccion 4.3): debe reproducir
# q = 1.666667 +/- 0.000001 en los 5 espectros originales, ambos delta.
# ----------------------------------------------------------------------
def test_contra_informe(dir_espectros="../datos/espectros"):
    NU = {100: 0.01, 200: 0.005, 500: 0.002, 1000: 0.001, 2000: 0.0005}
    peor_error = 0.0
    filas = []
    for delta in (0.10, 0.25):
        for re, nu in NU.items():
            k, E, _ = read_spectrum_csv(f"{dir_espectros}/spectrum_Re_{re}.csv")
            eta = EX.eta_of(nu)
            q, mu, k_c = medir_exponente_espectral(k, E, eta, delta=delta)
            err = abs(q - 5.0 / 3.0)
            peor_error = max(peor_error, err)
            filas.append((delta, re, q, mu, err))
    return filas, peor_error


if __name__ == "__main__":
    filas, peor = test_contra_informe()
    print(f"{'delta':>6} {'Re':>6} {'q_corregido':>13} {'mu_corregido':>13} "
          f"{'|error|':>10}")
    for delta, re, q, mu, err in filas:
        print(f"{delta:6.2f} {re:6d} {q:13.6f} {mu:13.6f} {err:10.2e}")
    print(f"\npeor error absoluto respecto de q=5/3: {peor:.2e}")
    assert peor < 1e-5, "el estimador deberia reproducir q=5/3 a 1e-5"
    print("OK: reproduce el informe v2, seccion 4.3.")
