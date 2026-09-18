#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sddf_exact.py — formas cerradas exactas para la curvatura espectral

    G[u] = INT [ d(ln E)/d(ln k) ]^2  d(ln k)

sobre espectros del tipo Pao/Pope generalizado

    E(k) = C eps^(2/3) k^(-q) exp[ -beta (k eta)^p ],    eta = (nu^3/eps)^(1/4)

con q = 5/3 + mu  (mu = correccion de intermitencia del exponente espectral)
y   p = 4/3       (forma de Pao; p=1 seria la forma exponencial pura).

Todo lo que sigue es analitico: no hay ajuste ni integracion numerica.

--------------------------------------------------------------------------
RESUMEN DE LOS RESULTADOS (demostraciones en INFORME.md, seccion 3)

  pendiente local      s(x) = -q - beta p x^p,          x = k eta

  G(x1,x2) = q^2 ln(x2/x1)
             + 2 q beta (x2^p - x1^p)
             + (beta^2 p / 2) (x2^(2p) - x1^(2p))

  corte por pendiente  |s + q| = delta   <=>   x_c = ( delta / (beta p) )^(1/p)

  ventana fisica [k_min, k_c] con k_min fijo, eps = 1, Re = 1/nu:

     x1 = k_min Re^(-3/4)

     G*(Re,delta) = (3/4) q^2 ln Re
                    + b(delta)
                    - 2 q beta k_min^p Re^(-3p/4)
                    - (beta^2 p/2) k_min^(2p) Re^(-3p/2)

     b(delta) = (q^2/p) ln[ delta/(beta p) ] - q^2 ln k_min
                + 2 q delta / p + delta^2/(2 p)

  Para q=5/3, p=4/3, k_min=1 esto colapsa a la forma compacta

     G*(Re,delta) = (25/12) ln[ 3 delta Re / (4 beta) ]
                    + (5/2) delta + (3/8) delta^2
                    - (10 beta/3) / Re
                    - (2 beta^2/3) / Re^2

  => PREFACTOR ASINTOTICO  a = (3/4) q^2 ,  que vale 25/12 = 2.08333 en K41
     y NO depende de beta, de delta, de C ni de k_min.
--------------------------------------------------------------------------
"""
import numpy as np

# Parametros del barrido original de Luciano
BETA_DEFAULT = 5.2      # valor de Pope para f_eta exponencial (ver INFORME sec. 6)
P_DEFAULT = 4.0 / 3.0   # exponente de Pao
Q_K41 = 5.0 / 3.0
A_ASYMP_K41 = 25.0 / 12.0


# ----------------------------------------------------------------------
# Bloque 1: pendiente, corte, G
# ----------------------------------------------------------------------
def slope_x(x, q=Q_K41, beta=BETA_DEFAULT, p=P_DEFAULT):
    """Pendiente log-log exacta s = dlnE/dlnk en funcion de x = k*eta."""
    return -q - beta * p * np.asarray(x, float) ** p


def x_c(delta, beta=BETA_DEFAULT, p=P_DEFAULT):
    """x=k*eta donde |s+q| = delta. Invariante de Re por construccion."""
    return (delta / (beta * p)) ** (1.0 / p)


def G_exact(x1, x2, q=Q_K41, beta=BETA_DEFAULT, p=P_DEFAULT):
    """G entre dos valores de x = k*eta. Forma cerrada, sin cuadratura."""
    x1 = np.asarray(x1, float)
    x2 = np.asarray(x2, float)
    return (q * q * np.log(x2 / x1)
            + 2.0 * q * beta * (x2 ** p - x1 ** p)
            + (beta ** 2 * p / 2.0) * (x2 ** (2 * p) - x1 ** (2 * p)))


def eta_of(nu, eps=1.0):
    return (nu ** 3 / eps) ** 0.25


# ----------------------------------------------------------------------
# Bloque 2: G* con ventana fisica, descompuesto termino a termino
# ----------------------------------------------------------------------
def b_delta(delta, q=Q_K41, beta=BETA_DEFAULT, p=P_DEFAULT, k_min=1.0):
    """Ordenada al origen exacta de la ley logaritmica."""
    return ((q * q / p) * np.log(delta / (beta * p))
            - q * q * np.log(k_min)
            + 2.0 * q * delta / p
            + delta ** 2 / (2.0 * p))


def a_asymptotic(q=Q_K41):
    """Prefactor asintotico dG*/dlnRe = (3/4) q^2."""
    return 0.75 * q * q


def finite_Re_terms(Re, q=Q_K41, beta=BETA_DEFAULT, p=P_DEFAULT, k_min=1.0):
    """
    Los dos terminos de correccion de Re finito, por separado y con signo:
    G*_exacto = a*lnRe + b(delta) + t1 + t2  con t1,t2 <= 0.
    """
    Re = np.asarray(Re, float)
    t1 = -2.0 * q * beta * k_min ** p * Re ** (-0.75 * p)
    t2 = -(beta ** 2 * p / 2.0) * k_min ** (2 * p) * Re ** (-1.5 * p)
    return t1, t2


def G_star(Re, delta, q=Q_K41, beta=BETA_DEFAULT, p=P_DEFAULT, k_min=1.0):
    """G*(Re,delta) exacto: ley logaritmica + correcciones de Re finito."""
    t1, t2 = finite_Re_terms(Re, q, beta, p, k_min)
    return (a_asymptotic(q) * np.log(np.asarray(Re, float))
            + b_delta(delta, q, beta, p, k_min) + t1 + t2)


def G_star_compact(Re, delta, beta=BETA_DEFAULT):
    """Version compacta valida solo para q=5/3, p=4/3, k_min=1."""
    Re = np.asarray(Re, float)
    return ((25.0 / 12.0) * np.log(Re * 3.0 * delta / (4.0 * beta))
            + 2.5 * delta + 0.375 * delta ** 2
            - (10.0 * beta / 3.0) / Re
            - (2.0 * beta ** 2 / 3.0) / Re ** 2)


# ----------------------------------------------------------------------
# Bloque 3: estimadores del exponente espectral a partir del prefactor
# ----------------------------------------------------------------------
def q_from_a(a):
    """Invierte a = (3/4) q^2  ->  q = sqrt(4a/3). Devuelve el exponente -q."""
    return np.sqrt(4.0 * np.asarray(a, float) / 3.0)


def mu_from_a(a):
    """Correccion de intermitencia implicada por un prefactor medido."""
    return q_from_a(a) - Q_K41


def a_from_mu(mu):
    """Prefactor que corresponde a una intermitencia dada."""
    return 0.75 * (Q_K41 + np.asarray(mu, float)) ** 2


def deviscosify(Re, G, beta=BETA_DEFAULT, p=P_DEFAULT, q=Q_K41, k_min=1.0):
    """
    Quita analiticamente las correcciones de Re finito de una serie G*(Re).
    El ajuste lineal de la salida contra ln Re da (3/4) q^2 SIN sesgo,
    incluso con ventanas de Re cortas.
    """
    t1, t2 = finite_Re_terms(Re, q, beta, p, k_min)
    return np.asarray(G, float) - t1 - t2


def fit_log(Re, G):
    """Ajuste G = a ln Re + b. Devuelve (a, b, R2)."""
    x = np.log(np.asarray(Re, float))
    y = np.asarray(G, float)
    a, b = np.polyfit(x, y, 1)
    r2 = 1.0 - ((y - (a * x + b)) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    return a, b, r2


def fit_power(Re, G):
    """Ajuste G = A Re^(-alpha). Devuelve (A, alpha, R2)."""
    x = np.log(np.asarray(Re, float))
    y = np.log(np.asarray(G, float))
    s, i = np.polyfit(x, y, 1)
    r2 = 1.0 - ((y - (s * x + i)) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    return np.exp(i), -s, r2


# ----------------------------------------------------------------------
# Bloque 4: prediccion para turbulencia 2D
# ----------------------------------------------------------------------
def a_2d_enstrofia():
    """
    Cascada de enstrofia 2D: E ~ k^-3, escala de corte k_d ~ (chi/nu^3)^(1/6)
    => ln k_d = (1/2) ln Re + cte  =>  dG*/dlnRe = 9 * (1/2) = 4.5
    """
    return 9.0 * 0.5


def a_2d_inversa():
    """
    Cascada inversa 2D: E ~ k^-5/3 acotada por el tamano de caja, no por nu
    => G* no crece con Re en ese rango.
    """
    return 0.0


def a_general(exponente_espectral, exponente_eta):
    """
    Caso general: E ~ k^-q en el rango resuelto y k_corte ~ Re^(exponente_eta).
    dG*/dlnRe = q^2 * exponente_eta.
    3D K41: q=5/3, exponente_eta=3/4 -> 25/12.
    2D enstrofia: q=3, exponente_eta=1/2 -> 4.5.
    """
    return exponente_espectral ** 2 * exponente_eta


if __name__ == "__main__":
    Re = np.array([100, 200, 500, 1000, 2000], float)
    for d in (0.1, 0.25, 0.5, 1.0):
        g = G_star(Re, d)
        gc = G_star_compact(Re, d)
        assert np.allclose(g, gc), "las dos formas cerradas deben coincidir"
        a, b, r2 = fit_log(Re, g)
        ad, bd, _ = fit_log(Re, deviscosify(Re, g))
        print(f"delta={d:4} a_cruda={a:.6f} a_desviscosada={ad:.6f} "
              f"b_exacto={b_delta(d):.6f} mu_eff={mu_from_a(a):.5f}")
    print(f"\nasintota 25/12 = {A_ASYMP_K41:.6f}")
