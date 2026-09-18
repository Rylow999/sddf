#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
modelos_espectrales.py — las distintas parametrizaciones de E(k) que aparecen
(o deberian aparecer) en el analisis de SDDF.

Punto importante que motiva este modulo
---------------------------------------
El barrido original usa

    E(k) = C_k eps^(2/3) k^(-5/3) exp[ -beta (k eta)^(4/3) ],  beta = 5.2

Eso es una QUIMERA de dos modelos distintos:

  * la forma funcional exp[-beta (k eta)^(4/3)] es la de PAO (1965), cuya
    constante NO es libre: beta = (3/2) C_k. Con C_k = 1.5 -> beta = 2.25
    (la literatura tambien usa C_k=1.6 -> beta=2.40).

  * beta = 5.2 es la constante de POPE (2000), pero corresponde a otra
    funcion de corte,
        f_eta(k eta) = exp{ -beta ( [ (k eta)^4 + c_eta^4 ]^(1/4) - c_eta ) },
        beta = 5.2, c_eta = 0.4,
    que es EXPONENCIAL en k eta (p=1), no un exponente 4/3.

Usar beta=5.2 con el exponente 4/3 corta la cola ~2.3 veces antes de lo que
corresponde y mueve el rango inercial efectivo. No invalida el resultado
principal (el prefactor 25/12 no depende de beta ni de la forma de f_eta),
pero si mueve la ordenada b(delta) y cualquier calibracion absoluta.

Este modulo permite comparar las tres versiones y verificar la universalidad.
"""
import numpy as np

C_K = 1.5
EPS = 1.0

# --- constantes de la literatura --------------------------------------
BETA_PAO = 1.5 * C_K       # = 2.25, consistente con C_k = 1.5
BETA_POPE = 5.2            # va con f_eta de Pope, p = 1 efectivo
C_ETA_POPE = 0.4
C_L_POPE = 6.78
P0_POPE = 2.0


def eta_of(nu, eps=EPS):
    return (nu ** 3 / eps) ** 0.25


# ----------------------------------------------------------------------
def E_pao(k, nu, beta=BETA_PAO, C=C_K, eps=EPS, mu=0.0):
    """Pao (1965). Con mu != 0 se le mete intermitencia en el exponente."""
    x = k * eta_of(nu, eps)
    return C * eps ** (2 / 3) * k ** (-(5 / 3 + mu)) * np.exp(-beta * x ** (4 / 3))


def E_quimera(k, nu, beta=BETA_POPE, C=C_K, eps=EPS, mu=0.0):
    """Lo que usa el barrido original: forma de Pao con beta de Pope."""
    return E_pao(k, nu, beta=beta, C=C, eps=eps, mu=mu)


def f_eta_pope(x, beta=BETA_POPE, c_eta=C_ETA_POPE):
    return np.exp(-beta * ((x ** 4 + c_eta ** 4) ** 0.25 - c_eta))


def f_L_pope(kL, c_L=C_L_POPE, p0=P0_POPE):
    return (kL / np.sqrt(kL ** 2 + c_L)) ** (5 / 3 + p0)


def E_pope(k, nu, L=1.0, C=C_K, eps=EPS, beta=BETA_POPE,
           c_eta=C_ETA_POPE, c_L=C_L_POPE, p0=P0_POPE, mu=0.0):
    """Modelo completo de Pope (Turbulent Flows, 2000), ec. 6.246."""
    x = k * eta_of(nu, eps)
    return (C * eps ** (2 / 3) * k ** (-(5 / 3 + mu))
            * f_L_pope(k * L, c_L, p0) * f_eta_pope(x, beta, c_eta))


def E_log_periodico(k, nu, beta=BETA_POPE, C=C_K, eps=EPS,
                    amp=0.02, omega=2.0, fase=0.0):
    """
    Espectro con modulacion log-periodica superpuesta al rango inercial.
    Es la firma que predice Migdal (arXiv:2511.02165, 2604.12207): las
    activaciones de la 'escalera de Stokes' generan oscilaciones periodicas
    en ln k. Sirve como banco de pruebas: si G[u] y su pendiente local no
    detectan una modulacion de amplitud 'amp', el estimador no sirve para
    contrastar la teoria de Migdal.
    """
    base = E_pao(k, nu, beta=beta, C=C, eps=eps)
    return base * (1.0 + amp * np.cos(omega * np.log(k) + fase))


def E_2d_enstrofia(k, nu, k_f=1.0, chi=1.0, C2=1.5):
    """
    Turbulencia 2D, rango de enstrofia: E ~ k^-3 con corte de Pao 2D.
    k_d2D = (chi/nu^3)^(1/6)  =>  k_d2D ~ Re^(1/2).
    Se usa para predecir que deberia dar G[u] en el paper 2D de Luciano.
    """
    k_d = (chi / nu ** 3) ** (1.0 / 6.0)
    return C2 * chi ** (2 / 3) * k ** (-3.0) * np.exp(-1.5 * C2 * (k / k_d) ** (4 / 3))


MODELOS = {
    "pao_beta2.25": lambda k, nu: E_pao(k, nu, beta=BETA_PAO),
    "quimera_beta5.2": lambda k, nu: E_quimera(k, nu, beta=BETA_POPE),
    "pope_completo": lambda k, nu: E_pope(k, nu, L=1.0),
    "pope_sin_fL": lambda k, nu: (C_K * EPS ** (2 / 3) * k ** (-5 / 3)
                                  * f_eta_pope(k * eta_of(nu))),
}
