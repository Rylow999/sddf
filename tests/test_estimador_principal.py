#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test de regresion: el estimador principal (q_corregido) debe reproducir
q=5/3 en los 5 espectros K41 originales y mu=0.033 en los espectros con
intermitencia inyectada, siempre que se le pase el beta correcto del
modelo usado para generar el espectro (ver AUDITORIA.md, hallazgo A6).

Uso: python3 tests/test_estimador_principal.py  (desde la raiz del repo)
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "codigo"))

from estimador_principal import medir_exponente_espectral, test_contra_informe
from sddf_core import read_spectrum_csv
import sddf_exact as EX
import modelos_espectrales as ME

DATA = os.path.join(os.path.dirname(__file__), "..", "datos", "espectros")
NU = {100: 0.01, 200: 0.005, 500: 0.002, 1000: 0.001, 2000: 0.0005}


def test_k41_puro():
    filas, peor = test_contra_informe(DATA)
    assert peor < 1e-5, f"peor error {peor} supera 1e-5 en los espectros K41"
    print(f"[OK] test_k41_puro: peor error = {peor:.2e}")


def test_intermitencia_con_beta_correcto():
    peor = 0.0
    for re, nu in NU.items():
        k, E, _ = read_spectrum_csv(f"{DATA}/intermitente_Re_{re}.csv")
        eta = EX.eta_of(nu)
        q, mu, _ = medir_exponente_espectral(k, E, eta, delta=0.10,
                                              beta=ME.BETA_PAO)
        peor = max(peor, abs(mu - 0.033))
    assert peor < 1e-5, f"peor error {peor} supera 1e-5 recuperando mu=0.033"
    print(f"[OK] test_intermitencia_con_beta_correcto: peor error = {peor:.2e}")


def test_beta_mal_especificado_degrada_la_medicion():
    """Documenta (no solo advierte) que usar el beta de Pope (5.2) sobre un
    espectro generado con el beta de Pao (2.25) devuelve un mu casi nulo
    en vez de 0.033 -- la eleccion de beta importa tanto como el sesgo de
    Re finito. Ver AUDITORIA.md, hallazgo A6."""
    k, E, _ = read_spectrum_csv(f"{DATA}/intermitente_Re_1000.csv")
    eta = EX.eta_of(NU[1000])
    _, mu_beta_correcto, _ = medir_exponente_espectral(
        k, E, eta, delta=0.10, beta=ME.BETA_PAO)
    _, mu_beta_incorrecto, _ = medir_exponente_espectral(
        k, E, eta, delta=0.10, beta=ME.BETA_POPE)
    assert abs(mu_beta_correcto - 0.033) < 1e-5
    assert abs(mu_beta_incorrecto - 0.033) > 0.02, (
        "se esperaba que un beta incorrecto degradara la medicion "
        "en un orden de magnitud comparable a la senal")
    print(f"[OK] test_beta_mal_especificado: mu(beta correcto)="
          f"{mu_beta_correcto:.6f}  mu(beta de Pope, incorrecto)="
          f"{mu_beta_incorrecto:.6f}")


if __name__ == "__main__":
    test_k41_puro()
    test_intermitencia_con_beta_correcto()
    test_beta_mal_especificado_degrada_la_medicion()
    print("\nTodos los tests pasaron.")
