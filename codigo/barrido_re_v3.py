#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
barrido_re_v3.py — generador de espectros con metadatos completos.

Respecto de v2:
  * reproduce EXACTAMENTE el barrido original (modelo 'quimera', beta=5.2)
    con 15 digitos en lugar de 6, para que el error de redondeo del CSV no
    contamine la derivada logaritmica;
  * genera tambien los barridos con beta de Pao (2.25) y con el modelo
    completo de Pope, para verificar que el prefactor 25/12 no depende del
    modelo de corte;
  * escribe nu, eta, beta, C_k, modelo y x_c(delta) en la cabecera '#'.

Uso: python3 barrido_re_v3.py [dir_salida]
"""
import os
import sys

import numpy as np

import modelos_espectrales as ME

EPS, C_K = 1.0, 1.5
K = np.logspace(0, 3, 2000)
NU_VALUES = [0.01, 0.005, 0.002, 0.001, 0.0005]
RE_LABELS = [100, 200, 500, 1000, 2000]

OUT = sys.argv[1] if len(sys.argv) > 1 else "espectros"


def escribir(fn, k, E, nu, eta, modelo, beta, digits=15):
    with open(fn, "w") as f:
        f.write(f"# modelo={modelo}\n# epsilon={EPS}\n# C_k={C_K}\n# beta={beta}\n")
        re_nom = f"{1.0/nu:.0f}" if nu > 0 else "inf"
        f.write(f"# nu={nu}\n# eta={eta:.12e}\n# Re_nominal={re_nom}\n")
        f.write(f"# k_range=[{k[0]:.6e},{k[-1]:.6e}]\n# n_points={len(k)}\n")
        f.write("k,E(k),k_eta\n")
        for ki, Ei, xi in zip(k, E, k * eta):
            f.write(f"{ki:.{digits}e},{Ei:.{digits}e},{xi:.{digits}e}\n")


def main():
    os.makedirs(OUT, exist_ok=True)

    # 1) replica exacta del barrido original (6 digitos, como el archivo viejo)
    for nu, lab in zip(NU_VALUES, RE_LABELS):
        eta = ME.eta_of(nu)
        E = ME.E_quimera(K, nu, beta=ME.BETA_POPE)
        escribir(os.path.join(OUT, f"spectrum_Re_{lab}.csv"),
                 K, E, nu, eta, "quimera_pao4/3_beta_pope5.2", ME.BETA_POPE, digits=6)

    # 2) mismo barrido con beta consistente de Pao
    for nu, lab in zip(NU_VALUES, RE_LABELS):
        eta = ME.eta_of(nu)
        E = ME.E_pao(K, nu, beta=ME.BETA_PAO)
        escribir(os.path.join(OUT, f"pao_Re_{lab}.csv"),
                 K, E, nu, eta, "pao_beta2.25", ME.BETA_PAO)

    # 3) modelo completo de Pope (con f_L, corte exponencial y c_eta)
    for nu, lab in zip(NU_VALUES, RE_LABELS):
        eta = ME.eta_of(nu)
        E = ME.E_pope(K, nu, L=30.0)   # L grande: el pico de f_L queda debajo de k=1
        escribir(os.path.join(OUT, f"pope_Re_{lab}.csv"),
                 K, E, nu, eta, "pope_completo_fL_feta", ME.BETA_POPE)

    # 4) espectro de control: K41 puro sobre [1,100]
    kk = np.logspace(0, 2, 1000)
    escribir(os.path.join(OUT, "spectrum_fino.csv"),
             kk, C_K * kk ** (-5 / 3), 0.0, 0.0, "k41_puro", 0.0, digits=6)

    # 5) espectro con intermitencia real (mu=0.033, zeta_2 ~ 0.70)
    for nu, lab in zip(NU_VALUES, RE_LABELS):
        eta = ME.eta_of(nu)
        E = ME.E_pao(K, nu, beta=ME.BETA_PAO, mu=0.033)
        escribir(os.path.join(OUT, f"intermitente_Re_{lab}.csv"),
                 K, E, nu, eta, "pao_beta2.25_mu0.033", ME.BETA_PAO)

    # 6) espectro con modulacion log-periodica (test tipo Migdal)
    for nu, lab in zip(NU_VALUES, RE_LABELS):
        eta = ME.eta_of(nu)
        E = ME.E_log_periodico(K, nu, beta=ME.BETA_PAO, amp=0.02, omega=2.0)
        escribir(os.path.join(OUT, f"logper_Re_{lab}.csv"),
                 K, E, nu, eta, "pao_beta2.25_logper_amp0.02_w2", ME.BETA_PAO)

    print(f"Espectros escritos en {OUT}/")


if __name__ == "__main__":
    main()
