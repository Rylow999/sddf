# Changelog

## v3.2 (2026-09-22) — null model Migdal + sync Rust + DNS real

### Added
- `codigo/exp_null_model_periodograma.py`: test de calibración con 500
  nulos instrumentales (beta/nu jittereados, ruido 1e-4), p-empíricos,
  Bonferroni sobre las 5 amplitudes y umbral global `max(null)`. Salida
  `datos/16_null_model_periodograma.csv`. **Confirma que el SNR~19
  reportado sin señal inyectada era artefacto del método** (max de los
  nulos = 19.63) y fija el umbral real de detección entre amp=0.005 y
  amp=0.02 (con la ventana inercial de Re=1e10 usada en bloque 14).
- `codigo/exp_dns_real_isotropo.py`: validación con DNS real. Baja el
  snapshot público `isotropic1024coarse` (256³, Re_λ≈433) de JHTDB desde el
  mirror HuggingFace, calcula el espectro isotropo por FFT 3D + conchas
  |k| y corre el pipeline SDDF. Salidas:
  `datos/espectros/spectrum_jhtdb_isotropic1024coarse_t<t>.csv` y
  `datos/17_dns_real_resumen.txt`.

### Fixed / Changed
- `codigo/rust/src/main.rs`: el corte por pendiente ahora usa
  **interpolación lineal** en (ln k, s) — idéntica a
  `sddf_core.truncate_by_slope_interp` (mismo k_corte, mismo G* al dígito
  sobre `spectrum_Re_1000.csv`). El binario pasa a v3.1.0.
- Fallback sintético del binario Rust: Pao con **β=2.25** (quimera
  acordada) en lugar de β=5.2 (Pope) — consistente con la advertencia A3/A6
  de `AUDITORIA.md`.
- `codigo/sddf_completo.py`: encabezado de bloques amplía a 16/17 (null
  model y DNS real).

### Recuerdos / pendientes abiertos
- Ventana larga en `exp_null_model` usa 50k puntos y 1500 frecuencias
  (bajado de 200k/4000 por costo) — es autosuficiente y numéricamente
  indistinguible del set original (error relativo del periodograma
  vectorizado <2e-14).
- Falta promediar los 10 timesteps del isotropic1024coarse y bajar
  channel4094 para tener inercial largo en DNS real.

## v3.1 (2026-09-17) — fixes post-auditoría de pares

Cambios propuestos por auditoría externa y validados con tests.

### Added

- `sddf_core.q_from_G_closed_form(G, k1, kc, beta, p)`: despeje **analítico
  directo** de q desde la forma cerrada (cuadrática explícita). Reemplaza la
  iteración de punto fijo; **lanza ValueError si el discriminante es
  negativo** (datos inconsistentes) en vez de devolver 0 en silencio.
- `sddf_core.inertial_window(k, e, delta, s_ref)`: detector de rango inercial
  de **dos lados** (k_low + k_high). El código original solo cortaba el lado
  de disipación; el de forzado quedaba invisibile.

### Fixed

- Documentación de `truncate_by_slope_interp` aclarada: el corte es contra
  `s_ref` (referencia K41 por defecto), **no** contra el q que se está
  estimando. La anterior ambigüedad era correcta matemáticamente pero
  peligrosa semánticamente.
- Figura 8 (comparación de modelos): el ratio de corte era ~1.87x, no 2.3x
  (título corregido).

### Deprecated / known issues

- El periodograma log-periódico de Migdal es débil: usa
  `peak-to-median` como score, no es un test estadístico calibrado. Queda
  pendiente un null-model empírico (ver `CONTRASTE_MIGDAL_DNS.md`).
- La implementación Rust queda atrasada respecto a Python: no tiene
  interpolación de corte ni el fallback beta correcto.

### Verificado

- La fórmula cerrada reproduce los 5 espectros K41: error relativo < 1e-13.
- El test de regresión pasa (`tests/test_estimador_principal.py`).
