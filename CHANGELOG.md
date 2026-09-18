# Changelog

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
