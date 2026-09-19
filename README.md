# SDDF — Curvatura espectral G[u]

Extensión de **LOGOS**, dentro de [`vega-vault`](https://github.com/Rylow999/vega-vault).

Estudio de la curvatura espectral

$$G[u] = \int \left(\frac{d\ln E(k)}{d\ln k}\right)^2 d\ln k$$

como diagnóstico del rango inercial de Kolmogorov, con forma cerrada
exacta, un estimador de exponente espectral corregido de sesgos, y un
test de la firma log-periódica propuesta por Migdal (arXiv:2604.12207,
arXiv:2511.02165) para turbulencia decayente.

**Resultado central:** $G^*[u] = \frac{25}{12}\ln Re + b(\delta) +
\mathcal{O}(Re^{-1})$, con prefactor $25/12=(25/9)\cdot(3/4)$ — el cuadrado
de la pendiente de Kolmogorov por el exponente de escala de $\eta\propto
Re^{-3/4}$. No es un ajuste.

## Estructura del repositorio

```
.
├── README.md              este archivo
├── AUDITORIA.md            historial completo de las rondas de revisión (v1→v2→v3→A7)
├── CONTRASTE_MIGDAL_DNS.md guía completa para contrastar la firma de Migdal en DNS real
├── paper/
│   └── PAPER.md            documento consolidado con todos los resultados
├── codigo/
│   ├── sddf_core.py         nucleo numerico (lectura, pendiente local, truncados)
│   ├── sddf_exact.py         formas cerradas analiticas
│   ├── estimador_principal.py  estimador recomendado (q_corregido), API reutilizable
│   ├── modelos_espectrales.py  parametrizaciones Pao / Pope / quimera / 2D
│   ├── barrido_re_v3.py       genera los espectros CSV de datos/espectros/
│   ├── sddf_completo.py       pipeline completo (regenera todos los CSV/PNG)
│   └── rust/                  implementacion de referencia en Rust (compilada)
├── datos/                  15 tablas de resultados + 27 espectros CSV (incluye spectrum_grueso.csv, recuperado)
├── figuras/                 10 figuras (PNG)
└── tests/
    └── test_estimador_principal.py   regresion contra el paper
```

## Reproducir todo desde cero

```bash
cd codigo
python3 barrido_re_v3.py ../datos/espectros
python3 sddf_completo.py ../datos/espectros ../datos ../figuras
./rust/bin/sddf ../datos/espectros/spectrum_Re_1000.csv 0.5
cd .. && python3 tests/test_estimador_principal.py
```

## Usar el estimador sobre un espectro propio

```python
from estimador_principal import medir_exponente_espectral

# k, E: arrays de tu espectro. eta: escala de Kolmogorov (nu^3/eps)^(1/4).
# beta: constante del modelo de corte que corresponda a TU espectro
# (no asumas 5.2 de Pope sin verificar — ver AUDITORIA.md, hallazgo A3/A6).
q, mu, k_corte = medir_exponente_espectral(k, E, eta, delta=0.10, beta=2.25)
```

## Estado

Ver `AUDITORIA.md` para el historial completo y `paper/PAPER.md` sección 11
para las limitaciones abiertas (recuperación de `spectrum_grueso.csv` y
recálculo del paper 2D, ambos bloqueados por falta de datos de entrada, no
de método).

## Novedades v3.1 (2026-09-17) — fixes post-auditoría

Cambios derivados de la auditoría externa, validados con tests
(ver `CHANGELOG.md` para el detalle completo):

- **`sddf_core.q_from_G_closed_form`** — despeje analítico directo de q
  desde la forma cerrada (cuadrática explícita, error ~1e-14). Reemplaza la
  iteración de punto fijo. **Lanza ValueError si el discriminante es
  negativo** (datos inconsistentes con el modelo) en vez de devolver 0 en
  silencio.
- **`sddf_core.inertial_window`** — detector de rango inercial de **dos
  lados** (k_low + k_high). El pipeline original solo cortaba el lado de
  disipación; el lado de forzado (k bajos) quedaba invisible y contaminaba
  el integral.
- Semántica de `truncate_by_slope_interp` documentada: el corte es contra
  `s_ref` (K41 por defecto), no contra el q que se está estimando.

### Pendiente (honesto)

- **Null model del periodograma de Migdal**: el score actual
  (`P.max()/median(P)`) no es significancia estadística — fabrica picos con
  "SNR"~19 sin señal inyectada. El null model (500 espectros sin modulación +
  p-values empíricos + Bonferroni) está en
  `NOUS/RHO_LAW/experiments/exp_null_model_periodograma.py`.
- **Sincronizar Rust con Python**: la implementación Rust no tiene
  interpolación de corte ni el fallback beta correcto (quimera). Deprecada
  hasta sincronizar.
- **Validación con DNS real**: todos los espectros actuales son sintéticos.
  Siguiente paso: JHTDB (Johns Hopkins Turbulence Database, públicos).

## Conexión con la ley ρ (RHO_LAW)

El observable G[u] comparte el patrón del "colapso del observador" con la
ley ρ del repo [`Rylow999/fhrr-rho-collapse`](https://github.com/Rylow999/fhrr-rho-collapse):

- Cuando el espectro se mide con grilla discreta, el error de truncamiento
  **no decae monótonamente** — oscila.
- El sustrato (el flujo) tiene la estructura; el instrumento (la grilla) es
  el que la distorsiona cerca de su límite de resolución.
- El experimento multi-observador (repo `Rylow999/rho-law`) muestra que hay
  cantidades invariants entre observadores (ratios entre ventanas, CV≈0.24)
  y cantidades que son artefacto (G_total, CV≈0.7).

Ver `NOUS/RHO_LAW/` para el marco unificador trans-dominio.

*Per Aspera, Ad Astra.*
