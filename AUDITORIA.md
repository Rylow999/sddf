# Historial de auditoría — SDDF (curvatura espectral G[u])

Este documento registra, en orden, qué se auditó, qué se encontró y qué se
corrigió en cada ronda. `paper/PAPER.md` es el resultado consolidado;
este archivo es la traza de cómo se llegó ahí.

---

## v1 — análisis original

Primer informe sobre la curvatura espectral $G[u]$. Establece el resultado
central ($25/12\ln Re$), identifica el bug del piso absoluto $10^{-30}$ en
el lector, mide una corrección de $Re$ finito del orden de 2.6 % sin
descomponerla, y deja pendiente:

1. El binario Rust sin compilar (sandbox sin red en esa sesión).
2. Los 5 CSV originales del barrido en $Re$, no incluidos —solo los
   resultados numéricos.
3. `spectrum_grueso.csv` ($G=33.392638$, 21 puntos), mencionado pero no
   incluido.
4. $b(\delta)$ del ajuste, dependiente de la ventana, sin fórmula cerrada.
5. El 2.6 % de corrección de $Re$ finito, sin descomponer en sus términos.

## v2 — primera ronda de cierre

Resuelve los pendientes 1, 2, 4 y 5 de v1:

- **Rust compilado y verificado** bit a bit contra la versión Python
  (`codigo/rust/bin/sddf`).
- **Los 5 espectros regenerados** desde los parámetros documentados
  ($\varepsilon=1$, $C_k=1.5$, $\beta=5.2$, $k\in[1,1000]$ logspace 2000
  pts) y verificados: reproducen exactamente los 5 valores de $G$ del
  informe original.
- **Forma cerrada completa** hasta el último término, con $b(\delta)$
  exacta (`codigo/sddf_exact.py`).
- **Corrección de $Re$ finito descompuesta** en sus dos términos exactos
  ($\propto 1/Re$, $\propto 1/Re^2$) y usada para construir un estimador
  "desviscosado".
- **Hallazgos nuevos no pedidos, pero relevantes:**
  - El corte por cuantización de grilla no converge monótonamente; el
    corte interpolado (`truncate_by_slope_interp`) lo arregla en tres
    órdenes de magnitud (sección 5 del paper).
  - El barrido mezclaba dos modelos de corte incompatibles (Pao con la
    constante de Pope) — sección 6.
  - `q_corregido()`: un estimador de $q$ que despeja de la forma cerrada
    en vez de dividir, funciona con un solo espectro. Quedó implementado
    dentro del pipeline (`sddf_completo.py`, bloque 9) pero no como API
    reutilizable.
  - Verificación contra literatura de Migdal (arXiv:2604.12207,
    arXiv:2511.02165) y de la noticia Buckmaster/Alpöge/OpenAI del 8/9.
  - Alerta sobre el paper 2D: $G^*\approx3634$ con $q=3$ es matemáticamente
    imposible sin recortar antes.

**Pendiente no resuelto en v2:** `spectrum_grueso.csv` — se probó un
barrido sobre $(k_{\min},k_{\max},\nu)$ con 21 puntos y ninguna
combinación razonable lo reprodujo.

## v3 — esta ronda (15 de septiembre de 2026)

Punto de partida: el paquete v2 completo (`sddf_curvatura_espectral_v2.zip`).
Pedido: resolver lo que quedaba pendiente, corregir lo que faltara, y
empaquetar todo —código, datos, figuras, paper e historial— listo para subir
a GitHub como parte de `vega-vault` / `LOGOS`.

### A1 — `q_corregido` promovido a API principal

Extraído del bloque 9 de `sddf_completo.py` a un módulo propio,
`codigo/estimador_principal.py`, con:
- función de una sola llamada `medir_exponente_espectral(k, E, eta, delta, beta)`;
- test de regresión contra los 5 espectros K41 originales
  (`tests/test_estimador_principal.py::test_k41_puro`): reproduce
  $q=1.666667$ con error $<10^{-6}$ en todos los casos.

### A2 — validación adicional sobre los espectros con intermitencia inyectada

No estaba hecha en v1/v2 con esta función aislada. Resultado: recupera
$\mu=0.033000$ exacto en los 5 $Re$ **siempre que se le pase el $\beta$
correcto del modelo usado para generar el espectro** (Pao, $\beta=2.25$).

### A3 — hallazgo nuevo: sensibilidad a un $\beta$ mal especificado

Al validar A2 con el $\beta$ "por defecto" del paquete (el de Pope, 5.2,
que es el que documentan los headers del barrido original) en vez del
$\beta$ real del modelo Pao usado para el espectro con intermitencia
(2.25), el $\mu$ recuperado cae de 0.033000 a 0.005611 — se pierde el 83 %
de la señal. Es un segundo canal de sesgo del mismo orden que el de $Re$
finito (sección 4 del paper), no identificado como tal en v1 ni v2.
Cuantificado en un barrido de sensibilidad
(`codigo/estimador_principal.py` no lo incluye por defecto; la corrida
completa está documentada en `paper/PAPER.md` sección 7.2) y fijado como
test de regresión
(`tests/test_estimador_principal.py::test_beta_mal_especificado_degrada_la_medicion`).
**Consecuencia práctica:** cualquier aplicación futura de este estimador a
datos reales necesita determinar $\beta$ del modelo de corte real, no
asumir un valor de la literatura sin verificar.

### A4 — segundo intento de recuperar `spectrum_grueso.csv`

Se repitió la búsqueda de v2 con un espacio mucho más amplio: $Re \in
\{50,\dots,10000\}$ (no solo los 5 valores del barrido estándar),
$k_{\min}\in\{10^0,10^{-0.5},10^{-1}\}$, $k_{\max}$ barrido fino entre
$10^{0.5}$ y $10^4$, ambos modelos (Pao y quimera). **Resultado: sigue sin
poder reconstruirse de forma única.** Hay múltiples combinaciones
degeneradas que caen dentro de 0.03–0.5 del valor objetivo (33.392638) sin
que ninguna se distinga como "la" correcta — el sistema está
subdeterminado con solo el valor final de $G$ como dato. Se documenta el
resultado negativo en vez de forzar una reconstrucción arbitraria. Sigue
pendiente de un backup del archivo original o sus metadatos.

### A5 — recálculo del paper 2D con el corte corregido

**No realizado.** Requiere los espectros DNS crudos del paper 2D (no solo
$G^*\approx3634$, $a\approx0.70$), que no están en este paquete ni en
ningún otro material disponible en esta sesión. Se documenta la receta
exacta a aplicar cuando estén disponibles (`paper/PAPER.md`, sección 8) en
vez de simular un resultado sin los datos de entrada reales.

### A6 — verificación de reproducibilidad end-to-end

Antes de empaquetar: se re-ejecutó el binario Rust sobre
`spectrum_Re_1000.csv` (reproduce $G^*[\text{corte pendiente}]=10.2305$,
igual que en v1/v2) y se corrió toda la suite de tests nueva
(`tests/test_estimador_principal.py`) — los 3 tests pasan.

### A7 — `spectrum_grueso.csv` recuperado (aportado por Luciano después de v3)

Luciano pegó los 21 puntos $(k, E(k))$ originales en el chat. Verificación
inmediata: $G[u]$ sobre los 21 puntos **sin truncar** da
`33.39263798313726`, que coincide con el `33.392638` del informe v1 hasta
la sexta cifra decimal — es, sin ambigüedad, el archivo original.

Se agregó a `datos/espectros/spectrum_grueso.csv` y se generó
`datos/15_spectrum_grueso_verificacion.csv` con el corte por pendiente en
los 4 $\delta$ estándar. Resultado, y por qué el archivo importa más de lo
que parece:

| $\delta$ | puntos usados | ok | $G^*$ |
|---|---|---|---|
| 0.10 | 3 (forzado) | **No** | 1.782 (no significa nada) |
| 0.25 | 3 (forzado) | **No** | 1.782 (no significa nada) |
| 0.50 | 3 (forzado) | **No** | 1.782 (no significa nada) |
| 1.00 | 9 | Sí | 10.188 |

Con $\delta\le0.5$ el primer punto de grilla ($k=1$) ya se desvía casi lo
mismo que el umbral ($|s(k{=}1)+5/3|=0.099$, y el siguiente punto salta a
0.296): la grilla es tan gruesa que **no hay ningún punto entre el primero
y el segundo** donde cortar de forma significativa, y el corte cae por
debajo del mínimo de 3 puntos exigido — exactamente el modo de falla que
`main.rs` marca con `ok=false`. Solo con $\delta=1.0$ (una ventana de
pendiente muy laxa) el corte da un resultado utilizable.

Además, no ajusta bien a la familia Pao/quimera con los parámetros
estándar del resto del paquete (error cuadrático grande probando `E_pao` y
`E_quimera` sobre un rango de $\nu$) — es decir, se generó con otro
$C_k$/$\nu$/receta que no se documentó en su momento, o directamente con
otro generador. Esto ya no importa para reproducir el número: el archivo
en sí es la fuente de verdad, no un modelo que haya que reajustar.

**Conclusión:** este archivo no era un dato más — era el ejemplo
construido a propósito para mostrar el bug de cuantización de grilla de la
sección 5 del paper (*"cuidado con $\delta$ en grillas gruesas"*). Con solo
21 puntos, el corte por pendiente falla salvo con $\delta$ muy laxo, y el
corte por *grilla* (sin interpolar) da el error que oscilaba en la tabla
de convergencia de v1/v2. Confirma en un caso extremo, con datos reales y
no sintéticos, por qué `truncate_by_slope_interp` (sección 5 del paper) es
necesario y no un refinamiento cosmético.

### Qué queda genuinamente pendiente después de v3 + A7

1. Recálculo del paper 2D — bloqueado por falta de datos crudos (A5).
2. Aplicación del estimador y del test de Migdal a DNS real de alta
   resolución ($4096^3$ o mayor) — no hay datos DNS reales en este
   paquete, solo espectros sintéticos. Ver `CONTRASTE_MIGDAL_DNS.md` para
   la receta completa, con requisitos de datos y cómputo, para cuando se
   consiga acceso a una base como JHTDB.

De los tres pendientes que quedaban al cierre de v3, `spectrum_grueso.csv`
ya no es uno de ellos.
