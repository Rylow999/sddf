# Contrastar la firma de Migdal en DNS real — guía completa

**Objetivo:** determinar si el espectro de energía de una turbulencia
decayente real muestra la modulación log-periódica que predicen
Migdal (arXiv:2604.12207, *Decaying Turbulence and the Riemann
Hypothesis*) y (arXiv:2511.02165, *Geometric Solution of Turbulence as
Diffusion in Loop Space*): oscilaciones ligadas a los ceros no triviales
$1/2+i\rho_n$ de $\zeta(s)$, con polos de Mellin en $p=-8+i\rho_n$ y
activaciones temporales $t_n\propto\rho_n^3$.

Esta guía no recorta por factibilidad: enumera todo lo que hace falta,
de punta a punta, y marca en cada paso qué tan realista es conseguirlo
hoy. Lo que ya está resuelto en el repo `SDDF` (`paper/PAPER.md`,
sección 9) es el diagnóstico de *qué* observable usar; esto es el plan
para aplicarlo a datos reales.

---

## 1. Qué predicción, exactamente, se está contrastando

Migdal predice que el **espectro residual** (lo que queda después de
restar el ajuste K41, $E(k)-E_{\rm K41}(k)$) o, de forma más directa, la
**pendiente logarítmica local** $s(\ln k) = d\ln E/d\ln k$, contiene una
superposición de modos con frecuencias angulares $\omega_n\propto\rho_n$
(las partes imaginarias de los ceros de $\zeta$: $\rho_1\approx14.13$,
$\rho_2\approx21.02$, $\rho_3\approx25.01$, ...). El observable a
construir es el **periodograma** (o la transformada de Mellin/Fourier)
de $s(\ln k)$, buscando picos en esas frecuencias — no picos genéricos,
sino picos en las posiciones *predichas de antemano* por los ceros de
$\zeta$, lo que hace al test falsable.

Sección 9 del paper de SDDF ya estableció, con espectros sintéticos, que
$G[u]$ (la curvatura espectral) es el observable equivocado para esto:
la integral suaviza la modulación a segundo orden y la vuelve
indistinguible del ruido de $Re$ finito. Este documento asume que se usa
el periodograma de $s(\ln k)$, no $G[u]$.

---

## 2. Requisito central: rango inercial resuelto suficientemente ancho

### 2.1 Cuánto rango hace falta

El test con espectros sintéticos (paper, sección 9) mostró que con solo
3.72 nats de $\ln k$ (la ventana real del barrido usado ahí, $Re=2000$)
el resultado es ruido puro (SNR=11, frecuencia detectada muy alejada de
la inyectada). Con 11.95 nats diseñados a propósito, el pico aparece
limpio (SNR=77–292). La frecuencia más lenta relevante es
$\omega_1\approx2\rho_1\approx28.3$ (factor 2 porque la relación exacta
entre $\rho_n$ y la frecuencia en $\ln k$ depende de la normalización de
Mellin que use cada paper — verificar contra la ecuación exacta antes de
fijar el número). Regla práctica conservadora: **se necesitan al menos
3–4 períodos completos de la oscilación más lenta esperada**, es decir
del orden de **10–15 nats de $\ln k$ resueltos en el rango inercial**
(no en el espectro completo — el rango inercial es la porción entre la
escala de forzado/integral y la escala de disipación).

### 2.2 Qué resolución de grilla implica esto

$10$–$15$ nats $=e^{10}$–$e^{15}\approx 22000$–$3.3\times10^6$ en razón
$k_{\max}/k_{\min}$ dentro del rango inercial. Como el rango inercial es
solo una fracción del rango total resuelto por la simulación (el resto
se lo llevan la escala de forzado grande y la disipación pequeña), la
resolución total de la grilla DNS ($N^3$) tiene que ser bastante mayor
que eso. La relación estándar es

$$k_{\max}\eta \sim 1 \;\Rightarrow\; N \sim k_{\max}/k_{\min}
   \sim Re^{3/4} \quad (\text{ley de Kolmogorov de grados de libertad}).$$

Para $Re^{3/4}\gtrsim 10^4$–$10^5$ (que es lo que un rango inercial de
10–15 nats normalmente exige, ya que el forzado ocupa las primeras
1–2 décadas), hace falta $Re\gtrsim10^{5.3}$–$10^{6.7}$, lo cual en
términos de $Re_\lambda$ (Reynolds de Taylor, la cifra que reportan las
bases de datos DNS) ronda **$Re_\lambda \gtrsim 500$–$1000$**, y en
términos de grilla, **$N\gtrsim 4096$**, idealmente $8192^3$.

### 2.3 Qué hay disponible hoy (verificado, septiembre 2026)

La **Johns Hopkins Turbulence Database (JHTDB)**, `turbulence.pha.jhu.edu`,
sigue operativa (con ventanas de mantenimiento programadas, la más
reciente anunciada para agosto 2026) y ofrece:

- Un dataset de turbulencia isotrópica forzada en grilla $1024^3$
  ($Re_\lambda\approx433$), con 5028 snapshots temporales — el dataset
  "de trabajo" estándar, accesible vía Web Services (Python, Matlab,
  Fortran, C) o descarga de subcubos en HDF5.
- Snapshots individuales (sin resolución temporal completa) de DNS de
  $4096^3$ y $8192^3$ a Reynolds más alto — estos son los que hacen
  falta para tener el rango inercial ancho que pide la sección 2.1; al
  ser snapshots aislados (no una serie temporal completa), sirven para
  el análisis espectral estático pero **no** para verificar directamente
  la predicción temporal $t_n\propto\rho_n^3$ (sección 5 de esta guía).
- El propio paper de Migdal (arXiv:2604.12207) reporta haber contrastado
  su predicción contra **DNS de $4096^3$** — el orden de magnitud de
  resolución que hace falta coincide con lo estimado en 2.2, y es
  evidencia de que el contraste es factible con datos ya existentes en
  algún lugar (no necesariamente públicos todavía).

**Conclusión de esta sección:** el dataset público más simple de JHTDB
($1024^3$, $Re_\lambda\approx433$) **no alcanza** — es el mismo problema
que ya se identificó para grillas de $128^3$ en el paper de SDDF, un
escalón más arriba pero con la misma limitación estructural. Hacen
falta los snapshots de $4096^3$/$8192^3$, o generar una simulación
propia a esa resolución (sección 4).

---

## 3. Pipeline de análisis paso a paso

Asumiendo que se consigue acceso a un campo de velocidad 3D (de JHTDB o
de una simulación propia) a la resolución necesaria:

1. **Calcular el espectro de energía isotrópico $E(k)$** a partir del
   campo de velocidad 3D: transformada de Fourier de $\mathbf{u}$,
   promediar $|\hat{\mathbf u}(\mathbf k)|^2$ en cáscaras esféricas de
   $|\mathbf k|$. (Herramientas: `numpy.fft`/`pyfftw`/`mpi4py-fft` para
   campos grandes; JHTDB entrega directamente utilidades para esto en su
   paquete `pyJHTDB`, o se puede pedir el cutout del campo y calcular
   localmente.)
2. **Ajustar y restar K41**: ajustar $E(k)\sim C k^{-5/3}$ (o, mejor, con
   el $q$ que dé `medir_exponente_espectral` de SDDF sobre el mismo
   espectro) en el rango inercial, y calcular el residuo
   $R(k) = \ln E(k) - \ln E_{\rm ajuste}(k)$, o directamente trabajar con
   la pendiente local $s(\ln k)$ (recomendado, porque no requiere fijar
   de antemano dónde "empieza" el ajuste).
3. **Periodograma / transformada de Mellin de $s(\ln k)$** sobre la
   ventana del rango inercial identificada con `truncate_by_slope_interp`
   de SDDF (sección 5 del paper) en ambos extremos —no solo en el de
   disipación, también en el de forzado/integral, que es donde el rango
   inercial "empieza" en un espectro real. Usar una ventana de Hann o
   similar antes de la FFT para reducir fuga espectral (*leakage*), dado
   que el rango es finito y los bordes no son periódicos.
4. **Buscar picos en las frecuencias predichas** $\omega_n\propto\rho_n$
   ($\rho_1\approx14.13, \rho_2\approx21.02, \rho_3\approx25.01,\dots$),
   verificando primero contra la ecuación exacta del paper de Migdal qué
   relación de escala convierte $\rho_n$ en una frecuencia de $\ln k$ (no
   asumir que es 1:1; el paper de SDDF usó $\omega=2.0$ para el $n=1$
   como ejemplo de test, no como la predicción real de Migdal).
5. **Significancia estadística**: no alcanza con "hay un pico ahí". Hace
   falta:
   - una estimación del ruido de fondo (nivel de piso del periodograma
     lejos de las frecuencias predichas, o un test de permutación
     aleatorizando la fase);
   - corrección por comparaciones múltiples si se testean varios $\rho_n$
     a la vez;
   - repetir el análisis sobre **más de un snapshot/realización**
     (si hay varias disponibles) para separar señal de fluctuación
     estadística de una sola muestra — un solo campo DNS es una sola
     realización del proceso turbulento, no un promedio de ensamble.

---

## 4. Si no hay acceso a DNS de resolución suficiente: generarla

Si JHTDB no ofrece un snapshot adecuado (o si se necesita la serie
temporal completa para la parte de la sección 5, que JHTDB no da a esa
resolución), la alternativa es correr una simulación pseudo-espectral
propia:

- **Método:** pseudo-espectral con des-aliasing (regla de los 2/3),
  integración temporal Runge-Kutta de 2º-4º orden, forzado de banda
  angosta en números de onda bajos para mantener energía constante
  (turbulencia forzada) o sin forzado tras una condición inicial (para
  contrastar específicamente el régimen de decaimiento libre, que es el
  que trata arXiv:2604.12207).
- **Costo computacional:** una simulación $4096^3$ necesita del orden de
  $4096^3\times3\text{ (componentes)}\times8\text{ bytes}\approx1.6$ TB
  de memoria solo para el campo de velocidad en doble precisión, más
  espacio de trabajo para FFTs (típicamente 3–5× eso) → **del orden de
  5–10 TB de RAM distribuida**, en un clúster HPC con red de baja
  latencia (InfiniBand) para las transposiciones globales de la FFT 3D.
  Para $8192^3$, multiplicar por 8. Esto está fuera del alcance de
  hardware personal o de un único nodo; requiere tiempo de cómputo en un
  centro HPC (allocation vía XSEDE/ACCESS en EE.UU., o el equivalente
  local en Argentina — SNCAD/CCAD, aunque probablemente con cupos muy por
  debajo de lo que esto exige).
- **Software:** códigos pseudo-espectrales open-source existen
  (`Dedalus`, `hit3d`, o extender el propio motor en Rust del repo SDDF a
  3D con FFT distribuida vía `mpi4py-fft`/`PFFT`) — pero escribir y
  validar un solver DNS 3D robusto es un proyecto en sí mismo, no un
  script de análisis.

**Evaluación honesta:** generar la propia DNS a la resolución necesaria
es, para un investigador independiente sin acceso a un clúster HPC
grande, la opción menos realista de las dos. Conseguir acceso a
snapshots de mayor resolución de JHTDB (o de otro repositorio público
equivalente, o directamente contactando al grupo de Migdal para pedir
los datos que usaron en su propio contraste) es el camino más corto.

---

## 5. La parte temporal: $t_n\propto\rho_n^3$

Esta predicción es más exigente que la espectral: requiere una
**serie temporal** de espectros durante el decaimiento libre —no un
snapshot único— para ver si las "activaciones" (posiblemente,
transiciones o eventos característicos en la evolución del espectro)
ocurren en los tiempos predichos. JHTDB tiene series temporales
completas solo en el dataset de $1024^3$ (forzado, no decayente), que
por la sección 2.3 no alcanza en resolución para el análisis espacial. Un
contraste completo de la parte temporal probablemente requiere una
simulación propia de decaimiento libre (sin forzado, condición inicial
turbulenta) guardando snapshots frecuentes — otro proyecto de cómputo
propio, sección 4.

---

## 6. Resumen de lo que hace falta, en orden de factibilidad

| Requisito | Disponibilidad hoy | Vía más corta |
|---|---|---|
| Pipeline de análisis (espectro → residuo → periodograma) | **Ya se puede escribir** | Extender `SDDF/codigo` con FFT 3D y periodograma; no depende de tener los datos todavía |
| Snapshot DNS $\gtrsim4096^3$, alto $Re_\lambda$ | Existe en JHTDB (snapshots de $4096^3$/$8192^3$) pero sin resolución temporal completa | Pedir acceso/descarga vía JHTDB Web Services o cutout service |
| Significancia estadística con múltiples realizaciones | Requiere varios snapshots independientes | Ver si JHTDB tiene más de un snapshot a esa resolución; si no, limitarse a un test de permutación sobre el único campo disponible y reportarlo como limitación |
| Serie temporal de decaimiento libre a resolución suficiente (para $t_n\propto\rho_n^3$) | **No disponible públicamente** que se sepa hoy | Simulación propia en HPC, o contactar directamente al grupo de Migdal por los datos de su propio contraste |
| Cómputo propio a $4096^3$–$8192^3$ si nada de lo anterior alcanza | Requiere ~5–10+ TB RAM distribuida, cluster HPC con InfiniBand, allocation de cómputo | SNCAD/CCAD (Argentina) u otro centro con cupos de ese tamaño; expectativa realista: meses de trámite, no un experimento de fin de semana |

**La recomendación concreta:** empezar por el pipeline de análisis
(sección 3) sobre el snapshot $1024^3$ de JHTDB, sabiendo de antemano
que va a dar un resultado no concluyente por falta de rango inercial —
sirve para validar el código de extremo a extremo sobre datos reales
(no solo sintéticos) antes de necesitar los datos de mayor resolución.
En paralelo, gestionar el acceso a los snapshots de $4096^3$/$8192^3$ de
JHTDB, que es la pieza que probablemente destrabe un contraste real.

---

*Este documento es la contraparte de `paper/PAPER.md`, sección 9, para
cuando haya datos reales disponibles. No reemplaza ese análisis —lo
extiende con lo que hace falta para pasar de "el observable correcto es
el periodograma de $s(\ln k)$" a "medido en DNS real, el resultado es
X".*
