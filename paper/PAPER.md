# SDDF — Curvatura espectral G[u]

**Extensión de LOGOS, dentro de [vega-vault](https://github.com/Rylow999/vega-vault)**

Autor: Luciano Benjamín Nieto (Delorien) · Auditoría técnica asistida: Claude (Anthropic)
Versión de este documento: v3 — 15 de septiembre de 2026

---

## Resumen

Se define un observable, la *curvatura espectral*

$$G[u] = \int \left(\frac{d\ln E(k)}{d\ln k}\right)^2 d\ln k,$$

sobre el espectro de energía $E(k)$ de un flujo turbulento, y se estudia su
comportamiento asintótico en el rango inercial de Kolmogorov (K41). Se
demuestra —numérica y analíticamente, con forma cerrada exacta— que

$$G^*[u] = \frac{25}{12}\ln Re + b(\delta) + \mathcal{O}(Re^{-1}),$$

con prefactor $25/12 = (25/9)\cdot(3/4)$, el cuadrado de la pendiente de
Kolmogorov por el exponente de escala de la microescala de disipación
($\eta \propto Re^{-3/4}$). El prefactor no depende de $\beta$, $\delta$,
$C$ ni $k_{\min}$. A partir de la forma cerrada se construye un estimador
del exponente espectral que corrige el sesgo de ventana de $Re$ finito y
que funciona con un único espectro (sección 7); se generaliza el resultado
a dos dimensiones (sección 8); y se contrasta con la propuesta de Migdal de
una firma log-periódica ligada a los ceros de $\zeta$ (sección 9),
concluyendo que $G[u]$ no es el observable adecuado para ese contraste.

Este documento consolida tres rondas de auditoría (v1, v2, v3) sobre el
mismo análisis. El historial completo de qué se corrigió en cada ronda está
en `AUDITORIA.md`; acá se presenta solo el resultado final.

---

## 1. Introducción

SDDF (curvatura espectral) es una extensión del programa LOGOS dentro de
`vega-vault`. El objetivo es doble: (a) caracterizar $G[u]$ como
diagnóstico del rango inercial —¿qué tan rápido y de qué forma crece con
$Re$?— y (b) evaluar si puede servir como observable para contrastar
predicciones recientes sobre estructura fina del espectro de turbulencia
decayente (Migdal, arXiv:2604.12207 y arXiv:2511.02165).

## 2. Definición y marco teórico

Se trabaja con espectros del tipo Pao/Pope generalizado,

$$E(k) = C\,\varepsilon^{2/3} k^{-q} \exp\!\big[-\beta (k\eta)^p\big],
\qquad \eta = (\nu^3/\varepsilon)^{1/4},$$

con $q = 5/3+\mu$ ($\mu$: corrección de intermitencia del exponente
espectral) y $p=4/3$ (forma de Pao). La pendiente logarítmica local es

$$s(x) = -q - \beta p\, x^p, \qquad x = k\eta.$$

El corte del rango inercial se define por desviación de pendiente:
$|s+q|=\delta$, lo que equivale a cortar en $x_c=(\delta/(\beta p))^{1/p}$.

## 3. Resultado central

Con ventana física $[k_{\min}, k_c]$, $k_{\min}$ fijo, $\varepsilon=1$,
$Re=1/\nu$:

$$G^*(Re,\delta) = \frac{3}{4}q^2\ln Re + b(\delta)
   - 2q\beta k_{\min}^p Re^{-3p/4} - \frac{\beta^2 p}{2}k_{\min}^{2p}Re^{-3p/2}.$$

Para K41 puro ($q=5/3$) el prefactor asintótico es
$a=(3/4)(5/3)^2=25/12=2.08333$. Verificado numéricamente contra el
integrador con error 0.003–0.13 % (`datos/01_resultados_G.csv`,
`datos/02_validacion_analitica.csv`).

## 4. Forma cerrada exacta

$$G^*(Re,\delta) = \frac{25}{12}\ln\!\left(\frac{3\delta\,Re}{4\beta}\right)
   + \frac{5}{2}\delta + \frac{3}{8}\delta^2
   - \frac{10\beta/3}{Re} - \frac{2\beta^2/3}{Re^2}$$

(caso $q=5/3$, $p=4/3$, $k_{\min}=1$; la forma general con $q,p$
arbitrarios está en `codigo/sddf_exact.py::b_delta`). Esto confirma que la
"corrección de Re finito" del orden del 2.6 % detectada en v1 es física —
no numérica— y tiene ahora fórmula cerrada (`datos/08_forma_cerrada_terminos.csv`).

## 5. El corte por cuantización de grilla

Cortar en el primer nodo de grilla que viola $\delta$ introduce un error
que **no decrece monótonamente con N** (de ahí las oscilaciones de la tabla
de convergencia de v1). Interpolando linealmente el punto de corte en
$(\ln k, s)$ en vez de tomar el nodo más cercano
(`truncate_by_slope_interp`), el error baja tres órdenes de magnitud con
los mismos 2000 puntos:

| N puntos | corte en grilla | corte interpolado |
|---|---|---|
| 100 | 1.461 % | 0.127 % |
| 1000 | 0.180 % | 0.001 % |
| 2000 | 0.245 % | 0.0003 % |
| 100000 | 0.003 % | 0.00002 % |

(`datos/07_convergencia_grilla.csv`)

## 6. La quimera Pao/Pope: por qué $\beta$ importa

El barrido original usa $E(k)\propto k^{-5/3}\exp[-\beta(k\eta)^{4/3}]$ con
$\beta=5.2$. Esa forma funcional (exponente $4/3$) es de **Pao (1965)**,
donde $\beta$ no es libre: $\beta=(3/2)C_K$, que con $C_K=1.5$ da
$\beta=2.25$, no 5.2. El valor 5.2 es la constante de **Pope (2000)**, pero
corresponde a otra función de corte —aproximadamente exponencial en
$k\eta$, no potencia $4/3$. Usar $\beta=5.2$ con exponente $4/3$ corta la
cola de disipación 1.87× antes de lo que corresponde, y desplaza $b(\delta)$
en $+1.745$ nats. El prefactor $25/12$ no se mueve —verificado sobre tres
modelos de corte (quimera, Pao consistente, Pope completo) en
`datos/11_modelos_espectrales.csv`— pero cualquier calibración absoluta de
$G^*$ sí.

**Hallazgo adicional de esta ronda (v3):** la elección de $\beta$ no solo
desplaza $b(\delta)$: si se usa para *medir* el exponente espectral con el
estimador de la sección 7, un $\beta$ mal especificado en un factor ~2
(5.2 en vez de 2.25) borra hasta el 85 % de una señal de intermitencia real
inyectada de $\mu=0.033$ (ver sección 7.2). Es un segundo canal de sesgo,
del mismo orden que el de $Re$ finito, y hay que controlarlo con el mismo
cuidado antes de reportar un $\mu$ medido en datos reales.

## 7. El estimador principal: `q_corregido`

### 7.1 Definición

En vez de dividir $G^*/\ln(k_c/k_{\min})$ —que converge al exponente
asintótico solo como $1/\ln Re$, con sesgo $\mu_{\rm eff}=0.084$ a
$Re=100,\delta=0.10$, cuatro veces la intermitencia física—, se despeja $q$
de la forma cerrada de la sección 4 de forma iterativa:

$$G = q^2\ln(x_c/x_1) + 2q\beta(x_c^p-x_1^p) + \frac{\beta^2 p}{2}(x_c^{2p}-x_1^{2p}),$$

iterando $q\leftarrow\sqrt{(G-t_2(q)-t_3)/L}$ desde $q_0=5/3$ (3–5
iteraciones alcanzan 6 cifras, porque $q$ entra linealmente solo en el
término $t_2$). Implementado en `codigo/estimador_principal.py`, función
`medir_exponente_espectral(k, E, eta, delta, beta)`.

### 7.2 Validación

Sobre los 5 espectros K41 originales, en ambos $\delta$ probados:
**$q=1.666667\pm0.000001$** en todos los $Re$ — sin necesitar el barrido
completo, con un solo espectro (`tests/test_estimador_principal.py`,
`test_k41_puro`).

Sobre espectros con intermitencia inyectada ($\mu=0.033$, modelo Pao
consistente, $\beta=2.25$): el estimador recupera **$\mu=0.033000$ exacto**
en los 5 $Re$, siempre que se le pase el $\beta$ correcto
(`test_intermitencia_con_beta_correcto`). Si se le pasa el $\beta$ de Pope
(5.2) en cambio de el de Pao (2.25) —error de especificación de modelo, no
de medición—, el $\mu$ recuperado cae a 0.0056: se pierde el 83 % de la
señal (`test_beta_mal_especificado_degrada_la_medicion`). Esta es la
verificación numérica del hallazgo de la sección 6.

**Recomendación operativa:** `q_corregido` es el estimador a usar de acá en
más para medir el exponente espectral efectivo en datos reales (DNS,
experimento), reemplazando cualquier ajuste log-lineal simple sobre un
barrido en $Re$. Requiere conocer $\beta$ del modelo de corte con la
precisión suficiente — no asumir el valor de la literatura sin verificar
qué forma funcional de corte tiene el espectro real.

## 8. Generalización y estado del paper 2D

Para $E(k)\sim k^{-q}$ con corte físico en $k_c\propto Re^\gamma$:

$$a \equiv \frac{dG^*}{d\ln Re} = q^2\gamma.$$

| régimen | $q$ | $\gamma$ | $a$ predicho |
|---|---|---|---|
| 3D, K41 (este documento) | 5/3 | 3/4 | **25/12 = 2.083** |
| 3D con intermitencia real ($\zeta_2\approx0.70$) | 1.700 | 3/4 | 2.1675 |
| 2D, cascada de enstrofía ($E\sim k^{-3}$) | 3 | 1/2 | **9/2 = 4.5** |
| 2D, cascada inversa (acotada por la caja) | 5/3 | 0 | **0** |

**Alerta sobre el paper 2D existente:** $G^*\approx3634$ con $q=3$ exigiría
$\ln(k_c/k_{\min})=3634/9\approx404$ nats, es decir $k_c/k_{\min}\sim
e^{404}$ — matemáticamente imposible en cualquier simulación real. Es
estructuralmente el mismo artefacto de "corte sin ventana física definida"
que en este documento daba $G=182017$ espurio a $Re=100$ antes de aplicar
el corte por pendiente (sección 5).

**Estado al cierre de esta ronda (v3): no se pudo recalcular.** El
recálculo requiere los espectros DNS crudos del paper 2D (no solo el
número final $G^*\approx3634$, $a\approx0.70$), que no forman parte de
este paquete. La receta a aplicar cuando estén disponibles es directa:
correr `truncate_by_slope_interp` + `medir_exponente_espectral` de la
sección 7 sobre esos espectros, con el $\beta$ correcto del modelo de corte
que se haya usado ahí (ver sección 6). Si el exponente resultante sigue
cerca de 0.70 tras esa corrección, hay una conexión real entre la
intermitencia 3D y el escalamiento 2D que vale la pena reportar; si cambia
como pasó acá con 0.175→log, era el mismo artefacto de umbral.

## 9. Conexión con Migdal: el observable correcto

arXiv:2604.12207 (*Decaying Turbulence and the Riemann Hypothesis*) y
arXiv:2511.02165 (*Geometric Solution of Turbulence as Diffusion in Loop
Space*), ambos verificados como preprints reales, predicen oscilaciones
log-periódicas en el espectro residual ligadas a los ceros no triviales de
$\zeta$. Se inyectó una modulación $E(k)\to E(k)[1+a\cos(\omega\ln k)]$ y se
midieron dos observables: el ajuste de $G[u]$, y el periodograma de la
pendiente local $s(\ln k)$.

| amplitud $a$ | ventana (nats de $\ln k$) | $\Delta q$ vía $G[u]$ | $\omega$ detectada | SNR |
|---|---|---|---|---|
| 0.020 | 11.95 (diseñada para el test) | $2\times10^{-5}$ | **2.064** (inyectado: 2.0) | 77 |
| 0.050 | 11.95 | $1\times10^{-3}$ | 2.006 | 292 |
| 0.020 | 3.72 (ventana real del barrido, $Re=2000$) | — | 3.356 (falla) | 11 |

**Conclusión:** $G[u]$ es casi ciego a la modulación log-periódica —la
integral la corrige solo a segundo orden, $(a\omega)^2/2$, indistinguible
del ruido de $Re$ finito. El observable correcto es el periodograma (o la
transformada de Fourier) de la pendiente local en $\ln k$, y hace falta una
ventana inercial de al menos ~10 nats para que el pico sea significativo
—con 3.7 nats (el barrido actual) el resultado es ruido puro, y con datos
DNS de $128^3$ es imposible por no tener rango inercial suficiente; con
$4096^3$ probablemente sí. El plan correcto para la Fase 2 del roadmap
("conexión con Migdal") es ir directo del espectro residual (tras restar
K41) al periodograma o a la transformada de Mellin, sin pasar por $G[u]$.

## 10. Contexto verificado: Buckmaster/Alpöge/OpenAI (8 de septiembre 2026)

Contrastado contra fuentes primarias: Tristan Buckmaster (NYU) y Levent
Alpöge (Anthropic) publicaron tres preprints con blowup en tiempo finito y
forzamiento suave (medio poroso incompresible, Boussinesq 2D, Euler 3D
incompresible), formalizados en Lean. OpenAI anunció por separado un
resultado interno de blowup para Navier-Stokes 3D forzado, sin prueba
publicada a la fecha de esa cobertura; reconoció la prioridad de
Buckmaster/Alpöge en el caso de Euler forzado y aclaró que no reclama el
premio Clay por este resultado. **Punto relevante para SDDF: todos estos
resultados son sobre versiones forzadas de las ecuaciones** — la página
oficial del problema Clay sigue listando Navier-Stokes sin forzamiento como
no resuelto, así que el enfoque de este programa en flujo libre sigue
siendo un ángulo distinto, no superado por ninguno de estos resultados.

## 11. Limitaciones abiertas

- **`spectrum_grueso.csv` — recuperado.** El archivo original (21 puntos,
  $G=33.392638$) apareció después del cierre de v3 y se verificó bit a
  bit: $G[u]$ sin truncar sobre esos 21 puntos da `33.39263798313726`.
  Resultó ser, más que un dato de más, el caso construido para mostrar el
  problema de la sección 5: con esta grilla tan gruesa, el corte por
  pendiente falla (cae por debajo de 3 puntos) para todo $\delta\le0.5$, y
  solo produce un resultado significativo con $\delta=1.0$
  ($G^*=10.188$). Ver `AUDITORIA.md`, hallazgo A7, para el detalle
  completo — incluye la verificación de por qué era, en sí mismo, el
  ejemplo del bug de cuantización.
- **El recálculo del paper 2D con el corte corregido queda pendiente** de
  los espectros DNS crudos de ese paper (sección 8).
- **La conexión con Migdal** solo se probó sobre espectros sintéticos; la
  aplicación a DNS real de $4096^3$ (o mayor) queda para cuando haya datos
  disponibles con suficiente rango inercial resuelto.

## 12. Reproducibilidad

```bash
cd codigo
python3 barrido_re_v3.py ../datos/espectros
python3 sddf_completo.py ../datos/espectros ../datos ../figuras
./rust/bin/sddf ../datos/espectros/spectrum_Re_1000.csv 0.5
python3 estimador_principal.py
cd .. && python3 tests/test_estimador_principal.py
```

## Referencias

- Pao, Y.-H. (1965). *Structure of turbulent velocity and scalar fields at
  large wavenumbers.*
- Pope, S. B. (2000). *Turbulent Flows.* Cambridge University Press.
- Migdal, A. (2026). *Decaying Turbulence and the Riemann Hypothesis.*
  arXiv:2604.12207.
- Migdal, A. (2025). *Geometric Solution of Turbulence as Diffusion in Loop
  Space.* arXiv:2511.02165.
- Buckmaster, T., Alpöge, L. (2026). Declaración pública del 8 de
  septiembre de 2026 y preprints asociados (medio poroso incompresible,
  Boussinesq 2D, Euler 3D incompresible).
- OpenAI (2026). *Navier-Stokes solution* (anuncio interno, blowup forzado).
  openai.com/index/navier-stokes-solution.

---

*Historial completo de auditoría (v1→v2→v3): ver `AUDITORIA.md`.*

*Per Aspera, Ad Astra.*
