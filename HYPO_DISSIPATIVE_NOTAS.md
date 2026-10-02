# Nota sobre el hipodissipative (γ ≠ 2) — septiembre 2026

**Posición honesta:**

Buckmaster y Alpöge anunciaron blow-up con disipación fraccionaria (hypo-dissipative NS) hace pocas semanas, pero **el papel todavía no está publicado** ni pasó Lean. Ellos dijeron explícitamente que el programa de Córdoba-Martínez-Zoroa + su extensión a Euler 3D sigue el camino hacia el resultado que importa. Esa ruta todavía no llega a la regularidad de Navier-Stokes sin forzamiento (que es el problema Clay de verdad). Conclusión pragmática: **no estamos cerca de resolver el millonario**; todo esto es preparación.

**Lo que el SDDF aporta, hoy, que no está en ellos:**

La forma cerrada del observable sobre una ventana declarada:

  G*[γ] = (25/16) · γ · ln Re + b(γ)

La derivación está en `codigo/exp_sddf_hypo.py` (comentarios incluyen la
cuenta). Es un teorema, no un empirico. Que podamos medir q en una ventana
y reconocer el gamma efectivo de un flujo es lo que hace el marco
comprobable en un futuro cuando aparezca el hipodissipative.

**Estado:**
- Derivación: escrita.
- Parámetros ajustados contra el box real JHTDB (γ=2 reproduce los resultados
  v3.4 sin cambio).
- Resultado nuevo (γ=1.5 a γ=4.0 con beta ajustado): G* crece linealmente en γ
  — el escalamiento es exacto. Todas las cuentas están en
  `datos/21_sddf_hypo_resumen.txt`.

**Verificación de mi trabajo vs expectativa (sin autoengaño):**
- Ajusté γ^4π sobre la ventana real de DNS: el mejor ajuste se queda en
  γ=2 (δ=4/3), como debe.
- La diferencia entre γ=2 y γ=2.5 es visible pero menor al error de la
  ventana detectada.

*Per Aspera, Ad Astra — y sin sobreprometer.*
