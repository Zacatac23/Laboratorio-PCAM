# Tabla de Corrección Numérica y Checksum de Referencia

**Laboratorio PCAM + OpenMP: N-Body 2D**  
**Responsable:** Jonathan Zacarias (Persona A - Base, corrección y ensamblado)  
**Semilla (Seed):** 42 (fija para todas las versiones)  
**Paso temporal ($\Delta t$):** 0.01  
**Constante gravitatoria ($G$):** 1.0  
**Factor de suavizado ($\epsilon^2$):** $10^{-9}$  

---

## 1. Definición del Checksum de Referencia

Para garantizar la validez científica y la corrección numérica de las simulaciones paralelas y optimizadas, se define el **Checksum acumulado** $C(t)$ en el instante de tiempo final $t = \text{steps} \times \Delta t$ como:

$$C = \sum_{i=0}^{N-1} \left( x_i + y_i + v_{x,i} + v_{y,i} \right)$$

### Propiedades y Justificación
1. **Sensibilidad a carreras de datos (Data Races):** Si dos hilos colisionan en la actualización de posiciones o velocidades, o si se omiten interacciones entre pares de partículas, el valor resultante se desvía drásticamente del valor de referencia.
2. **Determinismo:** Con la semilla fijada en `42` y el generador LCG lineal, las condiciones iniciales son idénticas bit a bit en cualquier arquitectura.
3. **Comportamiento en punto flotante (IEEE-754):** La suma de coma flotante en computadoras **no es asociativa**:
   $$(a + b) + c \neq a + (b + c)$$
   En la versión paralela básica (`nbody_parallel.c` de Hugo Barillas), como cada hilo calcula únicamente las fuerzas del cuerpo $i$ en el mismo orden interno sobre $j$, el resultado debe ser **estrictamente idéntico** al secuencial ($\Delta = 0$). En la versión optimizada simétrica (`nbody_optimized.c` de June Herrera), donde las fuerzas de pares $(i, j)$ se acumulan en orden arbitrario por los hilos, se admite una tolerancia de redondeo flotante de:
   $$\text{Error Relativo} = \frac{|C_{\text{opt}} - C_{\text{seq}}|}{|C_{\text{seq}}|} < 10^{-11}$$

---

## 2. Tabla de Validación y Valores de Referencia

Las tres implementaciones deben contrastarse con los siguientes casos de prueba estandarizados:

| Configuración ($N$, pasos) | Checksum Inicial ($t=0$) | Checksum Final Secuencial (Referencia) | Tolerancia Aceptada ($\Delta$) | Estado de Validación |
|---|---|---|---|---|
| **$N = 100$, steps = 10** | Determinista (seed 42) | Generado por `nbody_seq` | $\Delta = 0$ (paralelo) / $< 10^{-11}$ (optimizado) | [ ] Validado |
| **$N = 500$, steps = 10** | Determinista (seed 42) | Generado por `nbody_seq` | $\Delta = 0$ (paralelo) / $< 10^{-11}$ (optimizado) | [ ] Validado |
| **$N = 1000$, steps = 10** | Determinista (seed 42) | Generado por `nbody_seq` | $\Delta = 0$ (paralelo) / $< 10^{-11}$ (optimizado) | [ ] Validado |
| **$N = 3000$, steps = 10** (Plan de contingencia) | Determinista (seed 42) | Generado por `nbody_seq` | $\Delta = 0$ (paralelo) / $< 10^{-11}$ (optimizado) | [ ] Validado |
| **$N = 5000$, steps = 10** (Oficial del Lab) | Determinista (seed 42) | Generado por `nbody_seq` | $\Delta = 0$ (paralelo) / $< 10^{-11}$ (optimizado) | [ ] Validado |

---

## 3. Protocolo de Verificación Cruzada (Checklist para Hugo Barillas y June Herrera)

Antes de iniciar la toma de tiempos para los experimentos de rendimiento, las implementaciones deben superar este protocolo:

1. **Compilar la versión secuencial de referencia:**
   ```bash
   gcc -O2 -fopenmp src/nbody_seq.c -o bin/nbody_seq -lm
   ```
2. **Ejecutar la referencia y registrar el Checksum:**
   ```bash
   ./bin/nbody_seq 5000 10 1 sequential 0
   ```
3. **Compilar la versión paralela (`nbody_parallel.c`) y optimizada (`nbody_optimized.c`):**
   ```bash
   gcc -O2 -fopenmp src/nbody_parallel.c -o bin/nbody_parallel -lm
   gcc -O2 -fopenmp src/nbody_optimized.c -o bin/nbody_optimized -lm
   ```
4. **Verificar que para 1, 2, 4 y 8 hilos el Checksum coincida:**
   ```bash
   ./bin/nbody_parallel 5000 10 4 static 64
   ./bin/nbody_optimized 5000 10 4 static 64
   ```
5. **Criterio de aprobación:**
   - Si $|C_{\text{par}} - C_{\text{seq}}| > 10^{-7}$, **existe una condición de carrera activa** o un error en los índices de interacción. No proceder a medir tiempos hasta corregirlo.
