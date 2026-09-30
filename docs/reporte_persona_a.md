# Contribuciones al Informe PCAM - Jonathan Zacarias (Persona A)
**Documento base para ensamblar en `explicacion_pcam.pdf`**  
**Proyecto:** Simulación N-Body 2D con OpenMP  
**Rol:** Jonathan Zacarias (Base secuencial, corrección numérica y análisis de partición)  

---

## 1. Introducción y Baseline

### 1.1 Planteamiento del Problema
La simulación de $N$-cuerpos (*N-Body Problem*) es un problema fundamental en la física computacional y la astrofísica que consiste en calcular el movimiento dinámico y la evolución temporal de un conjunto de $N$ partículas que interactúan mutuamente a través de fuerzas gravitatorias.

En un espacio bidimensional, cada cuerpo $i \in \{0, 1, \dots, N-1\}$ posee una posición $\mathbf{r}_i = (x_i, y_i)$, una velocidad $\mathbf{v}_i = (v_{x,i}, v_{y,i})$ y una masa $m_i > 0$. De acuerdo con la Ley de Gravitación Universal de Newton, la fuerza neta que actúa sobre la partícula $i$ ejercida por el resto de partículas $j \neq i$ está dada por:

$$\mathbf{F}_i = \sum_{j=0, j \neq i}^{N-1} G \frac{m_i m_j}{( \|\mathbf{r}_j - \mathbf{r}_i\|^2 + \epsilon^2 )^{3/2}} (\mathbf{r}_j - \mathbf{r}_i)$$

donde:
- $G$ es la constante gravitatoria universal (normalizada a $1.0$ en la simulación).
- $\epsilon^2$ es un factor de suavizado (*softening parameter*, fijado en $10^{-9}$) incorporado para evitar singularidades y fuerzas numéricamente divergentes cuando la distancia euclidiana entre dos partículas tiende a cero ($\|\mathbf{r}_j - \mathbf{r}_i\| \to 0$).

A partir de la aceleración resultante $\mathbf{a}_i = \frac{\mathbf{F}_i}{m_i}$, el estado físico de cada cuerpo evoluciona a lo largo del tiempo discreto mediante un integrador numérico simple (Euler semi-implícito o explícito con paso de integración $\Delta t = 0.01$ s):

$$\mathbf{v}_i(t + \Delta t) = \mathbf{v}_i(t) + \mathbf{a}_i(t) \cdot \Delta t$$
$$\mathbf{r}_i(t + \Delta t) = \mathbf{r}_i(t) + \mathbf{v}_i(t + \Delta t) \cdot \Delta t$$

### 1.2 Complejidad Computacional y Motivación de la Paralelización
El cálculo directo (*all-pairs*) implica evaluar la interacción entre todos los pares ordenados de partículas. Para cada uno de los $N$ cuerpos, se deben computar las interacciones con los restantes $N-1$ cuerpos, lo que requiere un total de:

$$I_{\text{total}} = N(N - 1) = \mathcal{O}(N^2) \text{ evaluaciones de fuerza por paso de tiempo.}$$

Posteriormente, la fase de actualización cinemática de posiciones y velocidades demanda únicamente $\mathcal{O}(N)$ operaciones. Por consiguiente, para valores típicos de simulación ($N = 5000$ cuerpos, 10 pasos temporales), la fase de cómputo de fuerzas representa más del **99.5% del tiempo total de procesamiento**. Esta desproporción convierte al cálculo de fuerzas en el candidato óptimo para la descomposición y paralelización mediante la metodología **PCAM** y la API **OpenMP**.

### 1.3 Implementación de Referencia Secuencial (Baseline)
La implementación secuencial (`nbody_secuencial.c` / `src/nbody_seq.c`) sirve como línea base absoluta para cuantificar la aceleración (*speedup*) y validar la exactitud de los resultados.
- Se compila con optimización estándar: `gcc -O2 -fopenmp nbody_secuencial.c -o nbody_secuencial -lm`.
- Estructura desacoplada en dos fases:
  1. **Fase de Fuerzas:** Doble bucle anidado que acumula la fuerza en vectores auxiliares `fx[i]` y `fy[i]`.
  2. **Fase de Integración:** Bucle unidimensional que actualiza $v$ y $x, y$.
- Este desacoplamiento en dos fases es crucial, pues garantiza que ninguna partícula altere su posición durante el cómputo de fuerzas de un mismo paso temporal, eliminando dependencias de lectura y escritura concurrentes en la posterior paralelización.

---

## 2. Metodología PCAM: Fase de Partición (Partition)

La metodología PCAM (Foster, 1995) inicia con la fase de **Partition** (Particionamiento), cuyo propósito es exponer el máximo paralelismo posible dividiendo tanto los datos como las operaciones en tareas concurrentes elementales de grano fino.

```
            +---------------------------------------+
            |  Espacio del Problema: N Cuerpos      |
            +---------------------------------------+
                                |
             [Particionamiento Funcional y de Datos]
                                |
            +-------------------+-------------------+
            |                                       |
    [Dominio de Particulas]                 [Dominio de Interacciones]
   Grano medio: O(N) tareas                Grano fino: O(N^2) tareas
  Tarea i: Fuerza neta sobre i           Tarea (i,j): Fuerza mutua par
  (Base de nbody_paralelo.c)              (Base de nbody_optimized.c)
```

### 2.1 Descomposición de Datos y Dominio
En el problema N-Body se identifican dos tipos fundamentales de datos:
1. **Datos de Estado de los Cuerpos (Lectura global por paso):** Posiciones actuales $\mathbf{r}$ y masas $m$. En cada paso temporal son constantes durante la evaluación de fuerzas.
2. **Datos de Acumulación (Escritura):** Componentes de fuerza $\mathbf{F}_i$ y nuevos estados cinemáticos.

### 2.2 Estrategias de Particionamiento
Se evaluaron dos estrategias de partición:

1. **Partición por Partícula / Grano Medio ($O(N)$ tareas concurrentes):**
   - **Definición de Tarea:** La tarea elemental $T_i$ consiste en calcular la fuerza neta total ejercida sobre el cuerpo $i$ sumando la influencia de todos los cuerpos $j \in \{0, \dots, N-1\}$.
   - **Independencia:** Cada tarea $T_i$ escribe exclusivamente en su propio acumulador de fuerza `fx[i]` y `fy[i]`. Todas las lecturas sobre el arreglo `bodies` son compartidas y de solo lectura.
   - **Grado de Concurrencia:** Exactamente $N$ tareas independientes en paralelo.
   - **Condiciones de Carrera:** Cero. No existe riesgo de carreras de datos (*race-free by design*), siendo la base directa de `nbody_paralelo.c`.

2. **Partición por Pares de Interacción / Grano Fino ($O(N^2 / 2)$ tareas concurrentes):**
   - **Definición de Tarea:** La tarea elemental $T_{i,j}$ computa la fuerza mutua entre el par $(i, j)$ aprovechando la Tercera Ley de Newton ($\mathbf{F}_{ji} = -\mathbf{F}_{ij}$).
   - **Grado de Concurrencia:** $\frac{N(N-1)}{2}$ tareas elementales.
   - **Conflicto:** Múltiples tareas concurrentes intentan escribir simultáneamente en el cuerpo $i$ y en el cuerpo $j$, requiriendo mecanismos de sincronización o acumuladores privados (base de `nbody_optimized.c`).

---

## 3. Corrección Numérica y Validación de Resultados

### 3.1 Métrica de Validación: Checksum Acumulado
Para garantizar que las versiones paralela y optimizada reproduzcan fielmente la física de la versión secuencial, se diseñó la función de **Checksum Acumulado**:

$$C(t) = \sum_{i=0}^{N-1} \left( x_i(t) + y_i(t) + v_{x,i}(t) + v_{y,i}(t) \right)$$

### 3.2 Análisis de Asociatividad en Aritmética de Coma Flotante (IEEE-754)
En la computación secuencial, la suma de las fuerzas sobre un cuerpo $i$ sigue un orden estricto de índices:
$$\mathbf{F}_i = \sum_{j=0}^{N-1} \mathbf{f}_{ij} = (\dots((\mathbf{f}_{i,0} + \mathbf{f}_{i,1}) + \mathbf{f}_{i,2}) + \dots)$$

En el estándar IEEE-754 de coma flotante de doble precisión (`double`), la adición **no es matemáticamente asociativa** debido al truncamiento y redondeo de bits menos significativos en la mantisa cuando los sumandos difieren en órdenes de magnitud.

- En `nbody_paralelo.c` (Hugo Barillas), como cada hilo computa internamente el bucle $j$ en el orden idéntico $0 \dots N-1$, el resultado numérico es **idéntico bit a bit** al secuencial:
  $$\Delta C = |C_{\text{paralelo}} - C_{\text{secuencial}}| = 0.0$$
- En `nbody_optimized.c` (June Herrera), la acumulación de pares $(i, j)$ mediante reducción o reducción en variables privadas altera el orden cronológico de suma. Esto induce un error numérico residual propio del redondeo:
  $$\epsilon_{\text{rel}} = \frac{|C_{\text{optimizado}} - C_{\text{secuencial}}|}{|C_{\text{secuencial}}|} \le 10^{-11}$$
Cualquier desviación superior a este umbral ($> 10^{-7}$) denota inequívocamente una condición de carrera o un error en la formulación de los pares.

---

## 4. Respuestas a Preguntas de Análisis Asignadas (Preguntas 1 y 2)

### Pregunta 1: ¿Cómo influye el tamaño del problema ($N$) en la relación cómputo/comunicación y en el Speedup obtenido?
**Respuesta:**  
La simulación directa de N-Body presenta una complejidad de cálculo de fuerzas de orden $\mathcal{O}(N^2)$ y una sobrecarga de sincronización/comunicación de memoria compartida de orden $\mathcal{O}(N)$ (o $\mathcal{O}(p)$ con respecto a barreras y creación de hilos).

1. **Para valores pequeños de $N$ ($N \le 500$):**
   El volumen de operaciones aritméticas en punto flotante es reducido. El tiempo necesario para que el runtime de OpenMP cree y despache el equipo de hilos (*fork-join*), asigne iteraciones en la cola de planificación y sincronice la barrera implícita al final de cada paso temporal supera o iguala el tiempo de cálculo útil. En consecuencia, el *Speedup* es significativamente sub-lineal o incluso inferior a 1 ($S < 1$, desaceleración por sobrecarga de gestión).
2. **Para valores grandes de $N$ ($N \ge 3000$ y $N = 5000$):**
   La cantidad de operaciones crece cuadráticamente ($N^2 = 25 \times 10^6$ interacciones por paso). La relación entre cómputo útil y sobrecarga de gestión de hilos aumenta drásticamente en favor del cómputo. Siguiendo la formulación de la **Ley de Amdahl** y la **Ley de Gustafson-Barsis**, la fracción paralelizable $f$ se aproxima al 99.8%, amortizando completamente el *overhead* de OpenMP y permitiendo que el *Speedup* escale de forma casi ideal para 2, 4 y 8 hilos.

---

### Pregunta 2: ¿Qué impacto tiene la granularidad de la partición y el balance de carga en la fase de particionamiento del problema N-Body?
**Respuesta:**  
El particionamiento determina el tamaño de las tareas concurrentes y su regularidad computacional:

1. **Granularidad en `nbody_paralelo.c` (Partición sobre $i$):**
   - **Granularidad Uniforme:** Cada tarea $T_i$ evalúa exactamente $N-1$ interacciones. La carga de trabajo es homogénea y simétrica entre todas las tareas.
   - **Balanceo de Carga:** Al ser el costo computacional idéntico para cualquier $i$, la planificación estática (`schedule(static)`) distribuye bloques iguales de cuerpos ($N/p$) a cada hilo, garantizando un balance de carga óptimo sin incurrir en la sobrecarga que acarrearía la planificación dinámica (`dynamic`).
2. **Granularidad en `nbody_optimized.c` (Partición por pares simétricos $j > i$):**
   - **Granularidad No Uniforme:** Si se paraleliza el bucle externo $i$, el cuerpo $i=0$ debe computar $N-1$ interacciones con $j > 0$, mientras que el cuerpo $i=N-2$ solo computa 1 interacción con $j=N-1$. Esto genera una matriz triangular superior de iteraciones con carga altamente asimétrica.
   - **Impacto en el Mapeo:** Una planificación estática ingenua asignaría a los primeros hilos una carga de trabajo inmensamente superior a los últimos hilos (que terminarían prematuramente y quedarían ociosos en la barrera). En esta variante, una granularidad más fina combinada con `schedule(dynamic, chunk)` o `schedule(guided)` es indispensable para mitigar el desbalance computacional.
