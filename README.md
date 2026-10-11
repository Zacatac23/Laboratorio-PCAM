# Laboratorio PCAM + OpenMP: Simulación N-Body 2D

Proyecto de computación paralela aplicando la metodología **PCAM** (Partitioning, Communication, Agglomeration, Mapping) y la biblioteca **OpenMP** en lenguaje C para la simulación bidimensional de interacción gravitatoria entre $N$-cuerpos.

---

## 1. Estructura del Repositorio

El proyecto maneja una estructura interna modular para desarrollo colaborativo y una correspondencia directa con los archivos requeridos para la entrega en Canvas:

```
Laboratorio-PCAM/
├── bin/                          # Binarios compilados (generados)
├── docs/                         # Documentación técnica y borradores de informe
│   ├── reporte_persona_a.md      # Secciones del informe asignadas a Persona A (P, Intro, Preguntas 1 y 2)
│   └── tabla_correccion_numerica.md # Definición del Checksum y tabla de validación
├── src/                          # Código fuente de desarrollo interno
│   ├── common.h                  # Estructura Body, generador LCG (seed=42), checksum y utilidades
│   └── nbody_seq.c               # Implementación secuencial de referencia (Persona A)
├── nbody_secuencial.c            # Archivo secuencial autónomo para entrega en Canvas (Persona A)
├── Makefile                      # Automatización de compilación y pruebas rápidas
├── plan_lab_pcam_nbody_v3.md     # Plan de trabajo consensuado por el equipo
└── README.md                     # Documentación principal del proyecto
```

### Correspondencia de Archivos para Canvas

| Archivo de Trabajo Interno | Archivo Entregable en Canvas | Responsable |
|---|---|---|
| `src/nbody_seq.c` / `nbody_secuencial.c` | `nbody_secuencial.c` | **Persona A** |
| `src/nbody_parallel.c` | `nbody_paralelo.c` | **Persona B** |
| `src/nbody_optimized.c` | `nbody_optimized.c` (o en `.zip`) | **Persona C** |
| `docs/reporte_persona_a.md` + aportes B y C | `explicacion_pcam.pdf` | **Equipo** (Ensamblado por A) |
| Capturas + tablas experimentales | `evidencia.pdf` | **Equipo** (Ensamblado por A) |

---

## 2. Requisitos de Compilación y Entorno

- **Compilador:** `gcc` con soporte para OpenMP y C99 o superior.
- **Banderas de compilación estándar:** `-O2 -fopenmp -Wall` y enlace matemático `-lm`.
- **Plataformas compatibles:** Linux, Windows (MinGW / MSYS2 / WSL), macOS.

---

## 3. Instrucciones de Compilación

### Opción A: Usando el Makefile
```bash
# Compilar todos los ejecutables
make all

# Compilar solo el secuencial
make bin/nbody_seq
```

### Opción B: Compilación manual directa con GCC

1. **Compilar la versión de desarrollo interno:**
   ```bash
   mkdir -p bin
   gcc -O2 -fopenmp src/nbody_seq.c -o bin/nbody_seq -lm
   ```

2. **Compilar el archivo directo de Canvas:**
   ```bash
   gcc -O2 -fopenmp nbody_secuencial.c -o nbody_secuencial -lm
   ```

3. **Compilar la versión paralela (Persona B):**
   ```bash
   gcc -O2 -fopenmp src/nbody_parallel.c -o bin/nbody_parallel -lm   # desarrollo
   gcc -O2 -fopenmp nbody_paralelo.c -o nbody_paralelo -lm           # entregable Canvas
   ```

---

## 4. Instrucciones de Ejecución y Parámetros

Los ejecutables aceptan hasta 5 argumentos posicionales por línea de comandos:

```bash
./bin/nbody_seq [N] [steps] [threads] [schedule] [chunk]
```

### Parámetros:
1. `N` (entero): Número de cuerpos a simular. **Por defecto: `5000`**.
2. `steps` (entero): Número de pasos temporales de integración. **Por defecto: `10`**.
3. `threads` (entero): Número de hilos solicitados (fijado a 1 en secuencial). **Por defecto: `1`**.
4. `schedule` (cadena): Política de planificación (`sequential`, `static`, `dynamic`, `guided`). **Por defecto: `"sequential"`**.
5. `chunk` (entero): Tamaño de bloque para OpenMP. **Por defecto: `0`** (distribución por defecto).

### Ejemplos de Ejecución:

- **Ejecución con parámetros por defecto oficiales ($N=5000$, 10 pasos):**
  ```bash
  ./bin/nbody_seq
  ```

- **Ejecución rápida de prueba ($N=1000$, 10 pasos):**
  ```bash
  ./bin/nbody_seq 1000 10 1 sequential 0
  ```

- **Ejecución directa del entregable Canvas:**
  ```bash
  ./nbody_secuencial 5000 10 1 sequential 0
  ```

- **Versión paralela** (`schedule` = `static` | `dynamic` | `guided`; `chunk` = 0 usa el valor por defecto de OpenMP):
  ```bash
  ./bin/nbody_parallel 5000 10 8 dynamic 8
  ./nbody_paralelo 5000 10 8 static 0
  ```

### Experimentos de rendimiento (Persona B)

```bash
bash scripts/run_experiments.sh          # N=5000, 10 pasos, 5 repeticiones -> results.csv
python scripts/analyze_results.py        # -> results_summary.csv, docs/tabla_experimentos.{md,tex}
```

- `run_experiments.sh` compila desde cero, hace una corrida de calentamiento que se descarta y ejecuta la secuencial y las configuraciones paralelas (1/2/4/8 hilos con `static`; con 8 hilos `static` chunk 8/64, `dynamic` 8/64 y `guided`). Las repeticiones se **intercalan** (ronda 1 de todas las configuraciones, luego ronda 2, …) para que la deriva térmica afecte a todas por igual. Si existe `src/nbody_optimized.c`, también lo mide.
- `results.csv` guarda **cada corrida** (columna `rep`). `analyze_results.py` calcula el promedio, la desviación estándar, $S(p)=T_1/T_p$, $E(p)=S(p)/p$, la fracción serial de Karp-Flatt y el error del checksum frente al secuencial.
- En laptop: conectar el cargador y cerrar otros programas antes de medir.

---

## 5. Salida del Programa y Formato de Registro

El programa produce un resumen formateado en pantalla y una línea CSV para alimentar de forma directa los scripts de análisis y generación de `results.csv`:

```text
=================================================================
         Simulacion N-Body 2D - Laboratorio PCAM + OpenMP        
=================================================================
  Version              : nbody_secuencial (A)
  Cuerpos (N)          : 5000
  Pasos (Steps)        : 10
  Semilla (Seed)       : 42
  Threads solicitados  : 1
  Threads usados       : 1
  Schedule             : sequential
  Chunk size           : 0
  Procesadores disp.   : 8
  Max hilos OpenMP     : 8
  Tiempo de ejecucion  : 1.482000 segundos
  Checksum final       : -4.0266058332e+004
=================================================================
CSV_HEADER: version,N,steps,seed,threads_req,threads_used,schedule,chunk,procs,max_threads,time_sec,checksum
CSV_ROW: nbody_secuencial (A),5000,10,42,1,1,sequential,0,8,8,1.482000,-4.0266058332e+004
```

---

## 6. Validación de Corrección Numérica (Checksum)

Antes de realizar las corridas de rendimiento de las versiones paralela y optimizada, consulte el documento:
- [`docs/tabla_correccion_numerica.md`](file:///c:/Users/jonat/Desktop/Laboratorio-PCAM/docs/tabla_correccion_numerica.md)

### Criterio de Aceptación:
- **`nbody_paralelo.c`:** El checksum debe coincidir exactamente ($\Delta = 0$) con el secuencial.
- **`nbody_optimized.c`:** El checksum debe presentar un error relativo menor a $10^{-11}$ atribuible a la no asociatividad del punto flotante en la reducción simétrica.

---

## 7. Reparto de Responsabilidades

| Integrante | Rol | Entregables Principales |
|---|---|---|
| **Jonathan Zacarias** (Persona A) | Base, corrección y ensamblado | `nbody_seq.c`, `nbody_secuencial.c`, `common.h`, checksum de referencia, tabla de corrección numérica, `README.md`, ensamblado final de `evidencia.pdf` y `explicacion_pcam.pdf`. |
| **Hugo Barillas** (Persona B) | Paralelo y experimentos | `nbody_parallel.c`, script de experimentos, `results.csv`, tablas de speedup y eficiencia, respuestas a preguntas 3, 4, 7 y 8. |
| **June Herrera** (Persona C) | Optimización y concurrencia | `nbody_optimized.c`, análisis de race conditions, comparación parallel vs. optimized, respuestas a preguntas 5 y 6. |
