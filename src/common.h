/**
 * ============================================================================
 * common.h - Definiciones y utilidades compartidas para N-Body 2D
 * Laboratorio PCAM + OpenMP
 * 
 * Contiene:
 *  - Estructura Body (x, y, vx, vy, mass)
 *  - Constantes fisicas y valores por defecto
 *  - Generador determinista de condiciones iniciales (seed = 42)
 *  - Calculo de Checksum acumulado
 *  - Formato estandarizado de impresion de resumen y CSV
 * ============================================================================
 */

#ifndef COMMON_H
#define COMMON_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

/* Constantes fisicas de la simulacion */
#ifndef G_CONST
#define G_CONST 1.0            /* Constante gravitacional normalizada */
#endif

#ifndef SOFTENING
#define SOFTENING 1e-9         /* Factor de suavizado para evitar singularidades (r -> 0) */
#endif

#ifndef DT
#define DT 0.01                /* Tamano del paso temporal de integracion */
#endif

/* Valores por defecto acordados por el equipo */
#define DEFAULT_N 5000
#define DEFAULT_STEPS 10
#define DEFAULT_SEED 42

/**
 * Estructura que representa una particula/cuerpo en el espacio bidimensional.
 */
typedef struct {
    double x;      /* Posicion en el eje X */
    double y;      /* Posicion en el eje Y */
    double vx;     /* Velocidad en el eje X */
    double vy;     /* Velocidad en el eje Y */
    double mass;   /* Masa del cuerpo */
} Body;

/**
 * Generador pseudoaleatorio determinista LCG (Linear Congruential Generator).
 * Garantiza exactamente las mismas condiciones iniciales y checksums de
 * referencia en cualquier sistema operativo (Linux, Windows, macOS) y compilador.
 */
static inline double lcg_rand(unsigned long long *state) {
    *state = (*state * 6364136223846793005ULL + 1442695040888963407ULL);
    return (double)(*state >> 11) / (double)(1ULL << 53);
}

/**
 * Inicializa las condiciones de los N cuerpos de forma determinista usando la semilla provista.
 *
 * @param bodies Arreglo de estructuras Body con capacidad para n elementos.
 * @param n      Numero de cuerpos a inicializar.
 * @param seed   Semilla para el generador pseudoaleatorio (por defecto 42).
 */
static inline void init_bodies(Body *bodies, int n, unsigned int seed) {
    unsigned long long state = (unsigned long long)seed;
    if (state == 0) state = 42;

    for (int i = 0; i < n; i++) {
        /* Posiciones en el plano 2D acotadas en [-1000.0, 1000.0] */
        bodies[i].x = (lcg_rand(&state) * 2000.0) - 1000.0;
        bodies[i].y = (lcg_rand(&state) * 2000.0) - 1000.0;

        /* Velocidades iniciales moderadas en [-10.0, 10.0] */
        bodies[i].vx = (lcg_rand(&state) * 20.0) - 10.0;
        bodies[i].vy = (lcg_rand(&state) * 20.0) - 10.0;

        /* Masa positiva en [1.0, 100.0] */
        bodies[i].mass = (lcg_rand(&state) * 99.0) + 1.0;
    }
}

/**
 * Calcula el checksum de la simulacion como la suma acumulada de las componentes
 * de posicion y velocidad de todos los cuerpos: sum(x + y + vx + vy).
 *
 * Permite validar la correccion numerica entre las versiones secuencial,
 * paralela y optimizada.
 *
 * @param bodies Arreglo de cuerpos.
 * @param n      Numero total de cuerpos.
 * @return       Valor de coma flotante de doble precision con el checksum.
 */
static inline double compute_checksum(const Body *bodies, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += (bodies[i].x + bodies[i].y + bodies[i].vx + bodies[i].vy);
    }
    return sum;
}

/**
 * Imprime el resumen formateado de ejecucion requerido por los acuerdos del laboratorio:
 * N, pasos, seed, threads solicitados y usados, schedule, chunk, tiempo (omp_get_wtime()),
 * checksum, omp_get_num_procs() y omp_get_max_threads().
 * 
 * Tambien emite una linea CSV con prefijo 'CSV_ROW:' para facilitar la automatizacion
 * de la toma de datos por parte de la Persona B.
 */
static inline void print_summary(
    const char *version_name,
    int n,
    int steps,
    unsigned int seed,
    int threads_requested,
    int threads_used,
    const char *schedule,
    int chunk,
    double elapsed_time,
    double checksum
) {
    printf("=================================================================\n");
    printf("         Simulacion N-Body 2D - Laboratorio PCAM + OpenMP        \n");
    printf("=================================================================\n");
    printf("  Version              : %s\n", version_name);
    printf("  Cuerpos (N)          : %d\n", n);
    printf("  Pasos (Steps)        : %d\n", steps);
    printf("  Semilla (Seed)       : %u\n", seed);
    printf("  Threads solicitados  : %d\n", threads_requested);
    printf("  Threads usados       : %d\n", threads_used);
    printf("  Schedule             : %s\n", schedule ? schedule : "sequential");
    printf("  Chunk size           : %d\n", chunk);
    printf("  Procesadores disp.   : %d\n", omp_get_num_procs());
    printf("  Max hilos OpenMP     : %d\n", omp_get_max_threads());
    printf("  Tiempo de ejecucion  : %.6f segundos\n", elapsed_time);
    printf("  Checksum final       : %.10e\n", checksum);
    printf("=================================================================\n");

    /* Formato CSV para el script de recoleccion de resultados */
    printf("CSV_HEADER: version,N,steps,seed,threads_req,threads_used,schedule,chunk,procs,max_threads,time_sec,checksum\n");
    printf("CSV_ROW: %s,%d,%d,%u,%d,%d,%s,%d,%d,%d,%.6f,%.10e\n",
           version_name, n, steps, seed, threads_requested, threads_used,
           schedule ? schedule : "sequential", chunk, omp_get_num_procs(),
           omp_get_max_threads(), elapsed_time, checksum);
}

#endif /* COMMON_H */
