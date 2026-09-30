/**
 * ============================================================================
 * nbody_secuencial.c - Implementacion Secuencial (Baseline) para Canvas
 * Laboratorio PCAM + OpenMP: N-Body 2D
 * Responsable: Persona A
 *
 * Compilacion:
 *   gcc -O2 -fopenmp nbody_secuencial.c -o nbody_secuencial -lm
 *
 * Ejecucion:
 *   ./nbody_secuencial [N] [steps] [threads] [schedule] [chunk]
 *
 * Parametros por defecto:
 *   N = 5000, steps = 10, threads = 1, schedule = "sequential", chunk = 0
 *   Semilla fija = 42
 * ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

/* Constantes fisicas acordadas */
#define G_CONST 1.0            /* Constante gravitacional normalizada */
#define SOFTENING 1e-9         /* Factor de suavizado para evitar singularidades */
#define DT 0.01                /* Paso temporal de integracion */

/* Valores por defecto */
#define DEFAULT_N 5000
#define DEFAULT_STEPS 10
#define DEFAULT_SEED 42

/**
 * Estructura que representa un cuerpo en 2D.
 */
typedef struct {
    double x;      /* Posicion X */
    double y;      /* Posicion Y */
    double vx;     /* Velocidad X */
    double vy;     /* Velocidad Y */
    double mass;   /* Masa del cuerpo */
} Body;

/**
 * Generador pseudoaleatorio determinista LCG (Linear Congruential Generator).
 * Produce valores identicos en cualquier plataforma con la misma semilla.
 */
static inline double lcg_rand(unsigned long long *state) {
    *state = (*state * 6364136223846793005ULL + 1442695040888963407ULL);
    return (double)(*state >> 11) / (double)(1ULL << 53);
}

/**
 * Inicializa los cuerpos de manera determinista con seed = 42.
 */
static inline void init_bodies(Body *bodies, int n, unsigned int seed) {
    unsigned long long state = (unsigned long long)seed;
    if (state == 0) state = 42;

    for (int i = 0; i < n; i++) {
        bodies[i].x = (lcg_rand(&state) * 2000.0) - 1000.0;
        bodies[i].y = (lcg_rand(&state) * 2000.0) - 1000.0;
        bodies[i].vx = (lcg_rand(&state) * 20.0) - 10.0;
        bodies[i].vy = (lcg_rand(&state) * 20.0) - 10.0;
        bodies[i].mass = (lcg_rand(&state) * 99.0) + 1.0;
    }
}

/**
 * Calcula el checksum de la simulacion como sum(x + y + vx + vy).
 */
static inline double compute_checksum(const Body *bodies, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += (bodies[i].x + bodies[i].y + bodies[i].vx + bodies[i].vy);
    }
    return sum;
}

/**
 * Imprime el resumen formateado y linea CSV estandarizada.
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

    printf("CSV_HEADER: version,N,steps,seed,threads_req,threads_used,schedule,chunk,procs,max_threads,time_sec,checksum\n");
    printf("CSV_ROW: %s,%d,%d,%u,%d,%d,%s,%d,%d,%d,%.6f,%.10e\n",
           version_name, n, steps, seed, threads_requested, threads_used,
           schedule ? schedule : "sequential", chunk, omp_get_num_procs(),
           omp_get_max_threads(), elapsed_time, checksum);
}

int main(int argc, char *argv[]) {
    int n = DEFAULT_N;
    int steps = DEFAULT_STEPS;
    int threads_req = 1;
    const char *schedule = "sequential";
    int chunk = 0;
    unsigned int seed = DEFAULT_SEED;

    if (argc > 1) n = atoi(argv[1]);
    if (argc > 2) steps = atoi(argv[2]);
    if (argc > 3) threads_req = atoi(argv[3]);
    if (argc > 4) schedule = argv[4];
    if (argc > 5) chunk = atoi(argv[5]);

    int threads_used = 1;

    Body *bodies = (Body *)malloc(sizeof(Body) * n);
    double *fx = (double *)malloc(sizeof(double) * n);
    double *fy = (double *)malloc(sizeof(double) * n);

    if (!bodies || !fx || !fy) {
        fprintf(stderr, "Error al reservar memoria para los cuerpos.\n");
        free(bodies);
        free(fx);
        free(fy);
        return EXIT_FAILURE;
    }

    init_bodies(bodies, n, seed);

    double start_time = omp_get_wtime();

    /* Bucle temporal de simulacion */
    for (int step = 0; step < steps; step++) {
        /* Fase 1: Calculo de fuerzas O(N^2) */
        for (int i = 0; i < n; i++) {
            double fxi = 0.0;
            double fyi = 0.0;

            for (int j = 0; j < n; j++) {
                if (i == j) continue;

                double dx = bodies[j].x - bodies[i].x;
                double dy = bodies[j].y - bodies[i].y;
                double dist_sqr = dx * dx + dy * dy + SOFTENING;
                double dist = sqrt(dist_sqr);

                double inv_dist3 = 1.0 / (dist * dist_sqr);
                double force_mag = G_CONST * bodies[i].mass * bodies[j].mass * inv_dist3;

                fxi += force_mag * dx;
                fyi += force_mag * dy;
            }

            fx[i] = fxi;
            fy[i] = fyi;
        }

        /* Fase 2: Actualizacion de posiciones y velocidades O(N) */
        for (int i = 0; i < n; i++) {
            double ax = fx[i] / bodies[i].mass;
            double ay = fy[i] / bodies[i].mass;

            bodies[i].vx += ax * DT;
            bodies[i].vy += ay * DT;

            bodies[i].x += bodies[i].vx * DT;
            bodies[i].y += bodies[i].vy * DT;
        }
    }

    double end_time = omp_get_wtime();
    double elapsed_time = end_time - start_time;

    double final_checksum = compute_checksum(bodies, n);

    print_summary(
        "nbody_secuencial",
        n,
        steps,
        seed,
        threads_req,
        threads_used,
        schedule,
        chunk,
        elapsed_time,
        final_checksum
    );

    free(bodies);
    free(fx);
    free(fy);

    return EXIT_SUCCESS;
}
