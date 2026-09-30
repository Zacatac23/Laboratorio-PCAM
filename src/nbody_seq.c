/**
 * ============================================================================
 * nbody_seq.c - Implementacion Secuencial (Baseline) para N-Body 2D
 * Laboratorio PCAM + OpenMP
 * Responsable: Persona A
 *
 * Compilacion:
 *   gcc -O2 -fopenmp src/nbody_seq.c -o bin/nbody_seq -lm
 *
 * Uso:
 *   ./bin/nbody_seq [N] [steps] [threads] [schedule] [chunk]
 *
 * Parametros por defecto:
 *   N = 5000, steps = 10, threads = 1, schedule = "sequential", chunk = 0
 *   Semilla = 42
 * ============================================================================
 */

#include "common.h"

int main(int argc, char *argv[]) {
    /* 1. Lectura de argumentos con valores por defecto acordados */
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

    /* Forzar 1 hilo para la version secuencial */
    int threads_used = 1;

    /* 2. Reserva de memoria dinamica */
    Body *bodies = (Body *)malloc(sizeof(Body) * n);
    double *fx = (double *)malloc(sizeof(double) * n);
    double *fy = (double *)malloc(sizeof(double) * n);

    if (!bodies || !fx || !fy) {
        fprintf(stderr, "Error: No se pudo asignar memoria para %d cuerpos.\n", n);
        free(bodies);
        free(fx);
        free(fy);
        return EXIT_FAILURE;
    }

    /* 3. Inicializacion determinista de cuerpos con seed = 42 */
    init_bodies(bodies, n, seed);

    /* Checksum inicial de verificacion */
    double initial_checksum = compute_checksum(bodies, n);
    printf("Inicializacion completada. Cuerpos: %d, Pasos: %d, Semilla: %u\n", n, steps, seed);
    printf("Checksum inicial (t=0): %.10e\n", initial_checksum);

    /* 4. Medicion de tiempo de ejecucion con omp_get_wtime() */
    double start_time = omp_get_wtime();

    /* 5. Bucle principal de simulacion en el tiempo */
    for (int step = 0; step < steps; step++) {
        /* FASE 1: Calculo de fuerzas gravitacionales O(N^2)
         * Cada cuerpo i acumula la atraccion producida por todos los demas cuerpos j. */
        for (int i = 0; i < n; i++) {
            double fxi = 0.0;
            double fyi = 0.0;

            for (int j = 0; j < n; j++) {
                if (i == j) continue;

                double dx = bodies[j].x - bodies[i].x;
                double dy = bodies[j].y - bodies[i].y;
                double dist_sqr = dx * dx + dy * dy + SOFTENING;
                double dist = sqrt(dist_sqr);

                /* Fuerza gravitatoria: F = G * m_i * m_j / (dist^2)
                 * En forma vectorial normalizada: F_vec = (F / dist) * delta_vec
                 * Simplificado: G * m_i * m_j / (dist^3) * delta_vec */
                double inv_dist3 = 1.0 / (dist * dist_sqr);
                double force_mag = G_CONST * bodies[i].mass * bodies[j].mass * inv_dist3;

                fxi += force_mag * dx;
                fyi += force_mag * dy;
            }

            fx[i] = fxi;
            fy[i] = fyi;
        }

        /* FASE 2: Integracion de velocidades y posiciones O(N)
         * Se ejecuta en un bucle separado una vez calculadas todas las fuerzas del paso. */
        for (int i = 0; i < n; i++) {
            /* a = F / m */
            double ax = fx[i] / bodies[i].mass;
            double ay = fy[i] / bodies[i].mass;

            /* Integrador Euler simple */
            bodies[i].vx += ax * DT;
            bodies[i].vy += ay * DT;

            bodies[i].x += bodies[i].vx * DT;
            bodies[i].y += bodies[i].vy * DT;
        }
    }

    double end_time = omp_get_wtime();
    double elapsed_time = end_time - start_time;

    /* 6. Calculo del Checksum final */
    double final_checksum = compute_checksum(bodies, n);

    /* 7. Resumen de resultados */
    print_summary(
        "nbody_secuencial (A)",
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

    /* 8. Liberacion de recursos */
    free(bodies);
    free(fx);
    free(fy);

    return EXIT_SUCCESS;
}
