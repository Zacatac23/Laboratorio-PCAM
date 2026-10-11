/**
 * ============================================================================
 * nbody_parallel.c - Implementacion Paralela con OpenMP para N-Body 2D
 * Laboratorio PCAM + OpenMP
 * Responsable: Persona B
 *
 * Compilacion:
 *   gcc -O2 -fopenmp src/nbody_parallel.c -o bin/nbody_parallel -lm
 *
 * Uso:
 *   ./bin/nbody_parallel [N] [steps] [threads] [schedule] [chunk]
 *
 *   schedule : static | dynamic | guided   (por defecto static)
 *   chunk    : tamano de bloque; 0 = valor por defecto de OpenMP
 *
 * Parametros por defecto:
 *   N = 5000, steps = 10, threads = omp_get_max_threads(), schedule = static,
 *   chunk = 0, semilla = 42
 *
 * Estrategia (PCAM):
 *   - Partition: una tarea por cuerpo i (fuerza neta sobre i).
 *   - Communication: lectura compartida de bodies[]; cada tarea escribe solo
 *     fx[i] y fy[i], por lo que no hay escrituras concurrentes.
 *   - Agglomeration: las N tareas se agrupan en bloques de 'chunk' iteraciones.
 *   - Mapping: los bloques se asignan a los hilos segun 'schedule'
 *     (static/dynamic/guided), elegido en tiempo de ejecucion.
 * ============================================================================
 */

#include "common.h"

/**
 * Traduce el nombre del schedule recibido por linea de comandos al tipo de
 * OpenMP. Devuelve 0 si el nombre no es valido.
 */
static int parse_schedule(const char *name, omp_sched_t *kind) {
    if (strcmp(name, "static") == 0)  { *kind = omp_sched_static;  return 1; }
    if (strcmp(name, "dynamic") == 0) { *kind = omp_sched_dynamic; return 1; }
    if (strcmp(name, "guided") == 0)  { *kind = omp_sched_guided;  return 1; }
    return 0;
}

int main(int argc, char *argv[]) {
    /* 1. Lectura de argumentos con valores por defecto acordados */
    int n = DEFAULT_N;
    int steps = DEFAULT_STEPS;
    int threads_req = omp_get_max_threads();
    const char *schedule = "static";
    int chunk = 0;
    unsigned int seed = DEFAULT_SEED;

    if (argc > 1) n = atoi(argv[1]);
    if (argc > 2) steps = atoi(argv[2]);
    if (argc > 3) threads_req = atoi(argv[3]);
    if (argc > 4) schedule = argv[4];
    if (argc > 5) chunk = atoi(argv[5]);

    omp_sched_t sched_kind;
    if (n <= 0 || steps < 0 || threads_req <= 0 || chunk < 0 ||
        !parse_schedule(schedule, &sched_kind)) {
        fprintf(stderr, "Uso: %s [N] [steps] [threads] [static|dynamic|guided] [chunk]\n", argv[0]);
        return EXIT_FAILURE;
    }

    /* Numero de hilos y schedule del ciclo de fuerzas. Con chunk = 0 OpenMP
     * usa su valor por defecto (static: N/p bloques contiguos; dynamic: 1). */
    omp_set_num_threads(threads_req);
    omp_set_schedule(sched_kind, chunk);

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

    /* 3. Inicializacion determinista de cuerpos con seed = 42 (secuencial,
     * para que las condiciones iniciales sean identicas a la version base) */
    init_bodies(bodies, n, seed);

    double initial_checksum = compute_checksum(bodies, n);
    printf("Inicializacion completada. Cuerpos: %d, Pasos: %d, Semilla: %u\n", n, steps, seed);
    printf("Checksum inicial (t=0): %.10e\n", initial_checksum);

    int threads_used = 1;

    /* 4. Medicion de tiempo de ejecucion con omp_get_wtime() */
    double start_time = omp_get_wtime();

    /* 5. Region paralela unica para toda la simulacion.
     *
     * Se crea el equipo de hilos una sola vez (un solo fork-join) en lugar de
     * abrir una region por cada paso, para no pagar el costo de creacion en
     * cada iteracion temporal.
     *
     * default(none) obliga a declarar el alcance de cada variable:
     *   shared : bodies, fx, fy, n, steps, threads_used
     *            - bodies se lee por todos durante la fase de fuerzas y cada
     *              hilo escribe solo los cuerpos de sus iteraciones en la fase
     *              de integracion.
     *            - fx[i], fy[i] los escribe unicamente el hilo al que le toco i.
     *   private: step, i, j y todas las variables temporales (dx, dy, ...) se
     *            declaran dentro de la region, por lo que cada hilo tiene su
     *            propia copia en su pila.
     */
    #pragma omp parallel default(none) shared(bodies, fx, fy, n, steps, threads_used)
    {
        /* Un solo hilo registra cuantos hilos entrego realmente el runtime */
        #pragma omp single
        threads_used = omp_get_num_threads();

        /* Todos los hilos recorren el ciclo temporal; el trabajo de cada paso
         * se reparte con los 'omp for' internos. */
        for (int step = 0; step < steps; step++) {
            /* FASE 1: Calculo de fuerzas O(N^2), paralelizado sobre i.
             *
             * schedule(runtime) toma el schedule y el chunk definidos con
             * omp_set_schedule(), lo que permite probar static/dynamic/guided
             * sin recompilar.
             *
             * Sin condiciones de carrera: la iteracion i solo escribe fx[i] y
             * fy[i]; bodies[] es de solo lectura en esta fase. El ciclo
             * interno sobre j se recorre en el mismo orden que la version
             * secuencial, por lo que el resultado es identico bit a bit.
             *
             * La barrera implicita al final del 'omp for' garantiza que todas
             * las fuerzas esten listas antes de mover cualquier cuerpo. */
            #pragma omp for schedule(runtime)
            for (int i = 0; i < n; i++) {
                double fxi = 0.0;   /* acumuladores privados del hilo */
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

            /* FASE 2: Integracion O(N), ciclo aparte.
             *
             * Cada iteracion modifica solo bodies[i] usando fx[i] y fy[i],
             * asi que tampoco hay carreras. El costo por cuerpo es constante,
             * por lo que schedule(static) reparte bien y no tiene overhead.
             *
             * La barrera implicita al final garantiza que ningun hilo empiece
             * el calculo de fuerzas del siguiente paso con posiciones a medio
             * actualizar. */
            #pragma omp for schedule(static)
            for (int i = 0; i < n; i++) {
                double ax = fx[i] / bodies[i].mass;
                double ay = fy[i] / bodies[i].mass;

                bodies[i].vx += ax * DT;
                bodies[i].vy += ay * DT;

                bodies[i].x += bodies[i].vx * DT;
                bodies[i].y += bodies[i].vy * DT;
            }
        }
    }

    double end_time = omp_get_wtime();
    double elapsed_time = end_time - start_time;

    /* 6. Checksum final (secuencial, mismo orden de suma que la version base) */
    double final_checksum = compute_checksum(bodies, n);

    /* 7. Resumen de resultados */
    print_summary(
        "nbody_paralelo (B)",
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
