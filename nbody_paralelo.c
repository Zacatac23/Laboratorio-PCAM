/**
 * ============================================================================
 * nbody_paralelo.c - Implementacion Paralela con OpenMP (entrega Canvas)
 * Laboratorio PCAM + OpenMP
 * Responsable: Persona B
 *
 * Compilacion:
 *   gcc -O2 -fopenmp nbody_paralelo.c -o nbody_paralelo -lm
 *
 * Uso:
 *   ./nbody_paralelo [N] [steps] [threads] [schedule] [chunk]
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
