/**
 * ============================================================================
 * nbody_optimized.c - Version Optimizada por Simetria de Pares (N-Body 2D)
 * Laboratorio PCAM + OpenMP
 * Responsable: Persona C
 *
 * Compilacion:
 *   gcc -O2 -fopenmp src/nbody_optimized.c -o bin/nbody_optimized -lm
 *
 * Uso:
 *   ./bin/nbody_optimized [N] [steps] [threads] [schedule] [chunk] [estrategia]
 *
 *   schedule   : static | dynamic | guided   (por defecto static)
 *   chunk      : tamano de bloque; 0 = valor por defecto de OpenMP
 *   estrategia : privado | atomic | critical (por defecto privado)
 *
 * Idea: por la tercera ley de Newton F(j sobre i) = -F(i sobre j), asi que
 * cada par (i, j) con j > i se calcula UNA sola vez y se suma a los dos
 * cuerpos: N(N-1)/2 interacciones por paso en lugar de N(N-1).
 *
 * Race condition: la iteracion i escribe la fuerza de i y tambien la de
 * todos los j > i. Dos hilos con filas distintas (i1 < i2 < j) escriben la
 * fuerza del mismo cuerpo j al mismo tiempo, y "fx[j] -= f" es una
 * lectura-modificacion-escritura no atomica: se perderian contribuciones.
 * Tres formas de resolverlo, seleccionables para compararlas con mediciones:
 *   - privado : cada hilo acumula en su propia fila fxp[t][*]; al final del
 *               paso se suman las filas. No hay escrituras compartidas
 *               durante el calculo (estrategia elegida).
 *   - atomic  : #pragma omp atomic sobre cada escritura compartida.
 *   - critical: #pragma omp critical sobre cada escritura compartida.
 * ============================================================================
 */

#include "common.h"

typedef enum { ACUM_PRIVADO, ACUM_ATOMIC, ACUM_CRITICAL } Estrategia;

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

static int parse_estrategia(const char *name, Estrategia *e) {
    if (strcmp(name, "privado") == 0)  { *e = ACUM_PRIVADO;  return 1; }
    if (strcmp(name, "atomic") == 0)   { *e = ACUM_ATOMIC;   return 1; }
    if (strcmp(name, "critical") == 0) { *e = ACUM_CRITICAL; return 1; }
    return 0;
}

/**
 * Fuerza que el cuerpo b ejerce sobre a. La que a ejerce sobre b es la
 * opuesta. Misma formula y mismo orden de operaciones que la secuencial.
 */
static inline void fuerza_par(const Body *a, const Body *b, double *fx, double *fy) {
    double dx = b->x - a->x;
    double dy = b->y - a->y;
    double dist_sqr = dx * dx + dy * dy + SOFTENING;
    double dist = sqrt(dist_sqr);

    double inv_dist3 = 1.0 / (dist * dist_sqr);
    double force_mag = G_CONST * a->mass * b->mass * inv_dist3;

    *fx = force_mag * dx;
    *fy = force_mag * dy;
}

int main(int argc, char *argv[]) {
    /* 1. Lectura de argumentos con valores por defecto acordados */
    int n = DEFAULT_N;
    int steps = DEFAULT_STEPS;
    int threads_req = omp_get_max_threads();
    const char *schedule = "static";
    int chunk = 0;
    const char *estrategia_txt = "privado";
    unsigned int seed = DEFAULT_SEED;

    if (argc > 1) n = atoi(argv[1]);
    if (argc > 2) steps = atoi(argv[2]);
    if (argc > 3) threads_req = atoi(argv[3]);
    if (argc > 4) schedule = argv[4];
    if (argc > 5) chunk = atoi(argv[5]);
    if (argc > 6) estrategia_txt = argv[6];

    omp_sched_t sched_kind;
    Estrategia estrategia;
    if (n <= 0 || steps < 0 || threads_req <= 0 || chunk < 0 ||
        !parse_schedule(schedule, &sched_kind) ||
        !parse_estrategia(estrategia_txt, &estrategia)) {
        fprintf(stderr, "Uso: %s [N] [steps] [threads] [static|dynamic|guided] [chunk] "
                        "[privado|atomic|critical]\n", argv[0]);
        return EXIT_FAILURE;
    }

    omp_set_num_threads(threads_req);
    omp_set_schedule(sched_kind, chunk);

    /* 2. Reserva de memoria dinamica.
     * privado: una fila de N acumuladores por hilo (threads_req x N). El
     *          runtime nunca entrega mas hilos que los pedidos, asi que el
     *          indice omp_get_thread_num() siempre cabe.
     * atomic/critical: un solo par de arreglos fx, fy compartido. */
    int filas = (estrategia == ACUM_PRIVADO) ? threads_req : 1;
    Body *bodies = (Body *)malloc(sizeof(Body) * n);
    double *fx = (double *)malloc(sizeof(double) * (size_t)filas * n);
    double *fy = (double *)malloc(sizeof(double) * (size_t)filas * n);

    if (!bodies || !fx || !fy) {
        fprintf(stderr, "Error: No se pudo asignar memoria para %d cuerpos.\n", n);
        free(bodies);
        free(fx);
        free(fy);
        return EXIT_FAILURE;
    }

    /* 3. Inicializacion determinista de cuerpos con seed = 42 */
    init_bodies(bodies, n, seed);

    double initial_checksum = compute_checksum(bodies, n);
    printf("Inicializacion completada. Cuerpos: %d, Pasos: %d, Semilla: %u\n", n, steps, seed);
    printf("Checksum inicial (t=0): %.10e\n", initial_checksum);

    int threads_used = 1;

    /* 4. Medicion de tiempo de ejecucion con omp_get_wtime() */
    double start_time = omp_get_wtime();

    /* 5. Region paralela unica para toda la simulacion (un solo fork-join).
     *
     *   shared : bodies, fx, fy, n, steps, estrategia, threads_used
     *            - bodies: solo lectura en la fase de fuerzas; en la fase de
     *              integracion cada hilo escribe solo los cuerpos de sus
     *              iteraciones.
     *            - fx, fy: con 'privado' cada hilo escribe solo su fila;
     *              con atomic/critical todos escriben en el mismo arreglo,
     *              y cada escritura va protegida.
     *   private: t, mfx, mfy, step, i, j y los temporales, declarados dentro
     *            de la region (cada hilo tiene su copia en la pila).
     */
    #pragma omp parallel default(none) shared(bodies, fx, fy, n, steps, estrategia, threads_used)
    {
        int t = omp_get_thread_num();
        double *mfx = fx + (size_t)t * n;   /* fila propia (solo 'privado') */
        double *mfy = fy + (size_t)t * n;

        #pragma omp single
        threads_used = omp_get_num_threads();
        /* barrera implicita del single: threads_used ya es valido para todos */

        for (int step = 0; step < steps; step++) {
            /* FASE 0: poner en cero los acumuladores.
             * privado: cada hilo limpia solo su fila; no hace falta barrera,
             *          porque en la fase 1 cada hilo solo toca su fila.
             * compartido: hace falta la barrera implicita del 'omp for',
             *          porque en la fase 1 cualquier hilo escribe cualquier j. */
            if (estrategia == ACUM_PRIVADO) {
                memset(mfx, 0, sizeof(double) * n);
                memset(mfy, 0, sizeof(double) * n);
            } else {
                #pragma omp for schedule(static)
                for (int k = 0; k < n; k++) {
                    fx[k] = 0.0;
                    fy[k] = 0.0;
                }
            }

            /* FASE 1: fuerzas por pares (i, j), j > i. Ciclo triangular: la
             * fila i tiene N-1-i interacciones, asi que la carga NO es
             * uniforme y el schedule importa (static por defecto le da al
             * hilo 0 las filas mas largas). schedule(runtime) toma el
             * schedule y chunk de omp_set_schedule().
             *
             * Las tres estrategias son ciclos separados (no un 'if' dentro
             * del ciclo interno) para que la verificacion de la estrategia
             * no quede en el camino critico. Todos los hilos ven el mismo
             * valor de 'estrategia', asi que todos encuentran el mismo
             * 'omp for', como exige OpenMP. */
            if (estrategia == ACUM_PRIVADO) {
                #pragma omp for schedule(runtime)
                for (int i = 0; i < n; i++) {
                    double fxi = 0.0, fyi = 0.0;
                    for (int j = i + 1; j < n; j++) {
                        double gx, gy;
                        fuerza_par(&bodies[i], &bodies[j], &gx, &gy);
                        fxi += gx;          /* j atrae a i */
                        fyi += gy;
                        mfx[j] -= gx;       /* i atrae a j: fila propia, sin carrera */
                        mfy[j] -= gy;
                    }
                    mfx[i] += fxi;
                    mfy[i] += fyi;
                }
            } else if (estrategia == ACUM_ATOMIC) {
                #pragma omp for schedule(runtime)
                for (int i = 0; i < n; i++) {
                    double fxi = 0.0, fyi = 0.0;
                    for (int j = i + 1; j < n; j++) {
                        double gx, gy;
                        fuerza_par(&bodies[i], &bodies[j], &gx, &gy);
                        fxi += gx;
                        fyi += gy;
                        #pragma omp atomic
                        fx[j] -= gx;
                        #pragma omp atomic
                        fy[j] -= gy;
                    }
                    /* fx[i] tambien recibe escrituras de filas i' < i */
                    #pragma omp atomic
                    fx[i] += fxi;
                    #pragma omp atomic
                    fy[i] += fyi;
                }
            } else {
                #pragma omp for schedule(runtime)
                for (int i = 0; i < n; i++) {
                    double fxi = 0.0, fyi = 0.0;
                    for (int j = i + 1; j < n; j++) {
                        double gx, gy;
                        fuerza_par(&bodies[i], &bodies[j], &gx, &gy);
                        fxi += gx;
                        fyi += gy;
                        #pragma omp critical
                        {
                            fx[j] -= gx;
                            fy[j] -= gy;
                        }
                    }
                    #pragma omp critical
                    {
                        fx[i] += fxi;
                        fy[i] += fyi;
                    }
                }
            }
            /* Barrera implicita: todas las contribuciones estan escritas. */

            /* FASE 2: reduccion de filas + integracion, en un solo ciclo.
             * La fuerza total de k es la suma de las filas de todos los
             * hilos, en orden fijo t = 0..p-1. Fusionarla con la
             * integracion evita una barrera extra por paso. Cada iteracion
             * escribe solo bodies[k]: sin carreras. */
            #pragma omp for schedule(static)
            for (int k = 0; k < n; k++) {
                double sx, sy;
                if (estrategia == ACUM_PRIVADO) {
                    sx = 0.0;
                    sy = 0.0;
                    for (int r = 0; r < threads_used; r++) {
                        sx += fx[(size_t)r * n + k];
                        sy += fy[(size_t)r * n + k];
                    }
                } else {
                    sx = fx[k];
                    sy = fy[k];
                }

                double ax = sx / bodies[k].mass;
                double ay = sy / bodies[k].mass;

                bodies[k].vx += ax * DT;
                bodies[k].vy += ay * DT;

                bodies[k].x += bodies[k].vx * DT;
                bodies[k].y += bodies[k].vy * DT;
            }
            /* Barrera implicita: nadie limpia su fila ni empieza el paso
             * siguiente mientras otro hilo todavia la esta sumando. */
        }
    }

    double end_time = omp_get_wtime();
    double elapsed_time = end_time - start_time;

    /* 6. Checksum final (secuencial, mismo orden de suma que la version base) */
    double final_checksum = compute_checksum(bodies, n);

    /* 7. Resumen. La estrategia va en el nombre de la version para que
     * results.csv distinga las tres variantes. */
    char version_name[64];
    snprintf(version_name, sizeof version_name, "nbody_optimizado_%s (C)", estrategia_txt);
    printf("Estrategia de acumulacion: %s | Interacciones por paso: %lld\n",
           estrategia_txt, (long long)n * (n - 1) / 2);
    print_summary(
        version_name,
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
