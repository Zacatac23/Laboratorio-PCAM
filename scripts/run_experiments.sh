#!/usr/bin/env bash
# ==============================================================================
# run_experiments.sh - Bateria de experimentos de rendimiento (Persona B)
# Laboratorio PCAM + OpenMP (N-Body 2D)
#
# Ejecuta cada configuracion de la seccion 5 del plan REPS veces y guarda cada
# corrida (no promedios) en results.csv. Los promedios, speedup y efficiency
# se calculan despues con scripts/analyze_results.py.
#
# Uso (desde la raiz del repo, Linux o Git Bash en Windows):
#   bash scripts/run_experiments.sh [N] [steps] [reps] [salida.csv]
#   Por defecto: N=5000, steps=10, reps=5 (minimo exigido: 3), salida=results.csv
#
# Si existe bin/nbody_optimized (Persona C), tambien se mide con las mismas
# configuraciones para la comparacion parallel vs. optimized.
#
# Recomendacion: cerrar otros programas, usar siempre la misma maquina y, en
# laptop, conectarla al cargador con el plan de energia de alto rendimiento.
# ==============================================================================
set -euo pipefail

N=${1:-5000}
STEPS=${2:-10}
REPS=${3:-5}
OUT=${4:-results.csv}

BIN=bin
CC=${CC:-gcc}
CFLAGS="-O2 -fopenmp"

mkdir -p "$BIN"

# Compilacion desde cero para medir siempre el codigo actual
echo "== Compilando =="
$CC $CFLAGS src/nbody_seq.c      -o "$BIN/nbody_seq"      -lm
$CC $CFLAGS src/nbody_parallel.c -o "$BIN/nbody_parallel" -lm
if [ -f src/nbody_optimized.c ]; then
    $CC $CFLAGS src/nbody_optimized.c -o "$BIN/nbody_optimized" -lm
fi

# Configuraciones: "threads schedule chunk" (chunk 0 = default de OpenMP)
CONFIGS=(
    "1 static 0"
    "2 static 0"
    "4 static 0"
    "8 static 0"
    "8 static 8"
    "8 static 64"
    "8 dynamic 8"
    "8 dynamic 64"
    "8 guided 0"
)

echo "rep,version,N,steps,seed,threads_req,threads_used,schedule,chunk,procs,max_threads,time_sec,checksum" > "$OUT"

# Ejecuta un binario y agrega su linea CSV_ROW (con el numero de repeticion).
# Reintenta: en Windows el antivirus puede bloquear un .exe recien compilado
# ("Permission denied") durante unos segundos.
run() {
    local rep=$1; shift
    local row="" attempt
    for attempt in 1 2 3 4 5; do
        row=$("$@" | grep '^CSV_ROW:' | sed 's/^CSV_ROW: //') && [ -n "$row" ] && break
        echo "  (fallo la ejecucion, reintento $attempt)" >&2
        sleep 2
    done
    [ -n "$row" ] || { echo "No se pudo ejecutar: $*" >&2; exit 1; }
    echo "$rep,$row" >> "$OUT"
    echo "  rep $rep | $row"
}

echo "== N=$N steps=$STEPS reps=$REPS -> $OUT =="

# Calentamiento descartado: estabiliza frecuencia del CPU y cache de disco
echo "-- calentamiento (no se registra)"
"./$BIN/nbody_seq" "$N" "$STEPS" 1 sequential 0 > /dev/null || true
"./$BIN/nbody_parallel" "$N" "$STEPS" 8 static 0 > /dev/null || true

VERSIONS=(nbody_parallel)
if [ -x "$BIN/nbody_optimized" ] || [ -x "$BIN/nbody_optimized.exe" ]; then
    VERSIONS+=(nbody_optimized)
fi

# Repeticiones intercaladas: cada ronda recorre todas las configuraciones
# (secuencial incluida). Asi la deriva de frecuencia/temperatura del CPU
# afecta a todas por igual y no solo a las primeras que se miden.
for r in $(seq 1 "$REPS"); do
    echo "-- ronda $r de $REPS"
    run "$r" "./$BIN/nbody_seq" "$N" "$STEPS" 1 sequential 0
    for v in "${VERSIONS[@]}"; do
        for cfg in "${CONFIGS[@]}"; do
            # shellcheck disable=SC2086
            run "$r" "./$BIN/$v" "$N" "$STEPS" $cfg
        done
    done
done

echo "== Listo. Resumen: python scripts/analyze_results.py $OUT =="
