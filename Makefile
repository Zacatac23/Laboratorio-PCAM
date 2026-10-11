# ==============================================================================
# Makefile - Laboratorio PCAM + OpenMP (N-Body 2D)
# Persona A: Base, correccion y ensamblado
# ==============================================================================

CC = gcc
CFLAGS = -O2 -Wall -fopenmp
LDFLAGS = -lm

BIN_DIR = bin
SRC_DIR = src

# Targets
all: directories $(BIN_DIR)/nbody_seq nbody_secuencial $(BIN_DIR)/nbody_parallel nbody_paralelo

directories:
	@mkdir -p $(BIN_DIR)

# Compilacion para la estructura interna de trabajo
$(BIN_DIR)/nbody_seq: $(SRC_DIR)/nbody_seq.c $(SRC_DIR)/common.h
	$(CC) $(CFLAGS) $(SRC_DIR)/nbody_seq.c -o $@ $(LDFLAGS)

# Compilacion del archivo directo para entrega en Canvas
nbody_secuencial: nbody_secuencial.c
	$(CC) $(CFLAGS) nbody_secuencial.c -o $@ $(LDFLAGS)

# Version paralela OpenMP (Persona B)
$(BIN_DIR)/nbody_parallel: $(SRC_DIR)/nbody_parallel.c $(SRC_DIR)/common.h
	$(CC) $(CFLAGS) $(SRC_DIR)/nbody_parallel.c -o $@ $(LDFLAGS)

nbody_paralelo: nbody_paralelo.c
	$(CC) $(CFLAGS) nbody_paralelo.c -o $@ $(LDFLAGS)

# Ejecucion de prueba con los parametros por defecto (N=5000, steps=10)
run_seq: $(BIN_DIR)/nbody_seq
	./$(BIN_DIR)/nbody_seq 5000 10 1 sequential 0

# Ejecucion rapida de prueba para validacion (N=1000, steps=10)
run_quick: $(BIN_DIR)/nbody_seq
	./$(BIN_DIR)/nbody_seq 1000 10 1 sequential 0

run_par: $(BIN_DIR)/nbody_parallel
	./$(BIN_DIR)/nbody_parallel 5000 10 8 static 0

# Bateria completa de experimentos (5 repeticiones) + tabla de speedup/efficiency
experiments:
	bash scripts/run_experiments.sh 5000 10 5 results.csv
	python scripts/analyze_results.py results.csv

clean:
	rm -rf $(BIN_DIR) nbody_secuencial nbody_paralelo *.o

.PHONY: all directories run_seq run_quick run_par experiments clean
