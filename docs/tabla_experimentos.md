N = 5000, pasos = 10, seed = 42, T1 = 1.7878 s (promedio de 5 corridas secuenciales)

| Version | Threads | Schedule | Chunk | Tiempo prom. (s) | Desv. est. (s) | Speedup | Efficiency | Karp-Flatt | Error rel. checksum |
|---|---|---|---|---|---|---|---|---|---|
| nbody_secuencial | 1 | sequential | default | 1.7878 | 0.2449 | 1.00 | 1.00 | - | 0.0e+00 |
| nbody_paralelo | 1 | static | default | 1.7814 | 0.1678 | 1.00 | 1.00 | - | 0.0e+00 |
| nbody_paralelo | 2 | static | default | 1.0824 | 0.0872 | 1.65 | 0.83 | 0.211 | 0.0e+00 |
| nbody_paralelo | 4 | static | default | 0.7484 | 0.0545 | 2.39 | 0.60 | 0.225 | 0.0e+00 |
| nbody_paralelo | 8 | static | default | 0.5472 | 0.0253 | 3.27 | 0.41 | 0.207 | 0.0e+00 |
| nbody_paralelo | 8 | static | 8 | 0.5420 | 0.0310 | 3.30 | 0.41 | 0.204 | 0.0e+00 |
| nbody_paralelo | 8 | static | 64 | 0.5700 | 0.0476 | 3.14 | 0.39 | 0.222 | 0.0e+00 |
| nbody_paralelo | 8 | dynamic | 8 | 0.4740 | 0.0276 | 3.77 | 0.47 | 0.160 | 0.0e+00 |
| nbody_paralelo | 8 | dynamic | 64 | 0.5186 | 0.0253 | 3.45 | 0.43 | 0.189 | 0.0e+00 |
| nbody_paralelo | 8 | guided | default | 0.5368 | 0.0489 | 3.33 | 0.42 | 0.200 | 0.0e+00 |
