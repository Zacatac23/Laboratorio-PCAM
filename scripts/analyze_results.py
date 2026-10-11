"""
analyze_results.py - Promedios, speedup y efficiency a partir de results.csv
Laboratorio PCAM + OpenMP (N-Body 2D) - Persona B

Uso:
    python scripts/analyze_results.py [results.csv]

Genera:
    results_summary.csv          tabla resumen (una fila por configuracion)
    docs/tabla_experimentos.md   tabla en Markdown
    docs/tabla_experimentos.tex  tabla en LaTeX para el informe y evidencia.pdf

Definiciones:
    T1   = tiempo promedio de la version secuencial
    S(p) = T1 / Tp
    E(p) = S(p) / p   (p = threads usados)
    e(p) = (1/S - 1/p) / (1 - 1/p)   fraccion serial experimental (Karp-Flatt),
           solo para p > 1; si crece con p, el limite es overhead/hardware y
           no solo codigo serial
"""

import csv
import statistics
import sys
from collections import OrderedDict

src = sys.argv[1] if len(sys.argv) > 1 else "results.csv"

with open(src, newline="") as f:
    rows = list(csv.DictReader(f))

if not rows:
    sys.exit(f"{src} esta vacio")

# Agrupar repeticiones por (version, threads, schedule, chunk) en orden de aparicion
groups = OrderedDict()
for r in rows:
    key = (r["version"], int(r["threads_used"]), r["schedule"], int(r["chunk"]))
    groups.setdefault(key, []).append(r)

seq = [k for k in groups if k[2] == "sequential"]
if not seq:
    sys.exit("No hay corridas secuenciales para calcular T1")
seq_rows = groups[seq[0]]
t1 = statistics.mean(float(r["time_sec"]) for r in seq_rows)
ref_checksum = float(seq_rows[0]["checksum"])

summary = []
for (version, p, sched, chunk), rs in groups.items():
    times = [float(r["time_sec"]) for r in rs]
    mean = statistics.mean(times)
    sd = statistics.stdev(times) if len(times) > 1 else 0.0
    checksums = {float(r["checksum"]) for r in rs}
    max_err = max(abs(c - ref_checksum) / abs(ref_checksum) for c in checksums)
    summary.append({
        "version": version.split(" ")[0],
        "threads": p,
        "schedule": sched,
        "chunk": "default" if chunk == 0 else str(chunk),
        "reps": len(times),
        "t_mean": mean,
        "t_sd": sd,
        "speedup": t1 / mean,
        "efficiency": t1 / mean / p,
        "karp_flatt": (mean / t1 - 1 / p) / (1 - 1 / p) if p > 1 else None,
        "rel_err": max_err,
    })



def kf(s, fmt, empty):
    return empty if s["karp_flatt"] is None else fmt.format(s["karp_flatt"])


with open("results_summary.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["version", "threads", "schedule", "chunk", "reps",
                "time_mean_s", "time_sd_s", "speedup", "efficiency", "karp_flatt", "checksum_rel_err"])
    for s in summary:
        w.writerow([s["version"], s["threads"], s["schedule"], s["chunk"], s["reps"],
                    f"{s['t_mean']:.4f}", f"{s['t_sd']:.4f}", f"{s['speedup']:.3f}",
                    f"{s['efficiency']:.3f}", kf(s, "{:.4f}", ""), f"{s['rel_err']:.2e}"])

n, steps = rows[0]["N"], rows[0]["steps"]
md = [f"N = {n}, pasos = {steps}, seed = 42, T1 = {t1:.4f} s "
      f"(promedio de {len(seq_rows)} corridas secuenciales)\n",
      "| Version | Threads | Schedule | Chunk | Tiempo prom. (s) | Desv. est. (s) | Speedup | Efficiency | Karp-Flatt | Error rel. checksum |",
      "|---|---|---|---|---|---|---|---|---|---|"]
for s in summary:
    md.append(f"| {s['version']} | {s['threads']} | {s['schedule']} | {s['chunk']} | "
              f"{s['t_mean']:.4f} | {s['t_sd']:.4f} | {s['speedup']:.2f} | "
              f"{s['efficiency']:.2f} | {kf(s, '{:.3f}', '-')} | {s['rel_err']:.1e} |")
with open("docs/tabla_experimentos.md", "w", encoding="utf-8") as f:
    f.write("\n".join(md) + "\n")


def tex(s):
    return s.replace("_", r"\_")


lines = [r"\begin{table}[h!]", r"\centering",
         rf"\caption{{Tiempos promedio ({len(seq_rows)} repeticiones), speedup y eficiencia. "
         rf"$N={n}$, {steps} pasos, $T_1 = {t1:.3f}$ s}}",
         r"\label{tab:experimentos}", r"\small",
         r"\begin{tabular}{@{}lcccrrrrr@{}}", r"\toprule",
         r"\textbf{Versión} & \textbf{$p$} & \textbf{Schedule} & \textbf{Chunk} & "
         r"\textbf{$\bar{T}$ (s)} & \textbf{$\sigma$ (s)} & \textbf{$S(p)$} & \textbf{$E(p)$} & \textbf{$e_{KF}$} \\ \midrule"]
for s in summary:
    lines.append(f"{tex(s['version'])} & {s['threads']} & {s['schedule']} & {s['chunk']} & "
                 f"{s['t_mean']:.3f} & {s['t_sd']:.3f} & {s['speedup']:.2f} & {s['efficiency']:.2f} & {kf(s, '{:.3f}', '--')} \\\\")
lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
with open("docs/tabla_experimentos.tex", "w", encoding="utf-8") as f:
    f.write("\n".join(lines) + "\n")

print("\n".join(md))
print("\nEscrito: results_summary.csv, docs/tabla_experimentos.md, docs/tabla_experimentos.tex")
