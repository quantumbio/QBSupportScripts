#!/usr/bin/env python3
"""Summarize per-PDB qmechanic perf.data CPU samples across MDComparison jobs.

Requires Linux perf on the machine that recorded the profiles.  This analyzer
reads profiles; it does not run or modify the MD calculations.

Example:
    python summarize_qmechanic_perf.py . --top 25

Outputs (relative to --out-dir):
    qmechanic_perf_summary.txt       human-readable aggregate report
    qmechanic_perf_jobs.csv          per-PDB sample counts and coverage
    qmechanic_perf_hotspots.csv      cross-job hotspots (per-job and pooled %)
    qmechanic_perf_libraries.csv     DSO/library-level CPU sample shares
    qmechanic_perf_job_hotspots.csv  top functions for each individual job
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path
import re
import shutil
import subprocess
import sys

# perf script -F time,event,ip,sym,dso prints one sample header per event:
#  12345.678901: cpu-clock: 7f123456 my_function (/path/libfoo.so)
# Following callchain frames have no timestamp/event prefix and are ignored.
SAMPLE = re.compile(
    r"^\s*(?P<time>\d+\.\d+):\s+"
    r"(?P<event>[^:]+):\s+"
    r"(?P<ip>(?:0x)?[0-9a-fA-F]+)\s+"
    r"(?P<symbol>.*?)\s+\((?P<dso>[^\r\n()]*)\)\s*$"
)


def summarize_profile(path: Path, perf_bin: str):
    """Stream perf-script sample headers; never buffer full callchains in RAM."""
    cmd = [perf_bin, "script", "-i", str(path), "-F", "time,event,ip,sym,dso"]
    symbol_counts = Counter()
    library_counts = Counter()
    sample_count = 0
    unmatched_headers = 0
    error = None
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                text=True, errors="replace", bufsize=1)
        assert proc.stdout is not None
        for line in proc.stdout:
            m = SAMPLE.match(line)
            if m:
                # Ignore stack-trace frames: only the record with time/event is sampled.
                symbol = m.group("symbol").strip() or "[unknown]"
                dso = m.group("dso").strip() or "[unknown]"
                symbol_counts[(dso, symbol)] += 1
                library_counts[dso] += 1
                sample_count += 1
            elif re.match(r"^\s*\d+\.\d+:\s+", line):
                unmatched_headers += 1
        proc.stdout.close()
        stderr = proc.stderr.read() if proc.stderr is not None else ""
        rc = proc.wait()
        if rc != 0:
            error = "perf script failed (exit {}): {}".format(rc, " ".join(stderr.strip().split())[:350])
        elif sample_count == 0:
            error = "No sample headers recognized; inspect 'perf script -i FILE -F time,event,ip,sym,dso | head -30'"
        elif unmatched_headers:
            error = "{} unrecognized sample headers; check output format".format(unmatched_headers)
    except (OSError, ValueError) as exc:
        error = str(exc)
    return sample_count, symbol_counts, library_counts, error


def write_csv(path, fieldnames, rows):
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def pct(numerator, denominator):
    return round(100 * numerator / denominator, 4) if denominator else 0.0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", nargs="?", default=".", type=Path,
                    help="Directory containing <PDBID>/cpp/qmechanic.perf.data (default: .)")
    ap.add_argument("--out-dir", type=Path, default=None, help="Report output directory (default: root)")
    ap.add_argument("--top", type=int, default=25, help="Number of top hotspots in text report (default: 25)")
    ap.add_argument("--perf", default="perf", help="Path to perf executable (default: perf)")
    args = ap.parse_args(argv)
    root = args.root.expanduser().resolve()
    out = (args.out_dir or root).expanduser().resolve()
    if not root.is_dir():
        ap.error("root directory does not exist: {}".format(root))
    if not shutil.which(args.perf):
        ap.error("perf executable not found: {}".format(args.perf))
    if args.top < 1:
        ap.error("--top must be positive")
    out.mkdir(parents=True, exist_ok=True)

    paths = sorted(root.glob("*/cpp/qmechanic.perf.data"), key=lambda p: p.parent.parent.name)
    if not paths:
        print("No profiles found under {}/<PDBID>/cpp/qmechanic.perf.data".format(root), file=sys.stderr)
        return 2

    jobs = []
    symbols_by_job = {}
    libs_by_job = {}
    for path in paths:
        pdb = path.parent.parent.name
        print("Reading {} ...".format(path.relative_to(root)), file=sys.stderr, flush=True)
        n, symbols, libs, error = summarize_profile(path, args.perf)
        jobs.append({"pdbid": pdb, "profile": str(path.relative_to(root)),
                     "samples": n, "distinct_symbols": len(symbols),
                     "distinct_libraries": len(libs), "status": "OK" if not error else "ERROR",
                     "error": error or ""})
        if not error:
            symbols_by_job[pdb] = symbols
            libs_by_job[pdb] = libs
        if error:
            print("  WARNING: {}".format(error), file=sys.stderr)

    valid_jobs = len(symbols_by_job)
    total_samples = sum(sum(c.values()) for c in symbols_by_job.values())
    write_csv(out / "qmechanic_perf_jobs.csv",
              ["pdbid", "profile", "samples", "distinct_symbols", "distinct_libraries", "status", "error"], jobs)
    if not valid_jobs:
        print("No valid profiles parsed. See qmechanic_perf_jobs.csv.", file=sys.stderr)
        return 1

    total_symbols = Counter()
    total_libs = Counter()
    per_job_share = defaultdict(float)
    per_job_lib_share = defaultdict(float)
    present_jobs = Counter()
    lib_present_jobs = Counter()
    job_hotspots = []
    for pdb, counter in sorted(symbols_by_job.items()):
        n = sum(counter.values())
        total_symbols.update(counter)
        for key, count in counter.items():
            per_job_share[key] += count / n
            present_jobs[key] += 1
        for (dso, symbol), count in counter.most_common():
            job_hotspots.append({"pdbid": pdb, "dso": dso, "symbol": symbol,
                                 "samples": count, "sample_pct": pct(count, n)})
        total_libs.update(libs_by_job[pdb])
        for dso, count in libs_by_job[pdb].items():
            per_job_lib_share[dso] += count / n
            lib_present_jobs[dso] += 1

    # Equal-job mean: missing function in a job contributes 0% to that job.
    # Pooled mean: each sample is weighted equally, so long/CPU-heavy jobs dominate.
    hotspot_rows = []
    for (dso, symbol), count in total_symbols.items():
        hotspot_rows.append({"dso": dso, "symbol": symbol, "jobs_with_samples": present_jobs[(dso, symbol)],
                             "jobs_total": valid_jobs, "samples": count,
                             "mean_job_pct": round(per_job_share[(dso, symbol)] * 100 / valid_jobs, 4),
                             "pooled_sample_pct": pct(count, total_samples)})
    hotspot_rows.sort(key=lambda r: (-r["mean_job_pct"], -r["samples"], r["dso"], r["symbol"]))
    library_rows = []
    for dso, count in total_libs.items():
        library_rows.append({"dso": dso, "jobs_with_samples": lib_present_jobs[dso],
                             "jobs_total": valid_jobs, "samples": count,
                             "mean_job_pct": round(per_job_lib_share[dso] * 100 / valid_jobs, 4),
                             "pooled_sample_pct": pct(count, total_samples)})
    library_rows.sort(key=lambda r: (-r["mean_job_pct"], -r["samples"], r["dso"]))

    write_csv(out / "qmechanic_perf_hotspots.csv",
              ["dso", "symbol", "jobs_with_samples", "jobs_total", "samples", "mean_job_pct", "pooled_sample_pct"], hotspot_rows)
    write_csv(out / "qmechanic_perf_libraries.csv",
              ["dso", "jobs_with_samples", "jobs_total", "samples", "mean_job_pct", "pooled_sample_pct"], library_rows)
    write_csv(out / "qmechanic_perf_job_hotspots.csv",
              ["pdbid", "dso", "symbol", "samples", "sample_pct"], job_hotspots)

    lines = [
        "qmechanic CPU sampling profile summary",
        "====================================",
        "Root: {}".format(root),
        "Profiles found: {}; parsed OK: {}; failed: {}".format(len(paths), valid_jobs, len(paths)-valid_jobs),
        "Samples across valid profiles: {:,}".format(total_samples),
        "",
        "INTERPRETATION",
        "  mean_job_pct: each PDB job contributes equally (preferred for cross-system hotspots)",
        "  pooled_sample_pct: all samples pooled (CPU-heavy/long jobs have more influence)",
        "  Percentages represent ON-CPU samples, NOT seconds or wall-clock overhead.",
        "  A function appearing in a call chain does not count unless it was the sampled IP.",
        "  The profile includes the entire C++ MD qmechanic invocation (not prepare).",
        "  These data alone cannot distinguish time inside vs. outside MD stages.",
        "",
        "PER-JOB COVERAGE",
        "{:<14} {:>12} {:>12}  {}".format("PDBID", "Samples", "Symbols", "Status"),
    ]
    for job in jobs:
        lines.append("{:<14} {:>12,} {:>12,}  {}{}".format(
            job["pdbid"], job["samples"], job["distinct_symbols"], job["status"],
            (" - " + job["error"]) if job["error"] else ""))
    lines.extend(["", "TOP LIBRARIES (equal-job average)",
                  "{:>9} {:>9} {:>8}  {}".format("Mean %", "Pooled %", "Jobs", "Library")])
    for row in library_rows[:args.top]:
        lines.append("{:9.2f} {:9.2f} {:>3}/{:<4}  {}".format(
            row["mean_job_pct"], row["pooled_sample_pct"], row["jobs_with_samples"], valid_jobs, row["dso"]))
    lines.extend(["", "TOP FUNCTIONS (equal-job average)",
                  "{:>9} {:>9} {:>8}  {}".format("Mean %", "Pooled %", "Jobs", "Function [library]")])
    for row in hotspot_rows[:args.top]:
        lines.append("{:9.2f} {:9.2f} {:>3}/{:<4}  {} [{}]".format(
            row["mean_job_pct"], row["pooled_sample_pct"], row["jobs_with_samples"], valid_jobs,
            row["symbol"], row["dso"]))
    lines.extend(["", "NEXT STEP",
                  "  Use the existing per-stage wall-clock timers to isolate the ~200 s non-MD work.",
                  "  For phase attribution, capture perf for the same run with phase-specific timers",
                  "  or use perf script timestamps aligned with MD stage boundaries.", ""])
    text_path = out / "qmechanic_perf_summary.txt"
    text_path.write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    print("\nWrote reports to: {}".format(out))
    return 0 if len(paths) == valid_jobs else 1


if __name__ == "__main__":
    sys.exit(main())
