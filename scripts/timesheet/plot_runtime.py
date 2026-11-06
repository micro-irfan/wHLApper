#!/usr/bin/env python3
"""
Plot per-sample runtime of HLA*LA, SpecImmune and FuFiHLA on one panel
(box + points, lines join the same sample across tools, log y-axis).

Inputs (defaults are the files next to this script):
  trace.HLALA.txt  Nextflow trace; realtime of MINIMAP2_ALIGN + HLALA (EXTRACT_HLA_REGION excluded)
  SpecImmune.log   "[ts] SAMPLE | step | starting/done"; HLA step (add --include-align for the hg38 alignment)
  FuFiHLA.log      "[ts] SAMPLE | starting/done"

Usage:
  python plot_runtime.py [-o runtime_per_sample.png] [--include-align]
"""

import argparse
import csv
import re
import statistics
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

HERE = Path(__file__).resolve().parent

HLALA_STEPS = ("MINIMAP2_ALIGN", "HLALA")
TOOLS = ["HLA*LA", "SpecImmune", "FuFiHLA"]
COLORS = {"HLA*LA": "#4C72B0", "SpecImmune": "#DD8452", "FuFiHLA": "#55A868"}

LOG_LINE = re.compile(r"^\[(?P<ts>[\d-]+ [\d:]+)\]\s+(?P<sample>\S+)\s+\|\s+(?P<rest>.*)$")
TS_FMT = "%Y-%m-%d %H:%M:%S"


def parse_duration(s):
    """Nextflow duration ('1h 2m 7s', '42.1s', '565ms') -> seconds."""
    units = {"d": 86400, "h": 3600, "m": 60, "s": 1, "ms": 0.001}
    parts = re.findall(r"([\d.]+)\s*(ms|d|h|m|s)", s)
    if not parts:
        raise ValueError(f"Unrecognised duration: {s!r}")
    return sum(float(v) * units[u] for v, u in parts)


def read_trace(path):
    """Sum realtime of the HLA*LA steps per sample (cached tasks keep their original realtime)."""
    times = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["status"] not in ("COMPLETED", "CACHED"):
                continue
            m = re.match(r"^(\S+) \((.+)\)$", row["name"])
            if not m or m.group(1) not in HLALA_STEPS:
                continue
            sample = m.group(2)
            times[sample] = times.get(sample, 0.0) + parse_duration(row["realtime"])
    return times


def read_log(path, steps=None):
    """
    Wall time per sample from '[ts] SAMPLE | [step |] message' lines.
    The first line of a (sample, step) is its start, the line whose message
    starts with 'done' is its end. steps=None means lines have no step field.
    """
    start, total = {}, {}
    with open(path) as fh:
        for line in fh:
            m = LOG_LINE.match(line.strip())
            if not m:
                continue
            fields = [f.strip() for f in m.group("rest").split("|")]
            if steps is None:
                step, msg = "all", fields[0]
            else:
                if len(fields) < 2 or fields[0] not in steps:
                    continue
                step, msg = fields[0], fields[1]
            ts = datetime.strptime(m.group("ts"), TS_FMT)
            key = (m.group("sample"), step)
            start.setdefault(key, ts)
            if msg.startswith("done"):
                total[m.group("sample")] = total.get(m.group("sample"), 0.0) + (ts - start.pop(key)).total_seconds()
    return total


def fmt_time(sec):
    return f"{sec / 60:.1f} min" if sec >= 60 else f"{sec:.0f} s"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trace", default=HERE / "trace.HLALA.txt", type=Path)
    ap.add_argument("--specimmune", default=HERE / "SpecImmune.log", type=Path)
    ap.add_argument("--fufihla", default=HERE / "FuFiHLA.log", type=Path)
    ap.add_argument("--include-align", action="store_true",
                    help="Add SpecImmune's minimap2 alignment to hg38 to its runtime")
    ap.add_argument("-o", "--out", default=HERE / "runtime_per_sample.png", type=Path,
                    help="Output figure (.png/.pdf/.svg); a .tsv of the values is written alongside")
    args = ap.parse_args()

    spec_steps = ("align", "HLA") if args.include_align else ("HLA",)
    data = {
        "HLA*LA": read_trace(args.trace),
        "SpecImmune": read_log(args.specimmune, steps=spec_steps),
        "FuFiHLA": read_log(args.fufihla),
    }
    for tool in TOOLS:
        if not data[tool]:
            raise SystemExit(f"No runtimes parsed for {tool}")

    samples = sorted(set().union(*(d.keys() for d in data.values())))

    # Table of values used in the plot
    tsv = args.out.with_suffix(".tsv")
    with open(tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample_id"] + [f"{t}_seconds" for t in TOOLS])
        for s in samples:
            w.writerow([s] + [f"{data[t][s]:.1f}" if s in data[t] else "" for t in TOOLS])

    fig, ax = plt.subplots(figsize=(3.6, 4.2))
    xs = range(len(TOOLS))

    # Lines joining the same sample across tools
    for s in samples:
        pts = [(x, data[t][s]) for x, t in zip(xs, TOOLS) if s in data[t]]
        if len(pts) > 1:
            ax.plot(*zip(*pts), color="0.75", lw=0.8, zorder=1)

    for x, tool in zip(xs, TOOLS):
        vals = list(data[tool].values())
        ax.boxplot(vals, positions=[x], widths=0.5, showfliers=False, patch_artist=True,
                   boxprops=dict(facecolor=COLORS[tool], alpha=0.6, edgecolor="0.2"),
                   medianprops=dict(color="0.1", lw=1.5),
                   whiskerprops=dict(color="0.2"), capprops=dict(color="0.2"), zorder=2)
        ax.scatter([x] * len(vals), vals, s=22, color=COLORS[tool], edgecolor="white", lw=0.5, zorder=3)

    medians = {t: statistics.median(data[t].values()) for t in TOOLS}
    ax.set_xticks(list(xs))
    ax.set_xticklabels([f"{t}\n({fmt_time(medians[t])})" for t in TOOLS], fontsize=9)
    ax.set_yscale("log")
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.set_ylabel("Runtime per sample (s, log scale)")
    ax.set_title("Median runtime per sample", fontsize=10, style="italic")
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    fig.savefig(args.out, dpi=300)
    print(f"Wrote {args.out} and {tsv}")
    for t in TOOLS:
        print(f"  {t:<11} n={len(data[t])}  median={fmt_time(medians[t])}")


if __name__ == "__main__":
    main()
