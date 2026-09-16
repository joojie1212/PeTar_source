#!/usr/bin/env python3
"""Plot merger trajectories from actual petar.hard.debug output samples."""

import argparse
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np

PC_TO_RSUN = 3.0856775814913673e16 / 6.957e8


def read_event_pair(sidecar):
    with sidecar.open() as stream:
        for line in stream:
            f = line.split()
            if f and f[0] == "Dynamic_merge:":
                return int(f[1]), int(f[2]), "Dynamic_merge"
            if f and f[0] == "Coalescence":
                return int(f[-4]), int(f[-3]), "BSE_Coalescence"
    raise RuntimeError(f"No merger event in {sidecar}")


def event_properties(bse_path, id1, id2, petar):
    # Reuse the validated event parser used by the reconstructed plots.
    from plot_merger_orbits import load_events
    candidates = load_events(bse_path, petar)
    for event in candidates:
        if {event["id1"], event["id2"]} == {id1, id2}:
            return event
    raise RuntimeError(f"Cannot match merger IDs {id1}, {id2} in {bse_path}")


def load_hard_data(log, err, petar):
    text = err.read_text()
    match = re.search(r"Hard: n_ptcl: (\d+) n_group: (\d+)", text)
    if not match:
        raise RuntimeError(f"Cannot find particle count in {err}")
    n_particle = int(match.group(1))
    columns = np.loadtxt(log, skiprows=1, max_rows=1).size
    n_sd = None
    for trial in range(16):
        test = petar.HardData(interrupt_mode="bse", N_particle=n_particle,
                              N_sd=trial)
        if test.ncols == columns:
            n_sd = trial
            break
    if n_sd is None:
        raise RuntimeError(f"Cannot infer N_sd for {log}: {columns} columns")
    hard = petar.HardData(interrupt_mode="bse", N_particle=n_particle,
                          N_sd=n_sd)
    hard.loadtxt(log, skiprows=1)
    return hard, n_particle, n_sd


def extract_pair(hard, n_particle, id1, id2):
    times, pos1, pos2, vel1, vel2, mass1, mass2 = ([] for _ in range(7))
    for row in range(hard.size):
        found = {}
        for j in range(n_particle):
            particle = getattr(hard.particles, f"p{j}")
            pid = int(particle.id[row])
            if pid in (id1, id2) and particle.mass[row] > 0:
                found[pid] = (particle.pos[row].copy(), particle.vel[row].copy(),
                              float(particle.mass[row]))
        if id1 in found and id2 in found:
            times.append(float(hard.time_org[row]))
            pos1.append(found[id1][0]); pos2.append(found[id2][0])
            vel1.append(found[id1][1]); vel2.append(found[id2][1])
            mass1.append(found[id1][2]); mass2.append(found[id2][2])
    return tuple(np.asarray(x) for x in
                 (times, pos1, pos2, vel1, vel2, mass1, mass2))


def project_orbital_plane(pos1, pos2, vel1, vel2, mass1, mass2):
    rel = pos2 - pos1
    relv = vel2 - vel1
    # Define a fixed approximate plane from the last available real sample.
    xhat = rel[-1] / np.linalg.norm(rel[-1])
    h = np.cross(rel[-1], relv[-1])
    if np.linalg.norm(h) < 1e-30:
        _, _, vh = np.linalg.svd(rel - rel.mean(axis=0), full_matrices=False)
        zhat = vh[-1]
    else:
        zhat = h / np.linalg.norm(h)
    yhat = np.cross(zhat, xhat)
    yhat /= np.linalg.norm(yhat)
    total = mass1 + mass2
    cm = (mass1[:, None] * pos1 + mass2[:, None] * pos2) / total[:, None]
    q1 = pos1 - cm
    q2 = pos2 - cm
    xy1 = np.column_stack((q1 @ xhat, q1 @ yhat)) * PC_TO_RSUN
    xy2 = np.column_stack((q2 @ xhat, q2 @ yhat)) * PC_TO_RSUN
    return xy1, xy2


def draw(event, times, xy1, xy2, number, output, source, dpi):
    colors = "#277da1", "#f3722c"
    fig, ax = plt.subplots(figsize=(8.4, 7.0))
    ax.set_aspect("equal", adjustable="box")
    ax.axhline(0, color="0.88", lw=0.7, zorder=0)
    ax.axvline(0, color="0.88", lw=0.7, zorder=0)
    ax.plot(xy1[:, 0], xy1[:, 1], color=colors[0], lw=1.7,
            marker=".", ms=2.5, label=f"ID {event['id1']}")
    ax.plot(xy2[:, 0], xy2[:, 1], color=colors[1], lw=1.7,
            marker=".", ms=2.5, label=f"ID {event['id2']}")

    # Every arrow joins two actual consecutive debug output samples.  Limit
    # their number only to keep dense tracks readable.
    stride = max(1, (len(times) - 1) // 16)
    for track, color in ((xy1, colors[0]), (xy2, colors[1])):
        for j in range(stride, len(track), stride):
            ax.annotate("", xy=track[j], xytext=track[j-stride],
                        arrowprops=dict(arrowstyle="->", color=color,
                                        lw=0.9, alpha=0.72))

    final1, final2 = xy1[-1], xy2[-1]
    final_separation = float(np.linalg.norm(final2 - final1))
    ax.add_patch(Circle(final1, event["r1"], facecolor=colors[0],
                        edgecolor="black", alpha=0.8, zorder=5))
    ax.add_patch(Circle(final2, event["r2"], facecolor=colors[1],
                        edgecolor="black", alpha=0.8, zorder=5))
    ax.annotate(f"ID {event['id1']}\n$R_1={event['r1']:.4g}\,R_\odot$",
                final1, xytext=(-12, 24), textcoords="offset points",
                ha="right", arrowprops=dict(arrowstyle="->", color=colors[0]))
    ax.annotate(f"ID {event['id2']}\n$R_2={event['r2']:.4g}\,R_\odot$",
                final2, xytext=(12, -34), textcoords="offset points",
                ha="left", arrowprops=dict(arrowstyle="->", color=colors[1]))
    ax.scatter(0, 0, marker="+", s=70, color="black", label="pair centre of mass")
    ax.set_xlabel(r"approximate orbital-plane $x$ [$R_\odot$]")
    ax.set_ylabel(r"approximate orbital-plane $y$ [$R_\odot$]")
    ax.legend(loc="best", fontsize=9)
    ax.set_title(f"Actual hard.debug output samples: {len(times)};  "
                 f"last sampled separation = {final_separation:.6g} $R_\odot$")
    fig.suptitle(f"Merger {number:02d}: {event['id1']} + {event['id2']}   "
                 f"last two-body sample t = {times[-1]:.9g} Myr\n"
                 f"BSE event time = {event['time']:.9g} Myr   [{event['channel']}]",
                 fontsize=13)
    fig.text(0.5, 0.015,
             "Arrows connect real consecutive Hermite debug outputs; internal SDAR substeps are not printed.",
             ha="center", fontsize=9, color="0.35")
    fig.tight_layout(rect=(0, 0.045, 1, 0.93))
    filename = (f"debug_merger_{number:02d}_t{event['time']:.6f}_"
                f"id{event['id1']}_{event['id2']}.png")
    fig.savefig(output / filename, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return filename


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bse", default="data.bse")
    parser.add_argument("--logs", default="merger_debug_logs")
    parser.add_argument("--output-dir", default="merger_debug_orbits")
    parser.add_argument("--petar-tools", default="../tools")
    parser.add_argument("--dpi", type=int, default=180)
    args = parser.parse_args()
    sys.path.insert(0, str(Path(args.petar_tools).resolve()))
    import analysis as petar

    logs = Path(args.logs)
    output = Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    records = []
    for log in logs.glob("*.log"):
        if re.fullmatch(r"\d+\.\d+\.\d+\.\d+\.log", log.name):
            stem = log.name[:-4]
            sidecar = Path(f"dump_merger.{stem}.bse")
            id1, id2, _ = read_event_pair(sidecar)
            event = event_properties(Path(args.bse), id1, id2, petar)
            hard, n_particle, n_sd = load_hard_data(
                log, logs / f"{stem}.err", petar)
            track = extract_pair(hard, n_particle, id1, id2)
            if len(track[0]) == 0:
                raise RuntimeError(f"No simultaneous samples for {id1}, {id2}")
            xy1, xy2 = project_orbital_plane(*track[1:])
            records.append((event["time"], event, track[0], xy1, xy2,
                            stem, n_particle, n_sd))

    for number, (_, event, times, xy1, xy2, stem, n_particle, n_sd) in enumerate(
            sorted(records), 1):
        name = draw(event, times, xy1, xy2, number, output, stem, args.dpi)
        print(f"{name}: {len(times)} real output samples, "
              f"N={n_particle}, N_sd={n_sd}, dump={stem}")


if __name__ == "__main__":
    main()
