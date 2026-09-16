#!/usr/bin/env python3
"""Plot an approximate orbital-plane approach for every merger in data.bse.

PeTar's BSE event log stores pre-merger orbital elements, not an integration
step trajectory.  The curves drawn here are therefore Keplerian
reconstructions from the last logged a, e, masses and stellar radii.  Each
curve ends at first surface contact, r = R1 + R2.
"""

import argparse
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np


def load_events(path, petar):
    coalescence, dynamic = [], []
    with path.open() as stream:
        for line in stream:
            fields = line.split()
            if not fields:
                continue
            if fields[0] == "Coalescence":
                coalescence.append([float(x) for x in fields[1:]])
            elif fields[0] == "Dynamic_merge:":
                dynamic.append([float(x) for x in fields[1:]])

    events = []
    if coalescence:
        data = petar.BSETypeChange(np.asarray(coalescence))
        for i in range(data.size):
            events.append(dict(
                channel="BSE_Coalescence", time=float(data.init.time[i]),
                id1=int(data.id1[i]), id2=int(data.id2[i]),
                m1=float(data.init.m1[i]), m2=float(data.init.m2[i]),
                r1=float(data.init.rad1[i]), r2=float(data.init.rad2[i]),
                semi=float(data.init.semi[i]), ecc=float(data.init.ecc[i])))
    if dynamic:
        data = petar.BSEDynamicMerge(np.asarray(dynamic))
        for i in range(data.size):
            events.append(dict(
                channel="Dynamic_merge", time=float(data.init.p1.time[i]),
                id1=int(data.id1[i]), id2=int(data.id2[i]),
                m1=float(data.init.p1.mass[i]), m2=float(data.init.p2.mass[i]),
                r1=float(data.init.p1.rad[i]), r2=float(data.init.p2.rad[i]),
                semi=float(data.semi[i]), ecc=float(data.ecc[i])))
    return sorted(events, key=lambda x: x["time"])


def reconstructed_approach(event, n=600, extent_in_radii=18.0):
    """Return barycentric tracks from an inbound point to first contact."""
    a, e = event["semi"], event["ecc"]
    contact = event["r1"] + event["r2"]
    p = a * (1.0 - e * e)
    peri = a * (1.0 - e)
    apo = a * (1.0 + e)

    # r = p / (1 + e cos f).  Negative f is the inbound branch.
    def anomaly_at_radius(radius):
        cosine = np.clip((p / radius - 1.0) / e, -1.0, 1.0)
        return -np.arccos(cosine)

    # Every recorded event here reaches contact according to its osculating
    # orbit.  Guard the general case by ending at peri-centre when it does not.
    end_radius = contact if peri <= contact <= apo else peri
    f_end = anomaly_at_radius(end_radius) if end_radius > peri else 0.0
    start_radius = min(apo, max(extent_in_radii * contact,
                                1.15 * end_radius))
    f_start = anomaly_at_radius(start_radius)
    if f_start >= f_end:
        f_start = max(-np.pi + 1e-4, f_end - 0.8)
    f = np.linspace(f_start, f_end, n)
    radius = p / (1.0 + e * np.cos(f))
    relative = np.column_stack((radius * np.cos(f), radius * np.sin(f)))
    total_mass = event["m1"] + event["m2"]
    track1 = -event["m2"] / total_mass * relative
    track2 = event["m1"] / total_mass * relative
    return track1, track2, peri, contact


def draw_event(event, number, output_dir, dpi):
    track1, track2, peri, contact = reconstructed_approach(event)
    colors = ("#277da1", "#f3722c")
    fig, ax = plt.subplots(figsize=(8.4, 7.0))
    ax.set_aspect("equal", adjustable="box")
    ax.axhline(0, color="0.88", lw=0.7, zorder=0)
    ax.axvline(0, color="0.88", lw=0.7, zorder=0)
    ax.set_xlabel(r"orbital-plane $x$ [$R_\odot$]")
    ax.set_ylabel(r"orbital-plane $y$ [$R_\odot$]")

    ax.plot(track1[:, 0], track1[:, 1], color=colors[0], lw=2,
            label=f"ID {event['id1']}")
    ax.plot(track2[:, 0], track2[:, 1], color=colors[1], lw=2,
            label=f"ID {event['id2']}")
    step = max(2, len(track1) // 10)
    for track, color in ((track1, colors[0]), (track2, colors[1])):
        for j in range(step, len(track), step):
            ax.annotate("", xy=track[j], xytext=track[j-step],
                        arrowprops=dict(arrowstyle="->", color=color,
                                        lw=1.0, alpha=0.75))
    ax.scatter(0, 0, marker="+", s=70, color="black", label="centre of mass")
    final1, final2 = track1[-1], track2[-1]
    ax.add_patch(Circle(final1, event["r1"], facecolor=colors[0],
                        edgecolor="black", lw=1.2, alpha=0.78, zorder=5))
    ax.add_patch(Circle(final2, event["r2"], facecolor=colors[1],
                        edgecolor="black", lw=1.2, alpha=0.78, zorder=5))
    ax.plot([final1[0], final2[0]], [final1[1], final2[1]],
            color="0.25", ls="--", lw=1, zorder=4)
    ax.annotate(f"ID {event['id1']}\n$R_1={event['r1']:.4g}\,R_\odot$",
                xy=final1, xytext=(-12, 24), textcoords="offset points",
                ha="right", fontsize=10,
                arrowprops=dict(arrowstyle="->", color=colors[0]))
    ax.annotate(f"ID {event['id2']}\n$R_2={event['r2']:.4g}\,R_\odot$",
                xy=final2, xytext=(12, -34), textcoords="offset points",
                ha="left", fontsize=10,
                arrowprops=dict(arrowstyle="->", color=colors[1]))
    ax.set_title("Reconstructed inbound motion and first surface contact")
    ax.legend(loc="best", fontsize=9)

    fig.suptitle(
        f"Merger {number:02d}: {event['id1']} + {event['id2']}   "
        f"t = {event['time']:.6g} Myr   [{event['channel']}]\n"
        f"pre-merger orbit: a = {event['semi']:.6g} $R_\odot$, "
        f"e = {event['ecc']:.7f},  r_peri = {peri:.5g} $R_\odot$",
        fontsize=12)
    fig.text(0.5, 0.015,
             "Approximation reconstructed from the last logged osculating orbit; not raw integration samples.",
             ha="center", fontsize=9, color="0.35")
    fig.tight_layout(rect=(0, 0.045, 1, 0.91))
    name = (f"merger_{number:02d}_t{event['time']:.6f}_"
            f"id{event['id1']}_{event['id2']}.png")
    fig.savefig(output_dir / name, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return name, peri, contact


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bse", default="data.bse")
    parser.add_argument("--output-dir", default="merger_orbits")
    parser.add_argument("--petar-tools", default="../tools")
    parser.add_argument("--dpi", type=int, default=180)
    args = parser.parse_args()
    sys.path.insert(0, str(Path(args.petar_tools).resolve()))
    import analysis as petar

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    events = load_events(Path(args.bse), petar)
    rows = []
    for number, event in enumerate(events, 1):
        filename, peri, contact = draw_event(event, number, output_dir,
                                             args.dpi)
        rows.append([number, event["channel"], event["time"], event["id1"],
                     event["id2"], event["m1"], event["m2"], event["r1"],
                     event["r2"], event["semi"], event["ecc"], peri,
                     contact, filename])
        print(f"Wrote {output_dir / filename}")
    with (output_dir / "merger_orbits.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["number", "channel", "time_Myr", "id1", "id2",
                         "mass1_Msun", "mass2_Msun", "radius1_Rsun",
                         "radius2_Rsun", "semi_Rsun", "ecc", "peri_Rsun",
                         "contact_separation_Rsun", "figure"])
        writer.writerows(rows)
    print(f"Generated {len(events)} merger figures")


if __name__ == "__main__":
    main()
