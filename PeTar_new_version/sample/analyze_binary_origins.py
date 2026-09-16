#!/usr/bin/env python3
"""Infer the origins of binaries in PeTar snapshots.

A tidal capture is counted only when a Tide record changes an unbound orbit
(semi < 0) into a bound orbit (semi > 0, ecc < 1).  Dynamical subclasses are
inferred at the snapshot cadence: a new pair containing a member of a
different pair in the preceding snapshot is labelled an exchange candidate.
"""

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

import numpy as np


def ordered_pair(a, b):
    return min(int(a), int(b)), max(int(a), int(b))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix", default="data", help="PeTar output prefix")
    parser.add_argument("--G", type=float, default=0.00449830997959438)
    parser.add_argument("--rmax", type=float, default=0.1,
                        help="maximum separation used by findPair [pc]")
    parser.add_argument("--petar-tools", default="../tools")
    args = parser.parse_args()

    sys.path.insert(0, str(Path(args.petar_tools).resolve()))
    import analysis as petar

    prefix = Path(args.prefix)
    tide_data = np.loadtxt(str(prefix) + ".bse.tide")
    tide = petar.Tide(tide_data, interrupt_mode="bse")
    tide_pairs = {ordered_pair(a, b) for a, b in zip(tide.id1, tide.id2)}
    capture_mask = ((tide.semi0 < 0) & (tide.bin.semi > 0)
                    & (tide.bin.ecc < 1))
    capture_pairs = {
        ordered_pair(a, b)
        for a, b in zip(tide.id1[capture_mask], tide.id2[capture_mask])
    }

    snap_list = Path(str(prefix) + ".snap.lst")
    paths = [Path(line.strip()) for line in snap_list.read_text().splitlines()
             if line.strip()]
    first = {}
    seen = set()
    previous = set()
    final_candidates = set()
    final_hard = set()
    final_time = np.nan

    for index, path in enumerate(paths):
        header = petar.PeTarDataHeader(str(path), interrupt_mode="bse")
        particles = petar.Particle(interrupt_mode="bse")
        particles.loadtxt(path, skiprows=1)
        _, _, binaries = petar.findPair(
            particles, args.G, args.rmax, use_kdtree=True, simple_binary=True)

        candidates = set()
        hard = set()
        for id1, id2, state1, state2 in zip(
                binaries.p1.id, binaries.p2.id,
                binaries.p1.binary_state, binaries.p2.binary_state):
            pair = ordered_pair(id1, id2)
            candidates.add(pair)
            # Both stars naming each other excludes incidental wide pairs.
            if ((int(state1) >> 4) == int(id2)
                    and (int(state2) >> 4) == int(id1)):
                hard.add(pair)

        previous_members = {member for pair in previous for member in pair}
        for pair in hard - seen:
            if pair in capture_pairs:
                origin = "tidal_capture"
            elif pair[0] in previous_members or pair[1] in previous_members:
                origin = "exchange_candidate"
            else:
                origin = "three_body_candidate"
            first[pair] = (float(header.time), origin)
        seen.update(hard)
        previous = hard
        final_candidates, final_hard = candidates, hard
        final_time = float(header.time)

    output = Path(str(prefix) + ".binary_origins.csv")
    with output.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["id1", "id2", "first_snapshot_time_Myr",
                         "inferred_origin", "ever_tide_affected",
                         "present_in_final_snapshot"])
        for pair in sorted(seen):
            writer.writerow([*pair, *first[pair], pair in tide_pairs,
                             pair in final_hard])

    print(f"Snapshots: {len(paths)}; final time: {final_time:g} Myr")
    print(f"Unique tracked hard pairs: {len(seen)} {dict(Counter(x[1] for x in first.values()))}")
    print(f"Unique tide-affected pairs: {len(tide_pairs)}")
    print(f"Strict tidal captures: {len(capture_pairs)}")
    print(f"Final candidates: {len(final_candidates)}; tracked hard binaries: {len(final_hard)}")
    print(f"Final origins: {dict(Counter(first[x][1] for x in final_hard))}")
    print(f"Final tide-affected (not necessarily captured): {len(final_hard & tide_pairs)}")
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
