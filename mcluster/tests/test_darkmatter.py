#!/usr/bin/env python3
"""Compile opt-in/off variants and test numerical physics, formats and isolation.
Optionally set PETAR_DM to a DARKMATTER, interrupt=off, external=off executable
to run a real PeTar input/restart test. No third-party Python packages required.
"""
import math
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def run(args, cwd, ok=True):
    p = subprocess.run([str(a) for a in args], cwd=cwd, text=True,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT, errors="replace")
    if (p.returncode == 0) != ok:
        raise AssertionError(f"{args}\nexit={p.returncode}\n{p.stdout}")
    return p.stdout


def read(path):
    with open(path) as f:
        header = f.readline().split()
        rows = [[float(x) for x in line.split()] for line in f]
    assert len(rows) == int(header[1])
    assert all(math.isfinite(v) for row in rows for v in row)
    return header, rows


def main():
    with tempfile.TemporaryDirectory(prefix="mcluster-nfw-") as directory:
        d = Path(directory)
        cc = shlex.split(os.environ.get("CC", "gcc"))
        common = cc + ["-O2", "-DNOOMP"]
        for name, flags in [("plain", []), ("dm", ["-DDARKMATTER"])]:
            run(common + flags + [ROOT / "main.c", "-lm", "-o", d / name], d)
        run(common + ["-DDARKMATTER", ROOT / "tests/darkmatter_numeric.c",
                      "-lm", "-o", d / "dm_numeric_test"], d)
        print(run([d / "dm_numeric_test"], d))
        # Same seed, no new options: legacy results must be byte-identical.
        # Optional pre-change executable adds a true before/after comparison.
        variants = [d / "plain", d / "dm"]
        if os.environ.get("MCLUSTER_BASELINE"):
            variants.append(Path(os.environ["MCLUSTER_BASELINE"]).resolve())
        for profile, extra in [(0, []), (1, ["-W", "5"]), (0, ["-B", "8"])]:
            reference = None
            for j, executable in enumerate(variants):
                prefix = f"old{profile}_{len(extra)}_{j}"
                run([executable, "-N", "128", "-f", "0", "-m", "1", "-P",
                     str(profile), "-t", "0", "-s", "17", "-o", prefix]+extra, d)
                payload = (d / (prefix+".txt")).read_bytes()
                if reference is None:
                    reference = payload
                assert payload == reference, "legacy output changed"
        base = [d / "dm", "--dm-nfw", "--dm-number", "4096",
                "--dm-mass", "10000", "--dm-rs", "10", "--dm-rt", "100",
                "--dm-softening", "1", "-N", "1024", "-f", "0", "-m", "1",
                "-R", "1", "-t", "0", "-s", "42"]
        print(run(base + ["-o", "model"], d))
        hs, stars = read(d / "model.star")
        hd, dark = read(d / "model.dm")
        assert hs == ["0", "1024", "0"] and hd == ["0", "4096", "0", "1"]
        assert all(len(row) == 20 for row in stars)
        assert all(len(row) == 13 for row in dark)
        assert [int(row[9]) for row in stars] == list(range(1,1025))
        assert [int(row[7]) for row in dark] == list(range(1,4097))
        assert abs(sum(row[0] for row in stars)-1024) < 1e-9
        assert abs(sum(row[0] for row in dark)-10000) < 1e-8
        for j in range(1,7):
            assert abs(sum(row[0]*row[j] for row in stars+dark)/11024) < 1e-11
        # Reproducible output, and correct BSE/BH schema (34 columns).
        run(base+["-o","repeat"], d)
        for suffix in [".star",".dm",".dm.profile"]:
            assert (d/("model"+suffix)).read_bytes() == (d/("repeat"+suffix)).read_bytes()
        run(base+["--petar-interrupt","bse","--blackhole","-o","bse"], d)
        _, bse = read(d/"bse.star")
        assert all(len(row)==34 and row[12]==14 and row[16]==row[0] for row in bse)
        assert [int(row[23]) for row in bse] == list(range(1,1025))
        # Reject invalid options, unsupported physics and nonpositive DFs.
        for extra in [["-P","1"],["-B","1"],["-S","0.5"],["-D","2"],
                      ["-Q","0.4"],["-t","3"],["-e","1"],["-u","0"],
                      ["-C","5"],["--dm-number","1.5"],["--dm-mass","nan"],
                      ["--dm-softening","0"],["--dm-grid","100"]]:
            run(base+extra+["-o","bad"], d, ok=False)
            assert not (d/"bad.star").exists()
        run([d/"dm","--dm-number","10"],d,ok=False)
        result=run(base+["--dm-mass","1000000","--dm-rs","1","--dm-rt","10",
                         "--dm-softening","0.0001","-o","negative"],d,ok=False)
        assert "negative stellar DF" in result
        assert not (d/"negative.star").exists()
        # Large-scale density CDF sampling check, accounting for the common shift.
        info = dict(line.strip().split("=",1) for line in (d/"model.dm.info").read_text().splitlines())
        centre = [float(x) for x in info["subtracted_common_position"].split()]
        radii = sorted(math.sqrt(sum((row[j+1]+centre[j])**2 for j in range(3))) for row in dark)
        with open(d/"model.dm.profile") as f:
            next(f)
            table = [[float(x) for x in line.split()] for line in f]
        for quantile in [0.1,0.25,0.5,0.75,0.9]:
            r = min(table,key=lambda row:abs(row[8]-quantile))[0]
            measured = sum(x <= r for x in radii)/len(radii)
            assert abs(measured-quantile)<0.03, (measured,quantile)
        petar = os.environ.get("PETAR_DM")
        if petar:
            source = ROOT.parent/"PeTar_new_version"
            if source.is_dir():
                cxx = shlex.split(os.environ.get("CXX", "g++"))
                run(cxx+["-O2","-std=c++17","-DDARKMATTER",
                         "-DSTELLAR_EVOLUTION","-DBSE_BASE","-DBSEBBF",
                         "-I"+str(source/"src"),"-I"+str(source/"bse-interface"),
                         "-I"+str(ROOT.parent/"FDPS/src"),
                         "-I"+str(ROOT.parent/"SDAR/src"),
                         ROOT/"tests/petar_read.cxx","-o",d/"reader"],d)
                run([d/"reader","bse.star","bse.dm"],d)
                print("PeTar original BSE/DM readers passed")
            cmd = [Path(petar).resolve(), "-u","1","-i","1","-s","0.0009765625",
                   "-o","0.0009765625","--dm-input","model.dm"]
            run(cmd+["-t","0.0009765625","model.star"],d)
            sh,_=read(d/"data.1")
            dh,_=read(d/"dmdata.1")
            assert sh[:3] == dh[:3] or (sh[0]==dh[0] and float(sh[2])==float(dh[2]))
            assert float(dh[3])==1
            restart=[Path(petar).resolve(),"-u","1","-i","1","-s","0.0009765625",
                     "-o","0.0009765625","-t","0.001953125",
                     "--dm-input","dmdata.1","data.1"]
            run(restart,d)
            dh2,_=read(d/"dmdata.2")
            assert float(dh2[3])==1 and float(dh2[2])==0.001953125
            print("PeTar input and paired restart passed")
        print("legacy isolation, formats, sampling and rejection tests passed")


if __name__ == "__main__":
    main()
