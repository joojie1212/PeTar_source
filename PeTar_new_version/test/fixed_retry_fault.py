#!/usr/bin/env python3
"""Inject an error after fixed-step warmup using a temporary header, not production code."""
from pathlib import Path
import subprocess
import tempfile
root=Path(__file__).resolve().parents[2]
source=(root/'SDAR/src/AR/symplectic_integrator.h').read_text()
needle='                Float integration_error_rel_abs = abs(H-H_bk);'
assert source.count(needle)==1
with tempfile.TemporaryDirectory(prefix='petar-fixed-retry-') as tmp:
    tmp=Path(tmp);(tmp/'AR').mkdir()
    for case, injection in [('energy','integration_error_rel_abs=100*energy_error_rel_max'),
                             ('nonfinite','integration_error_rel_abs=std::numeric_limits<Float>::quiet_NaN()')]:
        (tmp/'AR/symplectic_integrator.h').write_text(source.replace(needle,needle+f'\n                if (step_count==8) {injection};'))
        binary=tmp/case
        cmd=['g++','-std=c++17','-O1','-fopenmp','-DEXPECT_FIXED_RETRY','-DAR_TTL',
             '-DAR_SLOWDOWN_TREE','-DAR_SLOWDOWN_TIMESCALE','-DHARD_CHECK_ENERGY',f'-I{tmp}']
        cmd += [f'-I{root/p}' for p in ['SDAR/src','FDPS/src','PeTar_new_version/src',
                                      'PeTar_new_version/bse-interface','PeTar_new_version/bhmerger']]
        cmd += [str(root/'PeTar_new_version/test/bugfix/fixed_retry.cxx'),'-o',str(binary)]
        subprocess.run(cmd,check=True)
        subprocess.run([str(binary)],check=True,timeout=30)
        print(case,'PASS',flush=True)
