#!/usr/bin/env python3
"""Replay a supplied large-step dump in temporary directories; never edit inputs.
Usage: python3 test/slowdown_replay.py DUMP_DIRECTORY [--binary PATH]
Requires input.par.hard, input.par.bse, dump_large_step.* from the same build layout.
"""
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('directory', type=Path)
parser.add_argument('--binary', type=Path, default=Path(__file__).resolve().parents[1]/'build/petar.hard.debug')
parser.add_argument('--output', type=Path)
args = parser.parse_args()
dump = next(p for p in sorted(args.directory.glob('dump_large_step.*')) if re.fullmatch(r'dump_large_step\.\d+\.\d+\.\d+\.\d+',p.name))
results=[]
for label, extra in [('default', []), ('tighter', ['-a','1e-9'])]:
    with tempfile.TemporaryDirectory(prefix='petar-slowdown-') as tmp:
        path=Path(tmp)
        for source in [dump,args.directory/'input.par.hard',args.directory/'input.par.bse']:
            shutil.copy2(source,path/source.name)
        random_state=args.directory/'bse.rand.par'
        if random_state.exists():shutil.copy2(random_state,path/random_state.name)
        env=os.environ.copy()
        env['OMP_NUM_THREADS']='1'
        env['PYTHONPATH']=str(Path(__file__).resolve().parents[1]/'bhmerger/precession-master')+os.pathsep+env.get('PYTHONPATH','')
        cmd=[str(args.binary.resolve()),'-e','1e-3','-p','input.par.hard','-b','input.par.bse','-s','1000000',*extra,dump.name]
        result=subprocess.run(cmd,cwd=path,env=env,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,timeout=30)
        log=result.stderr
        if args.output:
            args.output.mkdir(parents=True,exist_ok=True)
            (args.output/f'{label}.stderr').write_text(log)
        assert result.returncode==0, (label,result.returncode,log[-1000:])
        assert 'AR_ADAPTIVE_MULTIPLE' in log
        assert not any(x in log for x in ['AR_TIME_STALL','Negative_step:', 'Hard energy significant!', 'Warning: step count'])
        final=next(line for line in reversed(log.splitlines()) if line.startswith('Hard Energy:'))
        fields={key:float(value) for key,value in re.findall(r'(\w+):\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)',final)}
        e_initial=fields['Ekin']+fields['Epot']-fields['dE_change']-fields['dE']
        residual=abs(fields['dE']/e_initial)
        assert residual<1e-5,(label,residual)
        results.append(dict(case=label,relative_energy_residual=residual,AR_steps=int(fields['ARC_substep_sum']),recovery_events=log.count('AR_SLOWDOWN_RECOVERY ')))
print(json.dumps(results,indent=2))
