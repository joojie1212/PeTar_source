"""No simulations: exercise resource flags, queue refill, failure handling and scaling."""
import csv
import importlib.util
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('sweep', HERE/'sweep.py')
sweep = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sweep)


class SweepTests(unittest.TestCase):
    def test_scaling(self):
        rows = list(sweep.models())
        self.assertEqual(len(rows), 15)
        for ref in sweep.RH_REFS:
            group = [r for r in rows if r['rh_ref_pc'] == ref]
            self.assertEqual(len(group), 5)
            for row in group:
                self.assertTrue(math.isclose(row['rho0_Msun_pc3'], group[0]['rho0_Msun_pc3'], rel_tol=1e-12))
                self.assertAlmostEqual(row['end_Myr']/row['trh_Myr'], 300)

    def test_binary_restart(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/'input').write_text('ASCII initial model')
            (root/'petar').write_text("""#!/usr/bin/env python3
import sys
from pathlib import Path
args=sys.argv[1:]
mode=args[args.index('-i')+1]
if not Path('cycle1').exists():
    assert mode=='2' and args[-1]=='input', args
    Path('cycle1').touch()
    Path('data.1').write_bytes(b'\\x00binary snapshot')
    Path('input.par').touch()
    Path('data.hard_reopt').write_text('exit_code=42\\nsnapshot=data.1\\n')
    sys.exit(42)
assert mode=='0' and args[-1]=='data.1', args
assert '-p' in args and '-s' in args
Path('restarted').touch()
""")
            # Supervisor requires a nonempty parameter file.
            text = (root/'petar').read_text().replace("Path('input.par').touch()", "Path('input.par').write_text('parameters')")
            (root/'petar').write_text(text)
            (root/'petar.find.dt').write_text("""#!/usr/bin/env python3
import sys
from pathlib import Path
args=sys.argv[1:]
assert args[args.index('-i')+1]=='0' and args[-1]=='data.1', args
assert 'mpirun --mca plm isolated --mca ras ^slurm --bind-to none -np 1' in args[args.index('-r')+1]
Path('optimized').touch()
print('Best performance choice: tree step: 0.001')
""")
            (root/'mpirun').write_text("""#!/usr/bin/env python3
import os, sys
args=sys.argv[1:]
assert args[:10] == ['--mca','plm','isolated','--mca','ras','^slurm','--bind-to','none','-np','1'], args
os.execvp(args[10], args[10:])
""")
            for name in ('petar', 'petar.find.dt', 'mpirun'):
                (root/name).chmod(0o755)
            env = dict(os.environ, PATH=str(root)+os.pathsep+os.environ['PATH'],
                       WORK_DIR=tmp, INPUT_SNAPSHOT='input')
            result = subprocess.run(['bash', str(HERE/'supervisor.sh')], env=env,
                                    capture_output=True, text=True, timeout=10)
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
            self.assertTrue((root/'optimized').exists())
            self.assertTrue((root/'restarted').exists())

    def test_queue_and_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fake = root/'srun'
            fake.write_text('''#!/usr/bin/env python3
import fcntl, json, os, sys, time
from pathlib import Path
assert all(x in sys.argv for x in ['--exclusive','--exact','--ntasks=1','--cpus-per-task=4','--mem=4G','--cpu-bind=cores','--mpi=none'])
assert os.environ['OMP_NUM_THREADS'] == '4'
p=Path(os.environ['MOCK_STATE'])
def update(delta):
    with p.open('a+') as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        f.seek(0); text=f.read(); state=json.loads(text) if text else {'active':0,'peak':0,'count':0}
        state['active']+=delta
        state['peak']=max(state['peak'],state['active'])
        state['count']+=int(delta>0)
        f.seek(0); f.truncate(); json.dump(state,f); f.flush()
update(1)
time.sleep(1.4 if os.environ['TARGET_N']=='500' else 0.2)
update(-1)
sys.exit(7 if os.environ['TARGET_N']=='1000' and os.environ['REFERENCE_HALF_MASS_RADIUS_PC']=='0.013047672' else 0)
''')
            fake.chmod(0o755)
            env = dict(os.environ, PATH=str(root)+os.pathsep+os.environ['PATH'],
                       SLURM_JOB_ID='mock', ESCAPE_RADIUS_RH='200', OUTPUT_ROOT=str(root/'results'), MOCK_STATE=str(root/'state'))
            result = subprocess.run(['python3', str(HERE/'sweep.py'), '--run'], env=env,
                                    capture_output=True, text=True, timeout=40)
            self.assertEqual(result.returncode, 1, result.stdout+result.stderr)
            import json
            state = json.loads((root/'state').read_text())
            self.assertEqual(state, {'active':0,'peak':4,'count':15})
            batch = next((root/'results').iterdir())
            with (batch/'status.tsv').open() as f:
                rows = list(csv.DictReader(f, delimiter='\t'))
            self.assertEqual(sum(r['status']=='completed' for r in rows),14)
            errors = [r for r in rows if r['status']=='failed']
            self.assertEqual(len(errors),1)
            self.assertEqual(errors[0]['exit_code'],'7')
            configs=list(batch.glob('N*/config.txt'))
            self.assertEqual(len(configs),15)
            self.assertTrue(all('Escape radius (200 initial rh)' in p.read_text() for p in configs))


if __name__ == '__main__':
    unittest.main()
