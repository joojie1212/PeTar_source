from pathlib import Path
import subprocess, numpy as np
root=Path('/tmp/astra-mcluster-cases');root.mkdir(exist_ok=True)
exe=str(Path(__file__).resolve().parents[1]/'mcluster_sse')
cases={'mono':['-f','2','-m','2'], 'power':['-f','2','-m','1','-m','2','-a','-2.3'], 'optimal':['-f','3'], 'l3':['-f','4'], 'long_'+'x'*90:['-f','0','-m','1']}
for name,extra in cases.items():
 with (root/(name+'.log')).open('w') as out:
  p=subprocess.run([exe,('-M' if name=='optimal' else '-N'),'20','-B','0','-s','42','-C','3','-o',name,*extra],cwd=root,stdout=out,stderr=subprocess.STDOUT,timeout=30)
 assert p.returncode==0,(name,p.returncode)
 data=np.loadtxt(root/(name+'.txt'))
 assert np.all(np.isfinite(data)),name
 if name=='mono':assert np.all(data[:,0]==2.)
 if name=='power':assert np.all(data[:,11]>0.) and np.all(data[:,12]>0.)
 print(name[:12],data.shape,'passed')

for code in [4,5]:
 name='format'+str(code)+'_'+'x'*90
 with (root/(name+'.log')).open('w') as out:
  p=subprocess.run([exe,'-N','20','-f','0','-m','1','-s','42','-C',str(code),'-o',name],cwd=root,stdout=out,stderr=subprocess.STDOUT,timeout=30)
 assert p.returncode==0,(code,p.returncode)
 assert (root/(name+'.input')).is_file()
 assert (root/(name+('.fort.10' if code==4 else '.dat.10'))).is_file()
 print('long filename format',code,'passed')
