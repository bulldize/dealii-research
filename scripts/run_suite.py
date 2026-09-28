#!/usr/bin/env python3
"""Run independent paper cases; preserve failed runs and resume completed cases safely."""
import argparse,concurrent.futures,hashlib,json,os,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--suite',choices=['convergence','contraction','all'],required=True);p.add_argument('--jobs',type=int,default=1);p.add_argument('--output',type=Path,required=True);p.add_argument('--case',action='append',default=[]);a=p.parse_args()
if a.jobs<1:raise ValueError('jobs must be positive')
repo=Path(__file__).resolve().parents[1];root=a.output.resolve();root.mkdir(parents=True,exist_ok=True)
cases=json.loads((repo/'configs/paper/manifest.json').read_text());cases=[c for c in cases if (a.suite=='all' or c['suite']==a.suite) and (not a.case or c['name'] in a.case)]
if not cases:raise ValueError('No matching cases')
files=[repo/'CMakeLists.txt']+[f for d in ['include','src','apps'] for f in sorted((repo/d).rglob('*')) if f.is_file()]
source_digest=hashlib.sha256(b''.join(str(f.relative_to(repo)).encode()+f.read_bytes() for f in files)).hexdigest()
executable=repo/'build/giesekus';binary_digest=hashlib.sha256(executable.read_bytes()).hexdigest()
def work(case):
    dest=root/case['name'];config=repo/case['config']
    signature=dict(source_sha256=source_digest,binary_sha256=binary_digest,config_sha256=hashlib.sha256(config.read_bytes()).hexdigest())
    if case['suite']=='contraction': signature['mesh_sha256']=hashlib.sha256((repo/'meshes/contraction-paper.msh').read_bytes()).hexdigest()
    if dest.exists():
        if (dest/'summary.json').exists() and (dest/'provenance.json').exists():
            old=json.loads((dest/'provenance.json').read_text())
            if old['signature']==signature and json.loads((dest/'summary.json').read_text()).get('completed'):return dict(name=case['name'],status='already-completed')
        raise RuntimeError(f'{dest} is incomplete or has different inputs: preserve it and use a new output root')
    dest.mkdir();(dest/'used-config.prm').write_bytes(config.read_bytes())
    (dest/'provenance.json').write_text(json.dumps(dict(case=case,signature=signature,execution_environment={k:os.environ.get(k) for k in ["OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","DEAL_II_NUM_THREADS"]}),indent=2)+'\n')
    start=time.time();print('START',case['name'],flush=True)
    try:
        with (dest/'run.log').open('w') as log:
            result=subprocess.run([str(repo/'scripts/run.sh'),str(config),str(dest)],cwd=repo,stdout=log,stderr=subprocess.STDOUT)
        state=dict(name=case['name'],returncode=result.returncode,status='completed' if result.returncode==0 and (dest/'summary.json').exists() else 'failed',seconds=time.time()-start)
    except Exception as exc:state=dict(name=case['name'],status='failed',error=str(exc),seconds=time.time()-start)
    (dest/'status.json').write_text(json.dumps(state,indent=2)+'\n');print('END',state,flush=True);return state
states=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as pool:
    futures={pool.submit(work,c):c for c in cases}
    for f in concurrent.futures.as_completed(futures):
        try:states.append(f.result())
        except Exception as exc:states.append(dict(name=futures[f]['name'],status='not-run',error=str(exc)))
        temp=root/'suite-status.json.tmp';temp.write_text(json.dumps(states,indent=2)+'\n');temp.replace(root/'suite-status.json')
if any(s['status'] not in ['completed','already-completed'] for s in states):raise SystemExit(1)
