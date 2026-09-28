#!/usr/bin/env python3
"""Measure unchanged discrete problems before selecting CPU/memory concurrency."""
import argparse,json,os,re,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('output',type=Path);a=p.parse_args()
repo=Path(__file__).resolve().parents[1];root=a.output.resolve();root.mkdir(parents=True,exist_ok=True)
info=Path('/proc/meminfo').read_text();available=int(re.search(r'MemAvailable:\s+(\d+)',info)[1])*1024
rows=subprocess.check_output(['lscpu','-p=CORE,SOCKET'],text=True).splitlines()
cores=len(set(r for r in rows if not r.startswith('#')))
budget=int(available*.75);records=[];settings={}
for suite,name in [('convergence','mms-j7-ell1'),('contraction','contraction-wi1-a')]:
    candidates=[]
    for threads in [1,2,4]:
        if threads>cores:continue
        config=(repo/'configs/paper'/f'{name}.prm').read_text()
        config=re.sub(r'set end_time = .*','set end_time = 0.02',config)
        path=root/f'{suite}-threads{threads}.prm';path.write_text(config)
        dest=root/f'{suite}-threads{threads}';stats=root/f'{suite}-threads{threads}.time'
        env=os.environ.copy();env.update({k:str(threads) for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']})
        env['DEAL_II_NUM_THREADS']='1'
        start=time.monotonic()
        with (root/f'{suite}-threads{threads}.log').open('w') as log:
            r=subprocess.run(['/usr/bin/time','-f','%M', '-o',str(stats),str(repo/'scripts/run.sh'),str(path),str(dest)],cwd=repo,env=env,stdout=log,stderr=subprocess.STDOUT)
        if r.returncode:raise RuntimeError(f'Calibration failed: {suite}, {threads} threads; inspect log')
        peak=int(stats.read_text().strip())*1024;elapsed=time.monotonic()-start
        jobs=min(cores//threads,int(budget/(peak*1.35)))
        if jobs<1:raise RuntimeError(f'{suite} exceeds safe memory budget; do not launch full queue')
        record=dict(suite=suite,threads=threads,seconds=elapsed,peak_rss_bytes=peak,jobs=jobs,estimated_cases_per_second=jobs/elapsed)
        candidates.append(record);records.append(record)
        (root/'measurements.json').write_text(json.dumps(dict(physical_cores=cores,memory_budget_bytes=budget,measurements=records),indent=2)+'\n')
    # Independent cases benefit from throughput, not just single-case speed.
    best=max(candidates,key=lambda x:x['estimated_cases_per_second'])
    settings[suite]=best
(root/'settings.json').write_text(json.dumps(settings,indent=2)+'\n')
(root/'settings.sh').write_text(''.join(f'export {suite.upper()}_JOBS={v["jobs"]}\nexport {suite.upper()}_THREADS={v["threads"]}\n' for suite,v in settings.items()))
print(json.dumps(settings,indent=2))
