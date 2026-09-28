#!/usr/bin/env python3
"""Produce paper-style figures from actual completed/partial run data only."""
import argparse,csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
p=argparse.ArgumentParser();p.add_argument('root',type=Path);a=p.parse_args();root=a.root.resolve();out=root/'figures';out.mkdir(exist_ok=True)
repo=Path(__file__).resolve().parents[1];cases=json.loads((repo/'configs/paper/manifest.json').read_text())
def rows(path):
    with path.open() as stream:return [{k:float(v) for k,v in row.items()}for row in csv.DictReader(stream)]
def path_for(c):return root/c['suite']/c['name']
complete=[];audit=[]
for c in cases:
    d=path_for(c);done=(d/'summary.json').exists()
    audit.append(dict(name=c['name'],status='completed' if done else ('partial/failed' if (d/'history.csv').exists() else 'not-run')))
    if done and c['suite']=='convergence':complete.append((c,json.loads((d/'summary.json').read_text())))
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
for spatial in [False,True]:
    fig,axs=plt.subplots(1,3,figsize=(14,4.4),constrained_layout=True)
    for ax,key,label,order in zip(axs,['L2_spacetime_velocity','L2_spacetime_pressure','L2_spacetime_F'],['Velocity','Pressure','Deformation F'],[3,2,2]):
        grouping='ell' if spatial else 'j'
        for val in sorted(set(c[grouping]for c,s in complete)):
            data=sorted([(c['nominal_h']if spatial else c['dt'],s[key])for c,s in complete if c[grouping]==val])
            if data:ax.loglog(*zip(*data),marker='o',ms=3,label=f'{grouping}={val}')
        ax.set(title=label,xlabel='Nominal h'if spatial else 'Time step',ylabel='Space-time L2 error');ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.suptitle('Spatial convergence'if spatial else 'Temporal convergence')
    fig.savefig(out/('figure2-spatial.png'if spatial else 'figure1-temporal.png'),dpi=180);plt.close(fig)
# Slopes at fixed other parameter; no automatic claim that every point is asymptotic.
slope=[]
for ell in range(1,8):
    data=sorted([(c,s)for c,s in complete if c['ell']==ell],key=lambda cs:cs[0]['nominal_h'],reverse=True)
    for (c1,s1),(c2,s2) in zip(data,data[1:]):
        for field in ['velocity','pressure','F']:
            key='L2_spacetime_'+field
            slope.append(dict(ell=ell,j_from=c1['j'],j_to=c2['j'],field=field,observed_order=np.log(s1[key]/s2[key])/np.log(c1['nominal_h']/c2['nominal_h'])))
with (out/'spatial-orders.csv').open('w')as f:
    w=csv.DictWriter(f,fieldnames=['ell','j_from','j_to','field','observed_order']);w.writeheader();w.writerows(slope)
for quantity,num,title in [('energy',5,'Quadratic energy'),('log_energy',6,'Logarithmic energy diagnostic')]:
    fig,axs=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
    for ax,kind in zip(axs.flat,'abcd'):
        for c in cases:
            if c['suite']!='contraction'or c['case']!=kind:continue
            d=path_for(c)
            if not(d/'history.csv').exists():continue
            r=rows(d/'history.csv');partial=not(d/'summary.json').exists()
            x=np.array([v['time']for v in r]);y=np.array([v[quantity]for v in r]);ok=np.isfinite(y)
            ax.plot(x[ok],y[ok],label=f"Wi={c['wi']}"+(' (partial)'if partial else ''),ls='--'if partial else '-')
        ax.set(title=f'Case ({kind})',xlabel='Time',ylabel=title);ax.grid(alpha=.2)
        if ax.lines:ax.legend(fontsize=8)
    fig.savefig(out/f'figure{num}.png',dpi=180);plt.close(fig)
selections=[(.1,'a'),(1,'a'),(8,'a'),(8,'b'),(8,'c'),(8,'d')]
for stream,num in [(False,7),(True,8)]:
    fig,axs=plt.subplots(3,2,figsize=(12,9),constrained_layout=True)
    for ax,(wi,kind)in zip(axs.flat,selections):
        c=next(c for c in cases if c['suite']=='contraction'and c['wi']==wi and c['case']==kind);d=path_for(c)
        ax.set(title=f'Wi={wi}, case {kind}',xlim=(-1,1),ylim=(-1,1),aspect='equal')
        if not(d/'samples.csv').exists():ax.text(0,0,'No completed field output',ha='center');continue
        r=rows(d/'samples.csv');x=np.array([v['x']for v in r]);y=np.array([v['y']for v in r]);gx,gy=np.meshgrid(np.linspace(-1,1,300),np.linspace(-1,1,250));outside=(gx>0)&(abs(gy)>.5)
        if stream:
            u=griddata((x,y),[v['vx']for v in r],(gx,gy));v=griddata((x,y),[v['vy']for v in r],(gx,gy));u[outside]=np.nan;v[outside]=np.nan;ax.streamplot(gx,gy,u,v,density=1.8,linewidth=.55,color='k',arrowsize=.6)
        else:
            z=griddata((x,y),[v['stress']for v in r],(gx,gy));z[outside]=np.nan;im=ax.pcolormesh(gx,gy,z,shading='auto',cmap='magma');fig.colorbar(im,ax=ax,label='mu |F F^T - I|')
    fig.savefig(out/f'figure{num}.png',dpi=180);plt.close(fig)
report=['# Full numerical-experiment status','',f"Completed: {sum(v['status']=='completed'for v in audit)} / {len(audit)}",'', 'These figures use actual stored results. Partial/failed runs are not successful reproductions.','', '| Case | Status |','|---|---|']+[f"| {v['name']} | {v['status']} |"for v in audit]
report+=['','Mesh reconstruction: Gmsh h=0.05 plus five local bisection/conformity passes; exact original FEniCS mesh unavailable.','Explicit contraction coefficients imply alpha=0.9; paper text states 8/9.','No pointwise positivity guarantee is inferred from quadrature-point minima.']
(root/'REPORT.md').write_text('\n'.join(report)+'\n');(root/'audit.json').write_text(json.dumps(audit,indent=2)+'\n')
print('Wrote',out,'and',root/'REPORT.md')
