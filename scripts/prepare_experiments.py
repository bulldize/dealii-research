#!/usr/bin/env python3
"""Write the paper's complete numerical-experiment matrix without running it."""
import json
from pathlib import Path
repo=Path(__file__).resolve().parents[1]
base=(repo/'configs/paper');base.mkdir(parents=True,exist_ok=True)
cases=[]
def save(name,settings,metadata):
    path=base/(name+'.prm')
    path.write_text('# Generated from arXiv:2512.22831v1, Section 5\n'+''.join(f'set {k} = {v}\n' for k,v in settings.items()))
    cases.append(dict(name=name,config=str(path.relative_to(repo)),**metadata))
for j in range(2,8):
    for ell in range(1,8):
        # Literal vertex count from Section 5.2; distinguish nominal h from cell diameter.
        n=2**j-1;dt=.02*2**(-ell)
        save(f'mms-j{j}-ell{ell}',dict(problem='manufactured',subdivisions=n,dt=dt,end_time=.1,rho=1,nu=1,mu=1,**{'lambda':1,'newton_tolerance':1e-12,'max_newton':15}),dict(suite='convergence',j=j,ell=ell,nominal_h=2**(-j),cell_diameter=2**.5/n,dt=dt,steps=round(.1/dt)))
for wi,lam in [(0,0),(.1,.125),(1,1.25),(3,3.75),(5,6.25),(8,10)]:
    for kind,dt,diff in [('a',.01,.01),('b',.0025,.0025),('c',.01,.0001),('d',.01,0)]:
        save(f'contraction-wi{wi:g}-{kind}',dict(problem='contraction',mesh_file='meshes/contraction-paper.msh',dt=dt,end_time=10,rho=lam/18 if lam else .05,nu=lam/9 if lam else 1,mu=1,**{'lambda':lam,'diffusion':diff,'newton_tolerance':1e-12,'max_newton':15}),dict(suite='contraction',wi=wi,case=kind,dt=dt,diffusion=diff,steps=round(10/dt),note='Explicit paper coefficients imply alpha=0.9, despite text stating 8/9.'))
(repo/'configs/paper/manifest.json').write_text(json.dumps(cases,indent=2)+'\n')
print(f'Prepared {len(cases)} cases: 42 manufactured + 24 contraction. No simulations launched.')
