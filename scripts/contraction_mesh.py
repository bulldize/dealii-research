#!/usr/bin/env python3
"""Gmsh base mesh followed by five conforming local edge-bisection passes.

The original FEniCS mesh is unavailable: this is a documented mesh reconstruction,
not a claim that the exact triangulation matches the authors' mesh.
"""
import argparse, hashlib, json
from pathlib import Path
import gmsh
import meshio
import numpy as np

parser=argparse.ArgumentParser()
parser.add_argument('--h',type=float,default=.05)
parser.add_argument('--levels',type=int,default=5)
parser.add_argument('--output',type=Path,default=Path('meshes/contraction-paper.msh'))
a=parser.parse_args()
if a.h<=0 or a.levels<0: raise ValueError('invalid size/refinement count')
a.output.parent.mkdir(parents=True,exist_ok=True)
gmsh.initialize()
gmsh.option.setNumber('General.Terminal',0)
gmsh.model.add('contraction')
coords=[(-10,-2),(0,-2),(0,-.5),(20,-.5),(20,.5),(0,.5),(0,2),(-10,2)]
points=[gmsh.model.geo.addPoint(x,y,0,a.h) for x,y in coords]
lines=[gmsh.model.geo.addLine(points[i],points[(i+1)%8]) for i in range(8)]
loop=gmsh.model.geo.addCurveLoop(lines);gmsh.model.geo.addPlaneSurface([loop]);gmsh.model.geo.synchronize()
gmsh.option.setNumber('Mesh.MeshSizeMin',a.h);gmsh.option.setNumber('Mesh.MeshSizeMax',a.h)
gmsh.option.setNumber('Mesh.Algorithm',6);gmsh.option.setNumber('General.NumThreads',1)
gmsh.model.mesh.generate(2)
tags,xyz,_=gmsh.model.mesh.getNodes();xyz=np.array(xyz).reshape(-1,3)
lookup={int(t):i for i,t in enumerate(tags)}
_,nodes=gmsh.model.mesh.getElementsByType(2)
cells=np.array([lookup[int(x)] for x in nodes]).reshape(-1,3)
gmsh.finalize()
vertices=xyz[:,:2].tolist();history=[len(cells)]
for level in range(a.levels):
    xy=np.array(vertices);centers=xy[cells].mean(axis=1)
    distance=np.minimum(np.linalg.norm(centers-[0,.5],axis=1),np.linalg.norm(centers-[0,-.5],axis=1))
    edges=[tuple(sorted((int(c[i]),int(c[(i+1)%3])))) for c in cells for i in range(3)]
    tri_edges=[edges[i:i+3] for i in range(0,len(edges),3)]
    marked=set()
    for idx in np.flatnonzero(distance<.6*2**(-level)):
        es=tri_edges[idx]
        marked.add(max(es,key=lambda e:np.linalg.norm(xy[e[0]]-xy[e[1]])))
    # Two split edges are closed to red refinement; all neighboring shared edges agree.
    while True:
        additions=set()
        for es in tri_edges:
            if sum(e in marked for e in es)==2: additions.update(es)
        additions-=marked
        if not additions: break
        marked.update(additions)
    mids={}
    for edge in sorted(marked):
        mids[edge]=len(vertices);vertices.append(((xy[edge[0]]+xy[edge[1]])/2).tolist())
    new=[]
    for c,es in zip(cells,tri_edges):
        selected=[i for i,e in enumerate(es) if e in marked]
        if not selected:new.append(c.tolist())
        elif len(selected)==1:
            k=selected[0];u,v,w=map(int,[c[k],c[(k+1)%3],c[(k+2)%3]]);m=mids[es[k]]
            new.extend([[u,m,w],[m,v,w]])
        else:
            u,v,w=map(int,c);uv,vw,wu=[mids[e] for e in es]
            new.extend([[u,uv,wu],[uv,v,vw],[wu,vw,w],[uv,vw,wu]])
    cells=np.array(new,dtype=int);history.append(len(cells))
xy=np.array(vertices);p=xy[cells];u=p[:,1]-p[:,0];v=p[:,2]-p[:,0];signed=u[:,0]*v[:,1]-u[:,1]*v[:,0]
assert np.min(np.abs(signed))>1e-14
reverse=signed<0;cells[reverse]=cells[reverse][:,[0,2,1]]
area=np.sum(np.abs(signed))/2
assert abs(area-60)<1e-8,area
# Every edge must be either shared by two triangles or on the physical boundary.
counts={}
for c in cells:
    for i in range(3):
        edge=tuple(sorted((int(c[i]),int(c[(i+1)%3]))));counts[edge]=counts.get(edge,0)+1
assert all(n in (1,2) for n in counts.values())
for edge,n in counts.items():
    if n!=1:continue
    x,y=xy[list(edge)].mean(axis=0)
    on_boundary=abs(x+10)<1e-9 or abs(x-20)<1e-9 or (x<0 and abs(abs(y)-2)<1e-9) or (x>0 and abs(abs(y)-.5)<1e-9) or (abs(x)<1e-9 and .5-1e-9<=abs(y)<=2+1e-9)
    assert on_boundary,(x,y)
meshio.write(a.output,meshio.Mesh(np.column_stack([xy,np.zeros(len(xy))]),[('triangle',cells)]),file_format='gmsh22',binary=False)
meta={'base_h':a.h,'levels':a.levels,'vertices':len(xy),'cells':len(cells),'cells_after_each_pass':history,'area':float(area),'gmsh_algorithm':6,'refinement':'longest edge bisection plus conforming red closure; barycenter distance marking','original_fenics_mesh_available':False,'mesh_sha256':hashlib.sha256(a.output.read_bytes()).hexdigest()}
a.output.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n')
print(json.dumps(meta,indent=2))
