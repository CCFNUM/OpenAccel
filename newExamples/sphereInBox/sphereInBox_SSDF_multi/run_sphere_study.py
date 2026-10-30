#!/usr/bin/env python3
"""Run the ten separate-mesh SSDF cases and build sphere_refinement_study.xlsx.
Run with: /usr/bin/python3 run_sphere_study.py
Requires the OpenAccel environment (oa), numpy, scipy, h5py, PyYAML and VTK.
--analyze-only reuses completed runs. --threads defaults to 4; --repeats to 1.
Runtime is whole OpenAccel process wall time, including reading, masking and I/O.
Nearest triangle connectivity is independently recovered with VTK, not reported
by SSDF itself; equally near triangles may have different valid IDs.
"""
import argparse,csv,hashlib,json,math,os,shutil,subprocess,tempfile,time
from pathlib import Path
import numpy as np
import yaml,h5py,vtk
from scipy.io import netcdf_file
from scipy.spatial import cKDTree
from vtk.util.numpy_support import numpy_to_vtk,numpy_to_vtkIdTypeArray

ROOT=Path(__file__).resolve().parent
SIZES=[.020,.018,.016,.014,.012,.010,.009,.008,.006,.004]
HEADERS=['query_node_index','query_node_id','x_m','y_m','z_m',
         'triangle_node_id_1','triangle_node_id_2','triangle_node_id_3',
         'analytical_signed_distance_m','ssdf_signed_distance_m',
         'error_m','absolute_error_m']

class Exodus:
    def __init__(self,path):
        self.f=h5py.File(path,'r') if h5py.is_hdf5(path) else netcdf_file(path,'r',mmap=False)
        self.v=self.f if isinstance(self.f,h5py.File) else self.f.variables
    def get(self,k,s=slice(None)): return np.array(self.v[k][s],copy=True)
    def names(self,k): return [r.tobytes().split(b'\0')[0].decode().strip() for r in self.get(k)]
    def coords(self): return self.get('coord').T if 'coord' in self.v else np.column_stack([self.get('coord'+a) for a in 'xyz'])
    def ids(self,n): return self.get('node_num_map').astype(np.int64) if 'node_num_map' in self.v else np.arange(1,n+1)
    def close(self): self.f.close()

def digest(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''): h.update(b)
    return h.hexdigest()

def write_rows(path,header,rows):
    tmp=path.with_suffix(path.suffix+'.partial')
    with tmp.open('w',newline='') as f:
        w=csv.writer(f); w.writerow(header); w.writerows(rows)
    os.replace(tmp,path)

def analyze(folder,background,radius):
    sphere=next(folder.glob('sphere_*.exo'))
    m=Exodus(sphere); sxyz=m.coords(); sids=m.ids(len(sxyz))
    ss=m.names('ss_names').index('sphere_wall')+1
    e=m.get(f'elem_ss{ss}').astype(int)-1; sides=m.get(f'side_ss{ss}').astype(int)-1
    conn=m.get('connect1').astype(int)-1
    face=np.array([[0,1,3],[1,2,3],[0,3,2],[0,2,1]])
    if conn.shape[1]!=4 or np.any(sides<0) or np.any(sides>3): raise ValueError('Expected TET4 sphere.')
    tri=np.take_along_axis(conn[e],face[sides],axis=1)
    ntri=len(tri); ntet=len(conn); m.close()
    verts,inv=np.unique(tri,return_inverse=True); local=inv.reshape(-1,3)
    if np.max(np.abs(np.linalg.norm(sxyz[verts],axis=1)-radius))>1e-8: raise ValueError('Sphere radius/centre mismatch.')
    edges=np.sort(np.concatenate([local[:,[0,1]],local[:,[1,2]],local[:,[2,0]]]),axis=1)
    if not np.all(np.unique(edges,axis=0,return_counts=True)[1]==2): raise ValueError('Sphere surface is not closed.')
    poly=vtk.vtkPolyData(); points=vtk.vtkPoints()
    points.SetData(numpy_to_vtk(np.ascontiguousarray(sxyz[verts]),deep=True));poly.SetPoints(points)
    cells=vtk.vtkCellArray()
    cells.SetCells(ntri,numpy_to_vtkIdTypeArray(np.column_stack([np.full(ntri,3),local]).astype(np.int64).ravel(),deep=True))
    poly.SetPolys(cells);locator=vtk.vtkStaticCellLocator();locator.SetDataSet(poly);locator.BuildLocator()
    b=Exodus(background); xyz=b.coords(); ids=b.ids(len(xyz)); b.close()
    result=Exodus(folder/'results.e'); rxyz=result.coords()
    fields=result.names('name_nod_var')
    if 'signed_distance' not in fields: raise ValueError('Missing signed_distance; rebuild the masking-based SSDF executable.')
    values=result.get(f"vals_nod_var{fields.index('signed_distance')+1}",-1)
    d,ix=cKDTree(rxyz).query(xyz)
    if d.max()>1e-10 or len(np.unique(ix))!=len(xyz) or len(rxyz)!=len(xyz): raise ValueError('Background query points differ.')
    measured=values[ix]; result.close()
    if not np.isfinite(measured).all(): raise ValueError('Nonfinite signed distances.')
    analytical=np.linalg.norm(xyz,axis=1)-radius
    err=measured-analytical; ae=np.abs(err)
    rows=[]; triangle_error=[]
    closest=[0.,0.,0.];cid=vtk.reference(0);sub=vtk.reference(0);d2=vtk.reference(0.);cell=vtk.vtkGenericCell()
    for i,q in enumerate(xyz):
        locator.FindClosestPoint(q,closest,cell,cid,sub,d2)
        t=int(cid)
        if t<0 or not math.isfinite(float(d2)): raise ValueError('Nearest triangle search failed.')
        triangle_error.append(abs(abs(measured[i])-math.sqrt(float(d2))))
        rows.append([i+1,int(ids[i]),*map(float,q),*map(int,sids[tri[t]]),
                     float(analytical[i]),float(measured[i]),float(err[i]),float(ae[i])])
    write_rows(folder/'comparison.csv',HEADERS,rows)
    stats={'nodes':len(xyz),'triangles':ntri,'tetrahedra':ntet,
           'rmse':float(np.sqrt(np.mean(err**2))),'mae':float(ae.mean()),
           'bias':float(err.mean()),'stddev':float(err.std(ddof=0)),
           'max_abs':float(ae.max()),'median_abs':float(np.median(ae)),
           'p95_abs':float(np.percentile(ae,95)),
           'sign_mismatches':int(np.sum((analytical<0)!=(measured<0))),
           'triangle_max_error':float(max(triangle_error))}
    return rows,stats

def workbook(payload):
    # Use the desktop's bundled spreadsheet authoring runtime from WSL.
    node=Path('/mnt/c/Users/user/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node.exe')
    builder=ROOT/'build_study_workbook.mjs'
    with tempfile.TemporaryDirectory(prefix='sphere-study-report-') as temp:
        data=Path(temp)/'study.json';data.write_text(json.dumps(payload,allow_nan=False))
        win=lambda p:subprocess.check_output(['wslpath','-w',str(p)],text=True).strip()
        subprocess.run([str(node),'--max-old-space-size=8192',win(builder),win(data),win(ROOT/'sphere_refinement_study.xlsx')],check=True,cwd=ROOT)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--threads',type=int,default=4);ap.add_argument('--repeats',type=int,default=1)
    ap.add_argument('--analyze-only',action='store_true');ap.add_argument('--skip-workbook',action='store_true')
    ap.add_argument('--exe',type=Path,default=Path('/home/alimans/HSLU/OpenAccel/build-ssdf/openaccel-3D.exe'))
    a=ap.parse_args()
    if a.threads<1 or a.repeats<1: ap.error('Threads/repeats must be positive.')
    background=ROOT/'box.exo';base=yaml.safe_load((ROOT/'input.i').read_text());radius=.1675
    env=dict(os.environ,OMP_NUM_THREADS=str(a.threads),OMP_DYNAMIC='FALSE',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    cases=[]
    executable_hash=digest(a.exe)
    background_hash=digest(background)
    for i,h in enumerate(SIZES,1):
        label=f'sphere_{i:02d}';folder=ROOT/f'case_{i:02d}';folder.mkdir(exist_ok=True)
        sphere=folder/(label+'.exo')
        if not sphere.exists(): shutil.move(str(ROOT/(label+'.exo')),str(sphere))
        config=json.loads(json.dumps(base));config.pop('immersed_distance',None)
        config['mesh']['file_path']=str(background)
        pa=config['simulation']['physical_analysis']
        domain=next(d for d in pa['domains'] if d['name']=='fluid')
        domain['location']=['background']
        domain['boundaries']=[b for b in domain['boundaries'] if set(b.get('location',[]))<=set(['inlet','outlet','symmetry'])]
        pa.setdefault('overrides',{})['masking']=[{'name':'sphere','domains':['fluid'],'shape':{'option':'mesh','path':sphere.name}}]
        solver=config['simulation']['solver']
        solver['solver_control']['expert_parameters'].update(disable_physics=True,force_wall_distance_calculation=True,wall_distance_method='signed_distance_function')
        solver['output_control'].update(file_path='results.e',output_frequency=1,output_fields=['signed_distance'])
        content=yaml.safe_dump(config,sort_keys=False)
        inputfile=folder/'input.i'
        if not a.analyze_only: inputfile.write_text(content)
        elif inputfile.read_text()!=content: raise ValueError(f'{label}: input settings changed; rerun instead.')
        timing=folder/'runtime.csv'
        provenance=[digest(inputfile),background_hash,digest(sphere),executable_hash]
        times=[]
        if not a.analyze_only:
            if timing.exists(): raise FileExistsError(f'{timing} exists; use --analyze-only to reuse results.')
            for rep in range(1,a.repeats+1):
                print(f'Running {label}, repetition {rep}/{a.repeats}, OMP={a.threads}',flush=True)
                with (folder/f'run_{rep:02d}.log').open('w') as out:
                    start=time.perf_counter()
                    proc=subprocess.run([str(a.exe),'-i','input.i'],cwd=folder,env=env,stdout=out,stderr=subprocess.STDOUT)
                    elapsed=time.perf_counter()-start
                if proc.returncode:
                    raise RuntimeError(f'{label} exited {proc.returncode}: '+(folder/f'run_{rep:02d}.log').read_text()[-4000:])
                if not (folder/'results.e').exists(): raise RuntimeError(f'{label}: missing results.e')
                times.append(elapsed)
            write_rows(timing,['repetition','elapsed_seconds','omp_threads','mpi_ranks','input_sha256','background_sha256','sphere_sha256','executable_sha256'],
                       [[j+1,t,a.threads,1,*provenance] for j,t in enumerate(times)])
        else:
            with timing.open() as f: records=list(csv.DictReader(f))
            for r in records:
                if [r[k] for k in ['input_sha256','background_sha256','sphere_sha256','executable_sha256']]!=provenance: raise ValueError('Source files changed since run.')
                if int(r['omp_threads'])!=a.threads: raise ValueError('Thread count differs from recorded run.')
            times=[float(r['elapsed_seconds']) for r in records]
        rows,stats=analyze(folder,background,radius)
        stats.update(case=label,h=h,runtime_seconds=float(np.median(times)),runtime_min=min(times),runtime_max=max(times),repeats=len(times),threads=a.threads)
        cases.append({'name':label,'rows':rows,'stats':stats})
        print(f"{label}: RMSE={stats['rmse']:.9g} m; elapsed={stats['runtime_seconds']:.3f} s; triangles={stats['triangles']}",flush=True)
    payload={'radius':radius,'headers':HEADERS,'cases':cases}
    if not a.skip_workbook: workbook(payload)
    print('All ten cases complete.',flush=True)
if __name__=='__main__': main()
