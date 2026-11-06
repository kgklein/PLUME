"""Stage 5A characterization: accepted inputs, map grids, FPC readers, flags.

No production changes. All fixtures and actual outputs live under --artifacts.
"""
import argparse
import contextlib
import io
import json
from pathlib import Path
import subprocess
import sys
import numpy as np


def main():
    p=argparse.ArgumentParser();p.add_argument('--checkout',type=Path,required=True);p.add_argument('--artifacts',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    sys.path[:0]=[str(a.checkout),str(a.checkout/'tests')]
    from linfpclib.linfpc import loadlinfpccepar,load_plume_energy,load_plume_sweep
    import test_mixed_scans as scans
    import test_mixed_dispersion as mixed
    from test_normalization import input_text
    a.artifacts.mkdir(parents=True,exist_ok=True);records=[]
    def run(name,text):
        folder=a.artifacts/name;folder.mkdir(exist_ok=True)
        (folder/'case.in').write_text(text)
        for data in ['mixed','normalization']:(folder/'data'/data).mkdir(parents=True,exist_ok=True)
        r=subprocess.run([str(a.checkout/'plume.e'),'case.in'],cwd=folder,capture_output=True,text=True,timeout=25)
        (folder/'stdout.txt').write_text(r.stdout);(folder/'stderr.txt').write_text(r.stderr)
        return folder,r
    text=scans.source([scans.P,scans.E],[],option=0)+'''&maps
loggridw=.true.
loggridg=.true.
omi=0.1
omf=1.0
gami=-0.0001
gamf=-0.01
nr=2
ni=2
positive_roots=.true.
/
'''
    folder,r=run('logarithmic-map',text);file=next((folder/'data/mixed').glob('*.map'));grid=np.loadtxt(file)
    records.append(dict(case='logarithmic map endpoint semantics',returncode=r.returncode,requested_real=[.1,1.],actual_real=[float(grid[:,0].min()),float(grid[:,0].max())],
        requested_imag=[-.0001,-.01],actual_imag=[float(grid[:,1].max()),float(grid[:,1].min())],shape=list(grid.shape),
        source='src/disprels.f90::map_search; dr/di use log10 bounds but origins use omi/abs(gami)'))
    folder=a.artifacts/'missing-output-directory';folder.mkdir(exist_ok=True);(folder/'case.in').write_text(scans.source([scans.P,scans.E],[],option=1))
    missing=subprocess.run([str(a.checkout/'plume.e'),'case.in'],cwd=folder,capture_output=True,text=True,timeout=10)
    records.append(dict(case='fresh run without data subdirectory',returncode=missing.returncode,stderr_tail=missing.stderr[:450],expected_workflow='mkdir -p data/<dataName> before running; writeOut does not suppress required metadata files'))
    for first,second in [(False,True),(True,False),(True,True)]:
        axes=[scans.scan(0,.05,.051,steps=1,tensor=first),scans.scan(1,.06,.061,steps=1,tensor=second)]
        folder,r=run('grid-tensor-'+str(first)+'-'+str(second),scans.source([scans.P,scans.CE],axes,option=2))
        records.append(dict(case='2D tensor flags',first=first,second=second,returncode=r.returncode,tensor_files=[q.name for q in (folder/'data/mixed').glob('*tensor*')]))
    for name,header in [('uppercase','&SPECIES_1'),('comment','&species_1 ! protons')]:
        text=scans.source([scans.P,scans.E],[],option=1).replace('&species_1',header)
        folder,r=run('indexed-header-'+name,text)
        records.append(dict(case='indexed header '+name,returncode=r.returncode,stdout_tail=r.stdout[-350:],stderr_tail=r.stderr[:250],source='input_unit_exist compares the full stripped line to exact lowercase &species_N'))
    # Zero-rate forbidden geometry is not fixed by an arbitrary small kpar.
    for name,sp,mode,kwargs in [('warm-exact-parallel',[scans.P,scans.E],-1,dict(kperp=0)),('warm-exact-perpendicular',[scans.P,scans.E],-1,dict(kpar=0)),('cold-exact-perpendicular',[scans.CP,scans.CE],-1,dict(kpar=0)),('negative-density-warm',[scans.P,mixed.replace_species(scans.E,density=-1)],1,{})]:
        folder,r=run(name,scans.source(sp,[],option=mode,**kwargs))
        records.append(dict(case=name,returncode=r.returncode,stdout_tail=r.stdout[-900:],stderr_tail=r.stderr[-350:],note='accepted process exit is not proof of valid tensor/eigenfunction'))
    for flag in ['eigen','tensor']:
        axis=scans.scan(1,.06,.061,steps=1,**{flag:True})
        text=scans.source([scans.P,scans.E],[axis],option=4)+'''&maps
loggridw=.false.
loggridg=.false.
omi=.04
omf=.08
gami=-.001
gamf=0
nr=4
ni=4
positive_roots=.true.
/
'''
        folder,r=run('maps-'+flag,text);files=sorted(q.name for q in (folder/'data/mixed').glob('*'))
        records.append(dict(case='map requested '+flag,returncode=r.returncode,files=files,omitted_without_error=r.returncode==0 and not any('.mode' in f for f in files),note='map_search writes determinant/root maps, never supplementary eigen/tensor files'))
    # Actual gyro FPC output and intentional NaN marker in an isolated copy.
    sp=[(1,1,1.3,1,1,0),(.5,1836,1,-1,1,.02)]
    text=input_text('rho_ref',.4,1.3,.05,.06,sp).replace('option=-1','option=6')+'''&fpc
vperpmin=0
vperpmax=0.5
vparmin=-0.5
vparmax=0.5
delv=0.5
/
&guess_1
g_om=.08
g_gam=-.001
/
'''
    folder,r=run('fpc-gyro',text);file=next((folder/'data/normalization').glob('*.cpar.specie01.mode01'))
    actual=loadlinfpccepar(str(file));lines=file.read_text().splitlines();marker=folder/'nan.cpar.specie01.mode01'
    lines[5]='NaN '+' '.join(lines[5].split()[1:]);marker.write_text('\n'.join(lines)+'\n');loaded=loadlinfpccepar(str(marker))
    records.append(dict(case='gyro FPC reader missing values',returncode=r.returncode,actual_shape=list(actual['CEpar'].shape),injected='NaN',reader_value=float(loaded['CEpar'][0,0]),nan_preserved=bool(np.isnan(loaded['CEpar'][0,0])),source='linfpclib/linfpc.py::loadlinfpccepar'))
    from linfpclib.linfpc import plume_input
    fixture=a.artifacts/'roundtrip.in'
    fixture.write_text(scans.source([scans.P,scans.E],[scans.scan(1,.06,.061,tensor=True)],extra='low_n=.false.')+"&maps\nloggridw=.false.\nloggridg=.false.\nomi=.04\nomf=.08\ngami=-.001\ngamf=0\nnr=4\nni=5\npositive_roots=.true.\n/\n&fpc\nvxmax=1\nvxmin=-1\nvymax=1\nvymin=-1\nvzmax=1\nvzmin=-1\ndelv=.5\n/\n")
    obj=plume_input('audit')
    with contextlib.redirect_stdout(io.StringIO()):obj.load_from_file(str(fixture))
    records.append(dict(case='wrapper roundtrip coverage',missing_params=['low_n'] if 'low_n' not in obj.params else [],
        missing_maps=[k for k in ['loggridg','nr','ni'] if k not in obj.maps],
        missing_scan=[k for k in ['tensor'] if k not in obj.scan_inputs[0]],
        missing_cartesian_bounds=[k for k in ['vxmin','vxmax','vymin','vymax','vzmin','vzmax'] if k not in obj.fpc]))
    # Exact corrected-mode file layout and public loader association are retained
    # from actual example output, including interleaved/failed rows.
    a.output.write_text(json.dumps(records,indent=2)+'\n');print(json.dumps(records,indent=2))

if __name__=='__main__':main()
