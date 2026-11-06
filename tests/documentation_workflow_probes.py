"""Read-only Stage 5A examples and downstream reader characterization.

Run against an isolated checkout/executable; original inputs/readers are never
edited. Marker and real-output evidence is written under --artifacts only.
"""
import argparse
import contextlib
from concurrent.futures import ThreadPoolExecutor
import io
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import time
import numpy as np
from documentation_output_audit import parse_input


def example(checkout,path,artifacts,timeout):
    text=path.read_text();groups=parse_input(text);params=next(g['values'] for g in groups if g['group']=='params')
    started=time.perf_counter();name=path.stem;location=artifacts/'examples'/name;location.mkdir(parents=True,exist_ok=True)
    work=location/'run';work.mkdir(exist_ok=True);(work/'data'/str(params.get('dataname','example'))).mkdir(parents=True,exist_ok=True)
    inputpath=work/'case.in';inputpath.write_text(text)
    try:
        out=subprocess.run([str(checkout/'plume.e'),'case.in'],cwd=work,capture_output=True,text=True,timeout=timeout)
        code=out.returncode;stdout=out.stdout;stderr=out.stderr;timed_out=False
    except subprocess.TimeoutExpired as e:
        code=None;stdout=e.stdout or '';stderr=e.stderr or '';timed_out=True
        if isinstance(stdout,bytes):stdout=stdout.decode(errors='replace')
        if isinstance(stderr,bytes):stderr=stderr.decode(errors='replace')
    (location/'stdout.txt').write_text(stdout);(location/'stderr.txt').write_text(stderr)
    files=[]
    for p in sorted((work/'data').rglob('*')):
        if not p.is_file():continue
        f=dict(name=str(p.relative_to(work)),bytes=p.stat().st_size)
        if re.search(r'\.(?:mode\d+|roots|map)$',p.name):
            lines=[l.split() for l in p.read_text().splitlines() if l.strip() and not l.startswith('#')]
            lengths=sorted(set(map(len,lines)));f.update(rows=len(lines),column_counts=lengths)
            if lines and len(lengths)==1:
                try:
                    a=np.array(lines,float);f['nonfinite_cells']=int(np.count_nonzero(~np.isfinite(a)))
                    if p.suffix!='.map' and a.shape[1]>=6:f['first_root']=[float(x) if np.isfinite(x) else None for x in a[0,4:6]]
                except ValueError:f['parse_error']='nonnumeric row'
        if p.suffix=='.energy':
            from linfpclib.linfpc import load_plume_energy
            try:
                a=load_plume_energy(p);f.update(records=len(a['record']),energy_status_counts={str(x):int(sum(a['energy_status']==x)) for x in set(a['energy_status'])},
                    applicability_counts={str(x):int(sum(a['applicability']==x)) for x in set(a['applicability'])},
                    f_E_range=[float(np.nanmin(a['f_E'])),float(np.nanmax(a['f_E']))] if np.isfinite(a['f_E']).any() else None)
            except Exception as e:f['parse_error']=str(e)
        files.append(f)
    return dict(path=str(path.relative_to(checkout)),returncode=code,timeout=timed_out,seconds=time.perf_counter()-started,
        root_failures=stdout.count('scan root failed'),unavailable=stdout.lower().count('unavailable'),
        energy_warnings=stdout.count('WARNING: corrected heating/energy'),stderr_tail=stderr[-800:],stdout_tail=stdout[-800:],files=files,
        interpretation='execution/finite-row check only; branch identity is not inferred from exit code')


def markers(checkout,artifacts):
    sys.path.insert(0,str(checkout))
    from linfpclib.linfpc import load_plume_sweep,load_plume_sweep_debug,load_plume_sweep_orig,load_plume_energy,plume_input,loadeigen
    result=[];directory=artifacts/'readers';directory.mkdir(exist_ok=True)
    keys=['ttd_yy','ttd_yz','ld_zy','ld_zz','n_eq_0','cd_n_p','cd_n_m']
    for n in [1,2,3,4]:
        for eigen in [False,True]:
            for heating in [False,True]:
                front=18+8*n if eigen else 6;length=front+(8*n if heating else 0)+6*n
                row=np.arange(length,dtype=float)+1000.;path=directory/f'n{n}-e{eigen}-h{heating}.mode1';np.savetxt(path,[row])
                for label,reader in [('current',load_plume_sweep),('debug',load_plume_sweep_debug)]:
                    log=io.StringIO();errors=[]
                    try:
                        with contextlib.redirect_stdout(log):a=reader(path,nspec=n,eigen=eigen,heating=heating)
                        expected=row if label=='current' else np.arange(length,dtype=float)
                        for index,key in enumerate(['kperp','kpar','betap','vtp','w','g']):
                            if a[key][0]!=expected[index]:errors.append('prefix '+key)
                        if eigen:
                            for index,key in enumerate(['bxr','bxi','byr','byi','bzr','bzi','exr','exi','eyr','eyi','ezr','ezi'],6):
                                if a[key][0]!=expected[index]:errors.append('field '+key)
                        for s in range(n):
                            if heating:
                                if a[f'p{s+1}'][0]!=expected[front+s]:errors.append('species total '+str(s+1))
                                for c,k in enumerate(keys):
                                    if a[f'p{s+1}{k}'][0]!=expected[front+n+7*s+c]:errors.append(k+str(s+1))
                            params=front+(8*n if heating else 0)
                            for j,k in enumerate(['tau','mu','alph','q','D','vv']):
                                if a[f'p{s+1}{k}'][0]!=expected[params+6*s+j]:errors.append('parameter '+k+str(s+1))
                            if eigen:
                                for c,key in enumerate(['ux','uy','uz']):
                                    for part,offset in [('r',0),('i',1)]:
                                        if a[f'{key}{s+1}{part}'][0]!=expected[18+6*s+2*c+offset]:errors.append('moment '+key+part+str(s+1))
                                for part,offset in [('r',0),('i',1)]:
                                    if a[f'n{s+1}{part}'][0]!=expected[18+6*n+2*s+offset]:errors.append('density '+part+str(s+1))
                        result.append(dict(reader=label,nspec=n,eigen=eigen,heating=heating,columns=length,passed=not errors,errors=errors))
                    except Exception as e:result.append(dict(reader=label,nspec=n,eigen=eigen,heating=heating,passed=False,error=repr(e)))
    # The documented current loader assumes low_n=true: characterize instead
    # of asserting it should accept every ambiguous historical schema.
    for n in [2,3]:
        for eigen in [False,True]:
            length=(18+8*n if eigen else 6)+n+6*n+int(eigen)
            path=directory/f'low-false-n{n}-e{eigen}.mode1';np.savetxt(path,[np.arange(length)+1000.])
            log=io.StringIO()
            try:
                with contextlib.redirect_stdout(log):a=load_plume_sweep(path,nspec=n,heating=True,eigen=eigen)
                outcome='returned (incorrect schema is not explicitly rejected)'
            except Exception as e:outcome=repr(e)
            result.append(dict(reader='current',case='low_n=false',nspec=n,eigen=eigen,outcome=outcome))
    # Preserves unavailable cold heat / failed root NaNs.
    path=directory/'nan.mode1';row=np.arange(62,dtype=float);row[34]=np.nan;row[36:43]=np.nan
    np.savetxt(path,[row]);a=load_plume_sweep(path,nspec=2,eigen=True,heating=True)
    result.append(dict(reader='current',case='cold heat NaN',passed=bool(np.isnan(a['p1'][0]) and np.isnan(a['p1ttd_yy'][0]))))
    # Repeated 2D blocks have blank separators; the generic reader flattens them.
    path=directory/'grid.mode1';path.write_text(' '.join(map(str,row))+'\n\n'+' '.join(map(str,row))+'\n')
    a=load_plume_sweep(path,nspec=2,eigen=True,heating=True)
    result.append(dict(reader='current',case='2D blank separators',rows=len(a['w']),grid_indices_retained=False))
    inputs=[]
    for use_map in [False,True]:
        for writeOut in [False,True]:
            obj=plume_input('audit');obj.set_params(.4,.05,.06,.01,2,0,1,1,use_map,writeOut,normalization='d_ref')
            inputs.append(dict(case='set_params writeOut',use_map=use_map,requested=writeOut,actual=obj.params['writeOut'],correct=obj.params['writeOut']==('.true.' if writeOut else '.false.')))
    # Reader has a hard-coded two-species debug schema and fabricated Pe=0.
    path=directory/'fpc-debug.eigen';np.savetxt(path,[np.arange(62)+1000.]);a=loadeigen(path)
    result.append(dict(reader='loadeigen',case='unique marker',p_e=[a['Pe'].real,a['Pe'].imag],expected_markers=[1036,1037],hardcoded_two_species=True))
    return dict(markers=result,wrapper_parameters=inputs)


def additional_readers(checkout,artifacts,examples):
    from contours.plot_plume_resonances import read_plume_input,read_dispersion_output
    from linfpclib.linfpc import load_plume_energy
    extras=[]
    for name in ['mixed_fields_warm_reference','mixed_fields_cold_reference']:
        inp=read_plume_input(checkout/'inputs/example'/f'{name}.in')
        file=next((artifacts/'examples'/name/'run/data/example').glob('*.mode1'))
        disp=read_dispersion_output(file,inp,1)
        extras.append(dict(case='contours eigen-only schema',input=name,input_heating=inp.heating,inferred_heating=disp.heating_available,
            heat_start=disp.heat_start,actual_columns=disp.data.shape[1],false_heating_values=disp.data[0,disp.heat_start:disp.heat_start+len(inp.species)].tolist()))
    association=[];actual_readers=[]
    for e in examples:
        if e['timeout']:continue
        root=artifacts/'examples'/Path(e['path']).stem/'run'
        groups=parse_input((root/'case.in').read_text());params=next(g['values'] for g in groups if g['group']=='params')
        axes=[g['values'] for g in groups if g['group'].startswith('scan_input_')]
        eigen=axes[-1].get('eigen',False) if axes else False;heating=axes[-1].get('heating',False) if axes else False;n=params['nspec']
        from linfpclib.linfpc import load_plume_sweep
        for modefile in root.glob('data/*/*.mode*'):
            if '_tensor_' in modefile.name:continue
            raw=np.loadtxt(modefile,ndmin=2);loaded=load_plume_sweep(modefile,nspec=n,eigen=eigen,heating=heating)
            errors=[]
            for key,col in [('kperp',0),('kpar',1),('w',4),('g',5)]:
                if not np.allclose(loaded[key],raw[:,col],rtol=0,atol=0,equal_nan=True):errors.append(key)
            front=18+8*n if eigen else 6
            if heating:
                for species in range(n):
                    if not np.allclose(loaded[f'p{species+1}'],raw[:,front+species],rtol=0,atol=0,equal_nan=True):errors.append('Ps'+str(species+1))
                    for c,key in enumerate(['ttd_yy','ttd_yz','ld_zy','ld_zz','n_eq_0','cd_n_p','cd_n_m']):
                        if not np.allclose(loaded[f'p{species+1}{key}'],raw[:,front+n+7*species+c],rtol=0,atol=0,equal_nan=True):errors.append(key+str(species+1))
            actual_readers.append(dict(example=e['path'],file=modefile.name,rows=len(raw),eigen=eigen,heating=heating,errors=errors))
        for file in root.glob('data/*/*.energy'):
            energy=load_plume_energy(file);counters={};errors=[];failed=0;cache={}
            for i,source in enumerate(energy['source_file']):
                if source=='stdout':continue
                if source not in cache:cache[source]=np.loadtxt(root/source,ndmin=2)
                rows=cache[source];index=counters.get(source,0);counters[source]=index+1
                if index>=len(rows):errors.append('missing mode row');continue
                row=rows[index]
                if energy['root_valid'][i]==0:
                    failed+=1
                    if not(np.isnan(row[4:6]).all() and energy['energy_status'][i]==8):errors.append('invalid-root mismatch')
                elif not(np.allclose(row[:2],[energy['kperp'][i],energy['kpar'][i]],rtol=1e-6,atol=1e-13) and np.allclose(row[4:6],[energy['omega_r'][i],energy['gamma'][i]],rtol=1e-6,atol=1e-12)):
                    errors.append('frequency/coordinate mismatch')
            association.append(dict(example=e['path'],records=len(energy['record']),failed=failed,association_errors=len(errors),errors=errors[:3]))
    return dict(plotter_characterization=extras,energy_association=association,current_actual_output=actual_readers)


def main():
    p=argparse.ArgumentParser();p.add_argument('--checkout',type=Path,required=True);p.add_argument('--artifacts',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--timeout',type=float,default=20);a=p.parse_args()
    sys.path.insert(0,str(a.checkout));a.artifacts.mkdir(parents=True,exist_ok=True)
    probes=markers(a.checkout,a.artifacts)
    files=sorted((a.checkout/'inputs/example').glob('*.in'))
    with ThreadPoolExecutor(max_workers=3) as pool:r=list(pool.map(lambda path:example(a.checkout,path,a.artifacts,a.timeout),files))
    probes.update(additional_readers(a.checkout,a.artifacts,r))
    a.output.write_text(json.dumps(dict(readers=probes,examples=r,timeout_seconds=a.timeout),indent=2)+'\n')
    print('Example outcomes:',[(q['path'],q['returncode'],q['timeout'],q['root_failures']) for q in r])
    print('Marker failures:',[q for q in probes['markers'] if q.get('passed') is False])

if __name__=='__main__':main()
