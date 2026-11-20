"""Reproduce J1.B.1 measurements. No production files are written.

PYTHONPATH=/tmp/plume-z-python MPLCONFIGDIR=/tmp/plume-mpl \
 python3 tests/run_plasma_dispersion_audit.py --output /tmp/plume-z-results.json --high-precision

Install mpmath in any importable environment for optional --high-precision.
"""
import argparse
import csv
import json
import platform
import statistics
import sys
import tempfile
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import root
from scipy.special import ive
from plasma_dispersion_audit import *
from mixed_reference import species_tensor,wave_matrix
import test_mixed_dispersion as mixed


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--high-precision',action='store_true')
    parser.add_argument('--arguments-csv',type=Path)
    args=parser.parse_args()
    h=AuditHarness()
    report={'toolchain':{'compiler':subprocess.check_output(['gfortran','--version'],text=True).splitlines()[0],
                        'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,
                        'platform':platform.platform(),'flags':h.flags},
            'source_hashes':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in (ROOT/'src').glob('*') if p.is_file()},
            'z_summary':[],'stress':[],'roots':[]}
    grid=z_grid(); sampled=[];worst=None;worst_data=None
    for label,system in stress_systems():
        beta,alpha=normalization_state(system)
        for sign in [-1,1]:
            for kp in [.07,3.,30.]:
                for w in [.73+.13j,.73-.08j]:
                    for index,s in enumerate(system):
                        if s[0]==0:continue
                        kw=dict(kpar=sign*.31,kperp=kp,va=.02,beta=beta,alpha_scale=alpha,n=1024)
                        historical=h.tensor(s,w,**kw)
                        high=h.tensor(s,w,method=2,**kw)
                        ref=species_tensor(s,w,sign*.31,kp,.02,1024)
                        point={'system':label,'species_index':index+1,'species':s,'omega':w,'kpar':sign*.31,
                               'kperp':kp,'historical_error':metrics(historical['chi'],ref),
                               'wofz_substituted_error':metrics(high['chi'],ref),
                               'z_error':metrics(historical['z'],scipy_z(sign,historical['arguments'])),
                               'flags':historical['flags']}
                        report['stress'].append(point)
                        sampled.extend((sign*.31,z) for z in historical['arguments'])
                        if worst is None or point['historical_error']['scaled']>worst['historical_error']['scaled']:
                            worst=point;worst_data=(historical,high)
    report['worst_stress']=worst
    s=tuple(worst['species']);w=worst['omega'];kp=worst['kperp'];kpar=worst['kpar']
    beta,alpha=normalization_state(dict(stress_systems())[worst['system']])
    kw=dict(kpar=kpar,kperp=kp,va=.02,beta=beta,alpha_scale=alpha,n=1024)
    historical,high=worst_data
    c,b,extra=stix_terms(s,w,kpar,kp,.02,1024)
    delta_z=historical['z']-high['z'];predicted=np.sum(b*delta_z,axis=-1)
    real_delta=historical['chi']-high['chi']
    component=np.unravel_index(np.argmax(abs(real_delta)/np.maximum(1,abs(high['chi']))),real_delta.shape)
    row,column=component
    contributions=c+b*high['z']
    denom=abs(high['chi'][component])
    report['error_budget']={'component':[row+1,column+1],
         'historical_chi':historical['chi'],'wofz_chi':high['chi'],
         'delta_chi':real_delta,'predicted_delta_from_z':predicted,
         'linear_budget_closure':metrics(predicted,real_delta),
         'harmonic_sum_condition':(np.sum(abs(contributions[row,column]))+abs(extra[row,column]))/denom,
         'raw_constant_plus_z_condition':(np.sum(abs(c[row,column])+abs(b[row,column]*high['z']))+abs(extra[row,column]))/denom,
         'relative_z_condition':np.sum(abs(b[row,column]*high['z']))/denom,
         'dominant_harmonics':[]}
    for i in np.argsort(abs(b[row,column]*delta_z))[-8:][::-1]:
        z=historical['arguments'][i]
        entry={'n':int(i-1024),'zeta':z,'effective_zeta':np.sign(kpar)*z,
               'delta_z':delta_z[i],'delta_chi_component':b[row,column,i]*delta_z[i],
               'historical':h.z(z,k=kpar),'reference':scipy_z(kpar,z)}
        report['error_budget']['dominant_harmonics'].append(entry)
    report['tolerance_stress']=[];report['normalization_stress']=[];report['collision_stress']=[]
    ref=species_tensor(s,w,kpar,kp,.02,1024)
    for tol in TOLERANCES:
        d=h.tensor(s,w,method=1,tolerance=tol,**kw)
        report['tolerance_stress'].append({'tolerance':tol,'error':metrics(d['chi'],ref)})
    for bstar in [.01,.1,1.,10.]:
        kwargs={**kw,'beta':bstar,'alpha_scale':1.}
        old=h.tensor(s,w,**kwargs);accurate=h.tensor(s,w,method=2,**kwargs)
        report['normalization_stress'].append({'bstar':bstar,'historical_error':metrics(old['chi'],ref),
              'wofz_error':metrics(accurate['chi'],ref)})
    for nu in [0.,1e-4,.01,.1]:
        old=h.tensor(s,w,collision=nu,**kw);accurate=h.tensor(s,w,method=2,collision=nu,**kw)
        report['collision_stress'].append({'nu':nu,'historical_vs_wofz':metrics(old['chi'],accurate['chi'])})
    report['zero_collision_limit']=[]
    for species in [s,(s[0],s[1],1.,s[3],s[4],s[5])]:
        collisionless=h.tensor(species,w,method=2,**kw)
        col=h.tensor(species,w,method=2,collisional=True,**kw)
        tau,mu,al,q,density,u=species;g=mu/q;v=np.sqrt(mu*tau)
        orders=np.arange(1,1025)
        gn=ive(orders,kp*kp*q*q*al*tau/(2*mu))
        norm=2*density/(.02**2*w*w*q*q*kpar*kpar*al*tau)
        predicted=norm/(kpar*v)*np.sum(gn*orders**2*g*g*(1-al)*
             (1-orders*g*(collisionless['z'][1024+orders]-collisionless['z'][1024-orders])))
        report['zero_collision_limit'].append({'species':species,
             'collisional_minus_collisionless':col['chi']-collisionless['chi'],
             'scaled_difference':metrics(col['chi'],collisionless['chi']),
             # Historical counterfactual only: the corrected source should no
             # longer exhibit this old misgrouping. Keep the original JSON
             # snapshot immutable when comparing stages.
             'historical_parenthesis_only_difference':predicted,
             'historical_prediction_residual_on_current_source':
                 abs(col['chi'][2,2]-collisionless['chi'][2,2]-predicted)})
    report['near_real_continuity']=[]
    for x in [3.000000001,3.01,3.1,4.]:
        for y in [-1e-12,0.,1e-12]:
            z=complex(x,y);data=h.z(z)
            report['near_real_continuity'].append({'z':z,'historical':data,'reference':scipy_z(1.,z),
                  'error':metrics(data['value'],scipy_z(1.,z))})
    for k,z in dict.fromkeys(sampled):grid.append(('sampled_stress',np.sign(k)*z))
    if args.arguments_csv:
        with args.arguments_csv.open('w',newline='') as f:
            writer=csv.writer(f);writer.writerow(['kpar','zeta_real','zeta_imag'])
            writer.writerows((k,z.real,z.imag) for k,z in dict.fromkeys(sampled))
    report['z_point_count']=len(grid);report['sampled_argument_count']=len(set(sampled))
    worst_z=[]
    for method,tol in [(0,1e-7),*[(1,t) for t in TOLERANCES]]:
        buckets={}
        for regime,z in grid:
            d=h.z(z,method=method,tolerance=tol)
            with np.errstate(all='ignore'):refz=scipy_z(1.,z)
            bucket=buckets.setdefault(regime,{'count':0,'valid_count':0,'undefined_or_nonfinite':0,'hard_cap':0,
               'turning_stop':0,'overflow':0,'underflow':0,'absolute':0.,'relative':0.,'scaled':0.,'worst_z':0j})
            bucket['count']+=1
            if d['stop']==3:bucket['hard_cap']+=1
            if d['stop']==2:bucket['turning_stop']+=1
            bucket['overflow']+=d['flags'][2];bucket['underflow']+=d['flags'][3]
            if d['stop'] in (1,4) or d['nonfinite'] or not np.isfinite(refz):
                bucket['undefined_or_nonfinite']+=1;continue
            bucket['valid_count']+=1
            err=metrics(d['value'],refz)
            if err['scaled']>bucket['scaled']:bucket['worst_z']=z
            for name in ['absolute','relative','scaled']:bucket[name]=max(bucket[name],err[name])
            if method==0:worst_z.append((err['scaled'],z,d['value'],refz))
        for regime,values in buckets.items():report['z_summary'].append({'method':method,'tolerance':tol,'regime':regime,**values})
    report['worst_z_values']=sorted(worst_z,key=lambda r:-r[0])[:12]
    # Tolerance timing excludes Python/C crossing and includes unchanged branch logic.
    timings=[]
    timing_points=[.2+.1j,2.5+.1j,1.+1.j,3.01+.1j,10.+.1j,complex(-historical['arguments'][1027])]
    for method,tol in [(0,1e-7),*[(1,t) for t in TOLERANCES]]:
        seconds=[]
        for _ in range(3):seconds.append(sum(h.z(z,method=method,tolerance=tol,repeats=100000)['seconds'] for z in timing_points))
        timings.append({'method':method,'tolerance':tol,'evaluations':len(timing_points)*100000,'median_seconds':statistics.median(seconds)})
    report['z_timing']=timings
    report['tensor_timing']=[]
    for method,tol in [(0,1e-7),*[(1,t) for t in TOLERANCES]]:
        seconds=[h.tensor(s,w,method=method,tolerance=tol,repeats=100,**kw)['seconds'] for _ in range(3)]
        report['tensor_timing'].append({'method':method,'tolerance':tol,'evaluations':100,'median_seconds':statistics.median(seconds)})
    if args.high_precision:
        import mpmath as mp
        report['toolchain']['mpmath']=mp.__version__
        mp_z=[]
        points=[0j,1+.1j,3.01+.1j,1+1j,3.1+0j,10+.1j,-14j,
                *[d['effective_zeta'] for d in report['error_budget']['dominant_harmonics'][:4]]]
        for z in points:
            with mp.workdps(90):
                zmp=mp.mpc(z);value=complex(mp.j*mp.sqrt(mp.pi)*mp.exp(-zmp*zmp)*mp.erfc(-mp.j*zmp))
            mp_z.append({'z':z,'scipy_error':metrics(scipy_z(1.,z),value),'value':value})
        report['mpmath_z']=mp_z
        mp_tensor=high_precision_stix(s,w,kpar,kp,.02,n=128)
        gamma=ive(np.arange(130),kp*kp*s[3]**2*s[2]*s[0]/(2*s[1]))
        held=h.tensor(s,w,method=2,**{**kw,'n':128})
        promoted=high_precision_stix(s,w,kpar,kp,.02,n=128,gamma=gamma,zvalues=held['z'])
        report['mpmath_tensor']={'n':128,'dps':80,'value':mp_tensor,
             'scipy_stix_error':metrics(species_tensor(s,w,kpar,kp,.02,128),mp_tensor),
             'fortran_wofz_error':metrics(held['chi'],mp_tensor),
             'promoted_gamma_z_arithmetic':metrics(held['chi'],promoted)}
        # Characterize the largest discrepancy remaining after exact Z substitution;
        # this small-lambda electron case tests reference/assembly cancellation.
        residual=max(report['stress'],key=lambda p:p['wofz_substituted_error']['scaled'])
        es=tuple(residual['species']);ew=residual['omega'];ek=residual['kpar'];et=residual['kperp']
        eb,ea=normalization_state(dict(stress_systems())[residual['system']])
        exact=high_precision_stix(es,ew,ek,et,.02,n=16)
        held=h.tensor(es,ew,kpar=ek,kperp=et,va=.02,beta=eb,alpha_scale=ea,n=16,method=2)
        egamma=ive(np.arange(18),et*et*es[3]**2*es[2]*es[0]/(2*es[1]))
        promoted=high_precision_stix(es,ew,ek,et,.02,n=16,gamma=egamma,zvalues=held['z'])
        report['mpmath_remaining_residual']={'case':residual,'n':16,'dps':80,
            'exact_tensor':exact,'scipy_stix_error':metrics(species_tensor(es,ew,ek,et,.02,16),exact),
            'fortran_wofz_error':metrics(held['chi'],exact),
            'fortran_vs_promoted_gamma_z_arithmetic':metrics(held['chi'],promoted),
            'promoted_gamma_z_vs_exact':metrics(promoted,exact)}
    # Production evaluators/roots: link original objects, never replace their Z.
    subprocess.run(['make','-j2'],cwd=ROOT,check=True,capture_output=True)
    exe=h.directory/'production_probe'
    subprocess.run(['gfortran','-O3','-fdefault-real-8','-I',str(ROOT/'include'),
         str(ROOT/'tests/harmonic_probe.f90'),*[str(ROOT/f'{m}.o') for m in MODULES],'-o',str(exe)],check=True,capture_output=True)
    system=dict(stress_systems())['warm_pe_alpha']
    report['production_stress']=[]
    for mode,n in [('fixed',128),('fixed',1024),('adaptive',128)]:
        text=mixed.input_text(system,kpar=kpar,kperp=kp,va=.02,extra=f"bessel_mode='{mode}'\nbessel_nfixed={n}")
        result,files=mixed.run(exe,text,f'{w.real} {w.imag} 0 1\n')
        if result.returncode:raise RuntimeError(result.stderr)
        blocks,_=mixed.parse(result.stdout,len(system))
        prod=blocks[0]['chi'][2]
        report['production_stress'].append({'mode':mode,'n':n,'reference_error':metrics(prod,ref),
                  'held_bessel_error':metrics(prod,historical['chi'])})
    text=mixed.input_text(system,kpar=kpar,kperp=kp,va=.02,
         extra="collision_type=1\nKn=1e30\nbessel_nfixed=1024")
    result,_=mixed.run(exe,text,f'{w.real} {w.imag} 0 1\n')
    blocks,_=mixed.parse(result.stdout,len(system))
    if result.returncode or not blocks:raise RuntimeError(result.stderr+result.stdout)
    held=h.tensor(s,w,collisional=True,**kw)['chi']
    report['production_zero_collision_confirmation']={
         'nu':0.,'held_vs_production':metrics(held,blocks[0]['chi'][2]),
         'production_collisional_chi':blocks[0]['chi'][2]}
    roots=[('stable_cold' ,[(0.,1.,0.,1.,1.,0.),(0.,1836.,0.,-1.,1.,0.)],.06,.05,.001,.06+1e-5j),
           ('damped',[(0.,1.,0.,1.,1.,0.),(.2,1836.,.8,-1.,1.,0.)],.06,.05,.001,.059-1e-5j),
           ('growing',[(0.,1.,0.,1.,1.,0.),(0.,1836.,0.,-1.,1.1,.2/1.1),(.03,.25,1.,.5,.05,2.)],.3,.1,.01,.3+.08j)]
    for label,sp,kp_r,kt_r,va,guess in roots:
        text=mixed.input_text(sp,kpar=kp_r,kperp=kt_r,va=va,extra="bessel_nfixed=1024")
        result,_=mixed.run(exe,text,f'{guess.real} {guess.imag} 1 1\n')
        blocks,rs=mixed.parse(result.stdout,len(sp))
        if not rs or not 0<=rs[0][0]<75:raise RuntimeError((label,result.stdout,result.stderr))
        w0=rs[0][1]
        report['roots'].append(matrix_diagnostic(h,label,sp,w0,kp_r,kt_r,va))
    report['roots'].append(matrix_diagnostic(h,'cancellation_stress',system,w,kpar,kp,.02))
    h.close()
    args.output.write_text(json.dumps(json_clean(report),indent=2,allow_nan=False)+'\n')
    print('Wrote',args.output,'worst stress',worst['historical_error']['scaled'])


def matrix_diagnostic(h,label,species,w,kpar,kperp,va):
    def matrix_at(w,method):
        responses=[]
        beta,alpha=normalization_state(species) if any(s[0]>0 for s in species) else (1.,1.)
        for s in species:
            if s[0]==0:responses.append(species_tensor(s,w,kpar,kperp,va,1024))
            else:responses.append(h.tensor(s,w,kpar=kpar,kperp=kperp,va=va,beta=beta,alpha_scale=alpha,n=1024,method=method)['chi'])
        n=np.array([kperp,0.,kpar])/(va*w)
        return np.eye(3)*(1-n@n)+np.outer(n,n)+sum(responses)
    old=matrix_at(w,0);new=matrix_at(w,2);scale=np.max(abs(new))
    d0=np.linalg.det(old/scale);d1=np.linalg.det(new/scale)
    shifts=[]
    for step in [1e-6,5e-7,2e-6]:
        dx=step*max(1.,abs(w))
        derivative=(np.linalg.det(matrix_at(w+dx,2)/scale)-np.linalg.det(matrix_at(w-dx,2)/scale))/(2*dx)
        shifts.append((d0-d1)/derivative if derivative else complex(np.nan,np.nan))
    return {'label':label,'species':species,'omega':w,'kpar':kpar,'kperp':kperp,'va':va,
            'matrix_error':metrics(old,new),'normalized_det_historical':d0,'normalized_det_wofz':d1,
            'cond2_historical':np.linalg.cond(old),'cond2_wofz':np.linalg.cond(new),
            'linear_root_shift_estimates':shifts if label!='cancellation_stress' else None,
            'nonroot_iso_determinant_shifts':shifts if label=='cancellation_stress' else None,
            'is_root_case':label!='cancellation_stress'}

if __name__=='__main__':main()
