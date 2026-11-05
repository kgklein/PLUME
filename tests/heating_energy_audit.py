"""Isolated audit harness, measurements and source snapshots; production is read-only."""
import json
import hashlib
import subprocess
from pathlib import Path
import numpy as np
from scipy.optimize import root
from plasma_dispersion_audit import ROOT,MODULES
from faddeeva_build import c_objects
from heating_energy_reference import hermitian,contraction,absorption,channels,fd_derivative,energy_reference
from mixed_reference import wave_matrix
import test_mixed_scans as scans
import test_mixed_dispersion as mixed

P,E,CP,CE,A,CA=scans.P,scans.E,scans.CP,scans.CE,scans.A,scans.CA


def complex_values(f):
    a=np.array(list(map(float,f)));return a[::2]+1j*a[1::2]


def parse(text,nspec):
    points=[]
    for line in text.splitlines():
        f=line.split()
        if not f:continue
        if f[0]=='POINT':
            points.append(dict(status=int(f[1]),w=complex(float(f[2]),float(f[3])),e=np.empty(3,complex),b=np.empty(3,complex),
              chi=np.empty((nspec,3,3),complex),after=np.empty((nspec,3,3),complex),power=np.empty(nspec),split=np.empty((nspec,7)),
              low=np.empty((nspec,3,3,3),complex),after_low=np.empty((nspec,3,3,3),complex),recover=np.full((nspec,3,3,3),np.nan+1j*np.nan)))
        elif f[0]=='FIELD':points[-1]['e'][int(f[1])-1],points[-1]['b'][int(f[1])-1]=complex_values(f[2:])
        elif f[0]=='HEAT':points[-1]['power'][int(f[1])-1]=float(f[2]);points[-1]['split'][int(f[1])-1]=list(map(float,f[3:]))
        elif f[0] in ['CHI','AFTER']:points[-1][f[0].lower()][int(f[1])-1,int(f[2])-1,int(f[3])-1]=complex_values(f[4:])[0]
        elif f[0] in ['LOW','RECOVER']:
            index=(int(f[1])-1,int(f[2])+1,int(f[3])-1,int(f[4])-1)
            if f[0]=='LOW':points[-1]['low'][index],points[-1]['after_low'][index]=complex_values(f[5:])
            else:points[-1]['recover'][index]=complex_values(f[5:])[0]
    return points


class HeatingHarness:
    def __init__(self):
        self.hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'src').glob('*') if p.is_file()}
        scans.MixedScanTests.setUpClass()
        self.tmp=scans.MixedScanTests.tmp;self.directory=Path(self.tmp.name)
        self.checked=scans.MixedScanTests.checked;self.executables=scans.MixedScanTests.executables.copy();self.probes={}
        for checked in [False,True]:
            directory=self.checked if checked else ROOT
            flags=['-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow'] if checked else ['-O3']
            result=subprocess.run(['gfortran','-fdefault-real-8',*flags,'-fno-fast-math','-I',str(directory if checked else ROOT/'include'),
                 str(ROOT/'tests/heating_energy_probe.f90'),*[str(directory/f'{m}.o') for m in MODULES],*c_objects(directory),
                 '-o',str(self.directory/f'heating{checked}')],capture_output=True,text=True)
            if result.returncode:raise RuntimeError(result.stderr)
            self.probes[checked]=self.directory/f'heating{checked}'

    def close(self):
        for p,digest in self.hashes.items():
            assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==digest,f'Production source changed: {p}'
        self.tmp.cleanup()

    def evaluate(self,species,frequencies,modes=4,checked=False,**kwargs):
        if np.isscalar(modes):modes=[modes]*len(frequencies)
        stdin=''.join(f'{w.real:.17g} {w.imag:.17g} {m}\n' for w,m in zip(frequencies,modes))
        result,files=mixed.run(self.probes[checked],scans.source(species,[],**kwargs),stdin)
        if result.returncode:raise RuntimeError(result.stdout+result.stderr)
        result=parse(result.stdout,len(species))
        if len(result)!=len(frequencies):raise AssertionError('Missing probe output')
        return result

    def canary_executable(self):
        """Audit build only: deterministic NaN initialization exposes unread Ew."""
        flags=['-O0','-g','-fdefault-real-8','-fcheck=all','-finit-real=nan','-I',str(self.checked),'-J',str(self.directory)]
        obj=self.directory/'nan_disprels.o'
        result=subprocess.run(['gfortran',*flags,'-c',str(ROOT/'src/disprels.f90'),'-o',str(obj)],capture_output=True,text=True)
        if result.returncode:raise RuntimeError(result.stderr)
        objects=[str(obj) if m=='disprels' else str(self.checked/f'{m}.o') for m in MODULES]
        exe=self.directory/'nan_plume'
        result=subprocess.run(['gfortran',*flags,str(ROOT/'src/plume.f90'),*objects,str(self.checked/'fpc.o'),*c_objects(self.checked),'-o',str(exe)],capture_output=True,text=True)
        if result.returncode:raise RuntimeError(result.stderr)
        return exe


def audit_case(harness,species,guess=.06-.0001j,kpar=.06,kperp=.05,a=.01,checked=False,mode='adaptive',normalization='d_ref'):
    point=harness.evaluate(species,[guess],modes=1,checked=checked,kpar=kpar,kperp=kperp,va=a,mode=mode,normalization=normalization)[0]
    if point['status']!=0:raise AssertionError('Unavailable electromagnetic gauge')
    w=point['w'];wr=w.real;e=point['e'];b=point['b']
    reference=energy_reference(species,wr,kpar,kperp,a,e)
    steps=[1e-3,1e-4,1e-5,1e-6,1e-7]
    frequencies=[complex(wr)]+[complex(wr+sign*wr*h) for h in steps for sign in [-1,1]]+[complex(wr*1.000001)]
    snapshots=harness.evaluate(species,frequencies,checked=checked,kpar=kpar,kperp=kperp,va=a,mode=mode,normalization=normalization)
    lookup={p['w'].real:p['chi'] for p in snapshots}
    evaluate=lambda x:lookup[x]
    target=reference['derivative']
    center=[];forward=[];errors=[]
    norm=lambda d:float(np.linalg.norm(d-target)/max(1.,np.linalg.norm(target)))
    for h in steps:
        central=fd_derivative(evaluate,wr,h,True);one=fd_derivative(evaluate,wr,h,False)
        center.append(float(sum(contraction(e,central))));forward.append(float(sum(contraction(e,one))))
        errors.append(dict(step=h,centered_matrix=norm(central),forward_matrix=norm(one),
                           centered_energy=abs(center[-1]-sum(reference['contributions']))/max(1.,abs(sum(reference['contributions']))),
                           forward_energy=abs(forward[-1]-sum(reference['contributions']))/max(1.,abs(sum(reference['contributions'])))))
    raw=absorption(e,snapshots[0]['chi']);H=reference['energy']
    # Reproduce the actual legacy frequency product and sum-before-difference,
    # not just an algebraically equivalent w+w*h/per-species finite difference.
    legacy_derivative=(1.000001*hermitian(np.sum(snapshots[-1]['chi'],axis=0))-
                       hermitian(np.sum(snapshots[0]['chi'],axis=0)))/1e-6
    legacy_particle=float(contraction(e,legacy_derivative))
    reconstructed=legacy_particle+float(np.vdot(b,b).real)
    observed=None
    if all(s[0]>0 for s in species):
        nonzero=point['power']!=0.
        if np.any(nonzero):observed=float(np.median(raw[nonzero]/point['power'][nonzero]))
    historical=observed if observed is not None else reconstructed
    corrected=legacy_particle+reference['vacuum']+reference['magnetic']
    gamma=-wr*sum(raw)/H
    historical_gamma=-wr*sum(raw)/historical
    root_consistency=(abs(gamma-w.imag)/abs(w.imag) if abs(w.imag)>1e-12 else None)
    def function(x):
        matrix=wave_matrix(species,complex(*x),kpar,kperp,a,256);value=np.linalg.det(matrix/np.max(abs(matrix)))
        return [value.real,value.imag]
    independent=root(function,[w.real,w.imag],tol=1e-11)
    if not independent.success and np.linalg.norm(function(independent.x))>1e-12:raise AssertionError(independent.message)
    eigen_residual=np.linalg.norm(wave_matrix(species,w,kpar,kperp,a,256)@e)/(np.linalg.norm(wave_matrix(species,w,kpar,kperp,a,256))*np.linalg.norm(e))
    return dict(w=w,e=e,b=b,raw=raw,energy_reference=H,energy_historical=historical,energy_corrected_forward=corrected,
                vacuum=reference['vacuum'],magnetic=reference['magnetic'],species_energy=reference['contributions'],
                species_model=['kinetic' if s[0]>0 else 'cold-fluid (energy only)' for s in species],
                vacuum_fraction=reference['vacuum']/H,derivative_errors=errors,
                historical_normalizer_observed=observed,historical_normalizer_reconstructed=reconstructed,
                legacy_derivative_matrix_error=float(np.linalg.norm(legacy_derivative-np.sum(target,axis=0))/max(1.,np.linalg.norm(np.sum(target,axis=0)))),
                legacy_derivative_energy_error=abs(legacy_particle-sum(reference['contributions']))/max(1.,abs(sum(reference['contributions']))),
                predicted_gamma=gamma,historical_predicted_gamma=historical_gamma,
                closure_relative=root_consistency,weak_ratio=abs(w.imag/wr),
                independent_root=complex(*independent.x),independent_root_error=abs(w-complex(*independent.x)),null_residual=eigen_residual,
                reference_q=reference['q'],point=point,reference=reference)
