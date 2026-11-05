"""Reproduce Stage 4B.3A measurements without changing historical audit data."""
import argparse
import contextlib
import hashlib
import io
import json
import platform
import re
import subprocess
import time
import unittest
from pathlib import Path
import numpy as np
import scipy
from test_corrected_heating import CorrectedHeatingTests,HISTORICAL,BENCHMARKS,ROOT
from heating_energy_reference import mp


def clean(v):
    if isinstance(v,np.ndarray):return clean(v.tolist())
    if isinstance(v,(complex,np.complexfloating)):return [float(v.real),float(v.imag)]
    if isinstance(v,(float,np.floating)):return float(v) if np.isfinite(v) else None
    if isinstance(v,(bool,np.bool_)):return bool(v)
    if isinstance(v,(int,np.integer)):return int(v)
    if isinstance(v,dict):return {str(k):clean(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [clean(x) for x in v]
    return v


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--build-comparison',type=Path)
    parser.add_argument('--full-log',type=Path);parser.add_argument('--require-high-precision',action='store_true');args=parser.parse_args()
    if args.require_high_precision and mp is None:parser.error('mpmath must be importable for the independent high-precision reference')
    log=io.StringIO();started=time.perf_counter()
    with contextlib.redirect_stdout(log):
        result=unittest.TextTestRunner(stream=log,verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(CorrectedHeatingTests))
    if not result.wasSuccessful():raise RuntimeError(log.getvalue())
    cases={}
    for key,p in CorrectedHeatingTests.cases.items():
        name=key.split('/')[0];sp,kw=BENCHMARKS[name]
        c=HISTORICAL['cases'][f'{name}/False']
        cases[key]=dict(species=sp,kpar=kw.get('kpar',.06),kperp=kw.get('kperp',.05),a=kw.get('a',.01),mode=p,
            independent_H=c['energy_reference'],independent_H_species=c['species_energy'],
            independent_Q=c['reference_q'],audit_historical_H=c['energy_historical'])
    full=None
    if args.full_log:
        text=args.full_log.read_text();m=re.search(r'Ran (\d+) tests in ([\d.]+)s',text)
        full=dict(tests=int(m[1]),seconds=float(m[2]),passed=bool(re.search(r'^OK(?:\s|$)',text,re.M))) if m else dict(passed=False)
    protected=['src/kinetic_species.inc','src/kinetic_collisional_species.inc','src/kinetic_core.f90',
        'src/cold_susceptibility.f90','src/species_moments.f90','src/electromagnetic_fields.f90','src/adaptive_kinetic.f90',
        'src/bessel_cache.f90','src/harmonic_selection.f90','src/plasma_z_backend.f90','src/fpc.f90']
    output=dict(stage='4B.3A corrected collisionless heating and complete wave energy',
        toolchain=dict(gfortran=subprocess.check_output(['gfortran','--version'],text=True).splitlines()[0],platform=platform.platform(),
            numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__ if mp else None),
        tests=result.testsRun,skipped=len(result.skipped),validation_seconds=time.perf_counter()-started,full_suite=full,
        conventions=dict(energy='16*pi*W/abs(deltaEx)^2',Ps='Im(E† chi_s(wr) E)/H',
            complex_encoding='[real,imaginary]; null=unavailable/nonfinite',
            derivative='centered/Richardson; matrix and E contraction; two adjacent differences <=1e-6; best stable window',
            applicability='bit1:abs(gamma/wr)>0.05; bit2:independent weak-balance warning; formal H retained; Ps unavailable',
            checked_flags='-O0 -g -fdefault-real-8 -fcheck=all -ffpe-trap=invalid,zero,overflow',
            new_module_flags='production -O3 ... -fno-fast-math'),
        clean_build=json.loads(args.build_comparison.read_text()) if args.build_comparison else None,
        measurements=CorrectedHeatingTests.measurements,cases=cases,
        preserved_source_hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in protected},
        historical_evidence='MIXED_HEATING_ENERGY_AUDIT_RESULTS.json (unchanged)')
    args.output.write_text(json.dumps(clean(output),indent=2,allow_nan=False)+'\n')
    print(f'Wrote {args.output}; {result.testsRun} tests passed')

if __name__=='__main__':main()
