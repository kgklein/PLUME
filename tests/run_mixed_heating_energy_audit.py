"""Run read-only numerical audit and write observations, never golden outputs.

Optional high precision: install mpmath in any importable test environment;
--require-high-precision refuses to substitute the independent finite difference.
"""
import argparse
import contextlib
import io
import json
import platform
import subprocess
import unittest
from pathlib import Path
import numpy as np
import scipy
from heating_energy_reference import mp
from test_heating_energy_audit import HeatingEnergyAuditTests,BENCHMARKS
from plasma_dispersion_audit import ROOT


def clean(value):
    if isinstance(value,np.ndarray):return clean(value.tolist())
    if isinstance(value,(complex,np.complexfloating)):return [clean(float(value.real)),clean(float(value.imag))]
    if isinstance(value,(float,np.floating)):return float(value) if np.isfinite(value) else None
    if isinstance(value,(bool,np.bool_)):return bool(value)
    if isinstance(value,(int,np.integer)):return int(value)
    if isinstance(value,dict):return {str(k):clean(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [clean(v) for v in value]
    return value


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--require-high-precision',action='store_true');args=p.parse_args()
    if args.require_high_precision and mp is None:p.error('mpmath is an optional audit dependency; make it importable first')
    log=io.StringIO()
    with contextlib.redirect_stdout(log):
        suite=unittest.defaultTestLoader.loadTestsFromTestCase(HeatingEnergyAuditTests)
        result=unittest.TextTestRunner(stream=log,verbosity=2).run(suite)
    if not result.wasSuccessful():raise RuntimeError(log.getvalue())
    cases={}
    for key,d in HeatingEnergyAuditTests.cases.items():
        species,parameters=BENCHMARKS[key.split('/')[0]]
        cases[key]={k:v for k,v in d.items() if k not in ['point','reference']}
        cases[key]['species']=species;cases[key]['inputs']=dict(kpar=.06,kperp=.05,a=.01,**{k:v for k,v in parameters.items() if k not in ['kpar','kperp','a']})
        cases[key]['inputs'].update({k:v for k,v in parameters.items() if k in ['kpar','kperp','a']})
        cases[key]['historical_output']=dict(power=d['point']['power'],channels=d['point']['split'])
    report=dict(stage='4B.3A audit only; production unchanged',toolchain=dict(
        compiler=subprocess.check_output(['gfortran','--version'],text=True).splitlines()[0],platform=platform.platform(),
        numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__ if mp is not None else None),
        tests=result.testsRun,skipped=len(result.skipped),reference_derivative=('mpmath 65-digit analytic' if mp is not None else 'independent SciPy five-point real difference'),
        conventions=dict(complex_encoding='[real, imaginary]; null means unavailable/nonfinite',
                         energy='W*16*pi/abs(deltaEx)**2',raw_q='E^dagger chi_A(omega_r) E',
                         dimensional_power='Omega_ref*omega_r*abs(deltaEx)**2/(8*pi)*raw_q',
                         legacy_P='raw_q / (magnetic norm at complex mode + forward derivative of omega*sum(chi_H))'),
        measurements=HeatingEnergyAuditTests.measurements,cases=cases,
        exploratory_failure=dict(species=[(.01,1.,1.,1.,1.,0.),(.01,4.,1.,-1.,1.,0.)],kpar=.2,kperp=.1,a=.3,
            guess=8.+.01j,status='production secant did not refine in additional high-frequency trial',
            independent_real_frequency=7.461253,normalized_determinant_scale_approximately=1e-11,
            handling='Not used in weak-mode acceptance; no root or susceptibility changes made'),
        production_hashes=HeatingEnergyAuditTests.h.hashes)
    args.output.write_text(json.dumps(clean(report),indent=2,allow_nan=False)+'\n')
    print(f'Wrote {args.output}; tests={result.testsRun}, skipped={len(result.skipped)}')


if __name__=='__main__':main()
