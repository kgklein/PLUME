"""Reproduce native validation/production timing tables without editing old audits.

MPLCONFIGDIR=/tmp/plume-mpl python3 tests/measure_faddeeva_backend.py --output PATH
SciPy is a test reference only. Every timed loop runs in compiled Fortran/C.
"""
import argparse
import contextlib
import io
import json
import platform
import statistics
import subprocess
import unittest
from pathlib import Path
import numpy as np
from faddeeva_harness import FaddeevaHarness,audit_arguments
from plasma_dispersion_audit import ROOT,scipy_z,json_clean
from test_faddeeva_backend import FaddeevaBackendTests
from plasma_dispersion_audit import MODULES
from faddeeva_build import c_objects


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    report={'compiler':subprocess.check_output(['gfortran','--version'],text=True).splitlines()[0],
            'c_compiler':subprocess.check_output(['cc','--version'],text=True).splitlines()[0],
            'platform':platform.platform(),'regimes':{},'production_timing':[]}
    h=FaddeevaHarness()
    for regime,z in audit_arguments(h):
        row=report['regimes'].setdefault(regime,dict(count=0,range_failures=0,max_absolute=0.,max_relative=0.,max_scaled=0.))
        row['count']+=1
        out=h.accurate_z(z)
        with np.errstate(all='ignore'):expected=scipy_z(1.,z)
        if not np.isfinite(expected):
            assert out['status']!=0;row['range_failures']+=1;continue
        assert out['status']==0
        absolute=abs(out['value']-expected);scaled=absolute/max(1.,abs(expected))
        row['max_absolute']=max(row['max_absolute'],absolute)
        row['max_scaled']=max(row['max_scaled'],scaled)
        if abs(expected)>1e-280:row['max_relative']=max(row['max_relative'],absolute/abs(expected))
        assert scaled<1e-13
    h.close()
    # Collect deterministic production accuracy/root measurements from assertions.
    stream=io.StringIO()
    with contextlib.redirect_stdout(stream):
        suite=unittest.defaultTestLoader.loadTestsFromTestCase(FaddeevaBackendTests)
        result=unittest.TextTestRunner(stream=stream,verbosity=2).run(suite)
    if not result.wasSuccessful():raise RuntimeError(stream.getvalue())
    report['production_validation']=FaddeevaBackendTests.measurements
    # Separate setup keeps all benchmark objects alive after unittest cleanup.
    FaddeevaBackendTests.setUpClass()
    runner=FaddeevaBackendTests()
    executable=FaddeevaBackendTests.h.directory/'benchmark_uninstrumented'
    subprocess.run(['gfortran','-O3','-funroll-loops','-ffast-math','-fdefault-real-8','-I',str(ROOT/'include'),
                    str(ROOT/'tests/harmonic_probe.f90'),*[str(ROOT/f'{m}.o') for m in MODULES],
                    *c_objects(ROOT),'-lm','-o',str(executable)],check=True,capture_output=True)
    runner.probes={False:executable}
    report['timing_flags']='-O3 -funroll-loops -ffast-math -fdefault-real-8; strict C/wrapper objects'
    report['timing_instrumentation']='CPU_TIME native loops; no Z counter or Python callback'

    p=(1.,1.,1.,1.,1.,0.);e=(.2,1836.,1.,-1.,1.,0.)
    try:
        for name,species,kp,n in [('small',[p,e],.05,128),('moderate',[p,e],3.,128),
                                 ('large',[p,e],50.,512),
                                 ('cold_reference',[(0.,1.,0.,1.,1.,0.),e,(.03,.25,1.2,.5,.05,.3)],3.,128)]:
            item={'case':name,'mode':'fixed','harmonics':n,'repeats':250,'median_samples':3}
            for backend in ['legacy','faddeeva']:
                times=[]
                for _ in range(3):
                    _,_,stdout,_=runner.evaluate(species,[.73+.13j],backend=backend,mode='fixed',fixed=n,
                                                kpar=.31,kperp=kp,va=.02,repeats=250)
                    times.append(float(next(line.split()[1] for line in stdout.splitlines() if line.startswith('CPU_SECONDS'))))
                item[backend+'_seconds']=statistics.median(times)
            item['ratio']=item['faddeeva_seconds']/item['legacy_seconds']
            report['production_timing'].append(item)
    finally:FaddeevaBackendTests.doClassCleanups()
    args.output.write_text(json.dumps(json_clean(report),indent=2)+'\n')
    print('Wrote',args.output)


if __name__=='__main__':main()
