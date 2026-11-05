from faddeeva_build import c_objects
"""Reproducible J.1 root-search benchmark; timings are not pass/fail assertions.

Run: MPLCONFIGDIR=/tmp/plume-mpl python3 tests/benchmark_bessel.py
Writes measurements to stdout as JSON. Compile the same production probe used
by accuracy tests. Include process wall time and solver CPU time separately.
"""
import json
from pathlib import Path
import platform
import statistics
import subprocess
import tempfile
import time
import numpy as np
from scipy.optimize import root
from test_adaptive_dispersion import AdaptiveDispersionTests,ROOT,MODULES
from mixed_reference import wave_matrix


def main():
    subprocess.run(['make','-j2'],cwd=ROOT,check=True,capture_output=True)
    p=(1.,1.,1.,1.,1.,0.); e=(.2,1836.,1.,-1.,1.,0.)
    cases=[('small',[p,e],.06,.05,.055-.001j),
           ('moderate',[p,e],.01,3.,1.17-.0001j),
           ('large',[p,e],.01,50.,1.01-.005j),
           ('three_species',[p,(.2,1836.,.8,-1.,1.1,.02727272727),(.03,.25,1.2,.5,.05,.3)],
            .06,.05,.055-.001j)]
    report={'compiler':subprocess.check_output(['gfortran','--version'],text=True).splitlines()[0],
            'platform':platform.platform(),'repeated_roots':100,'median_runs':3,'cases':[]}
    with tempfile.TemporaryDirectory() as directory:
        exe=Path(directory)/'probe'
        subprocess.run(['gfortran','-O3','-fdefault-real-8','-I',str(ROOT/'include'),
                        str(ROOT/'tests/harmonic_probe.f90'),*[str(ROOT/f'{m}.o') for m in MODULES],*c_objects(ROOT),
                        '-o',str(exe)],check=True,capture_output=True)
        runner=AdaptiveDispersionTests();runner.probes={False:exe}
        for name,species,kpar,kp,guess in cases:
            args=dict(kpar=kpar,kperp=kp,va=.02)
            def residual(z):
                m=wave_matrix(species,complex(*z),**args,harmonics=1024)
                d=np.linalg.det(m/np.max(abs(m)));return [d.real,d.imag]
            r=root(residual,[guess.real,guess.imag],tol=1e-10)
            if not r.success: raise RuntimeError((name,r.message))
            expected=complex(*r.x)
            item={'name':name,'kpar':kpar,'kperp':kp,'independent_root':[expected.real,expected.imag],
                  'measurements':{}}
            for mode,count,repeats,label in [('fixed',128,100,'fixed128'),
                                             ('fixed',1024,100,'fixed1024'),
                                             ('adaptive',128,100,'adaptive'),
                                             ('adaptive',128,-100,'adaptive_uncached')]:
                wall=[];cpu=[]
                for _ in range(3):
                    start=time.perf_counter()
                    blocks,roots,stdout,files=runner.evaluate(species,[expected*.999+1e-6j],mode=mode,
                           fixed=count,repeats=repeats,roots=True,cap=512,**args)
                    wall.append(time.perf_counter()-start)
                    cpu.append(float(next(l.split()[1] for l in stdout.splitlines() if l.startswith('CPU_SECONDS'))))
                flag,w=roots[0]
                if not 0 <= flag < 75: raise RuntimeError((name,label,'root did not converge',roots))
                matrix=blocks[0]['lam']
                data={'wall_seconds':statistics.median(wall),'solver_cpu_seconds':statistics.median(cpu),
                      'root_iterations':flag,'root':[w.real,w.imag],
                      'root_error':[abs(w.real-expected.real),abs(w.imag-expected.imag)],
                      'normalized_residual':float(abs(blocks[0]['det'])/np.max(abs(matrix))**3),
                      'metadata':files['test.bessel'].splitlines()}
                item['measurements'][label]=data
            item['adaptive_speedup_vs_fixed128']=item['measurements']['fixed128']['solver_cpu_seconds']/item['measurements']['adaptive']['solver_cpu_seconds']
            item['cache_speedup']=item['measurements']['adaptive_uncached']['solver_cpu_seconds']/item['measurements']['adaptive']['solver_cpu_seconds']
            report['cases'].append(item)
    print(json.dumps(report,indent=2))

if __name__=='__main__': main()
