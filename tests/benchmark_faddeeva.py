"""Native-code standalone checkpoint timings; no Python callback per Z call."""
import argparse,json,statistics,subprocess
from pathlib import Path
from faddeeva_harness import FaddeevaHarness,ROOT
from plasma_dispersion_audit import json_clean


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    h=FaddeevaHarness()
    report={'compiler':subprocess.check_output(['gfortran','--version'],text=True).splitlines()[0],
            'c_compiler':subprocess.check_output(['cc','--version'],text=True).splitlines()[0],
            'samples':3,'z':[],'susceptibility':[]}
    distributions={'small':[.1+.01j,1+.1j,2-.05j],
                   'transition':[3.01+1e-12j,3.01-1e-12j,-2.947038584816192+.5930573003500077j],
                   'large':[10+.1j,1e3+.01j,1e6-.001j]}
    for label,points in distributions.items():
        old=[];new=[]
        for _ in range(3):
            old.append(sum(h.z(z,repeats=200000)['seconds'] for z in points))
            new.append(sum(h.accurate_z(z,repeats=200000)['seconds'] for z in points))
        a=statistics.median(old);b=statistics.median(new)
        report['z'].append({'name':label,'arguments':points,'calls':600000,'legacy_seconds':a,
                            'faddeeva_seconds':b,'ratio':b/a})
    cases=[('small',(.4,1.,1.7,1.,1.,0.),.07,8),
           ('stress',(2.,.25,1.3,.5,.05,.4),3.,128),
           ('large',(1.,1.,1.,1.,1.,0.),50.,512)]
    for label,s,kp,n in cases:
        old=[];new=[]
        for _ in range(3):
            old.append(h.tensor(s,.73+.13j,kpar=-.31,kperp=kp,n=n,repeats=200)['seconds'])
            new.append(h.accurate_tensor(s,.73+.13j,kpar=-.31,kperp=kp,n=n,repeats=200)['seconds'])
        a=statistics.median(old);b=statistics.median(new)
        report['susceptibility'].append({'name':label,'species':s,'kperp':kp,'n':n,'calls':200,
                                         'legacy_seconds':a,'faddeeva_seconds':b,'ratio':b/a})
    args.output.write_text(json.dumps(json_clean(report),indent=2)+'\n')
    h.close();print('Wrote',args.output)

if __name__=='__main__':main()
