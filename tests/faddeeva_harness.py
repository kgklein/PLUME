"""Standalone Faddeeva checkpoint and native-code performance adapter."""
import ctypes as ct
import json
from pathlib import Path
import subprocess
import numpy as np
from scipy.special import ive
from plasma_dispersion_audit import AuditHarness,ROOT,MODULES,DARRAY,IARRAY,z_grid,stress_systems,normalization_state


class FaddeevaHarness(AuditHarness):
    def __init__(self,checked=False):
        super().__init__(checked)
        cflags=['-O0','-g'] if checked else ['-O3']
        self._run(['cc','-std=c99','-fPIC',*cflags,'-fno-fast-math','-x','c','-c',
                   str(ROOT/'third_party/faddeeva/Faddeeva.cc'),'-o',str(self.directory/'faddeeva.o')])
        self._run(['cc','-std=c99','-fPIC',*cflags,'-fno-fast-math','-I',str(ROOT/'third_party/faddeeva'),
                   '-c',str(ROOT/'src/faddeeva_bridge.c'),'-o',str(self.directory/'faddeeva_bridge.o')])
        base=['gfortran' ,'-fPIC','-fdefault-real-8','-I',str(self.directory),'-J',str(self.directory)]
        self._run(base+self.flags+['-fno-fast-math','-c',str(ROOT/'src/plasma_z_faddeeva.f90')])
        self._run(base+self.flags+['-fno-fast-math','-c',str(ROOT/'tests/faddeeva_probe.f90')])
        library=self.directory/'native.so'
        objects=[self.directory/f'{m}.o' for m in MODULES]
        # Later production integration includes the standalone objects already.
        for m in ['variant','plasma_dispersion_probe','faddeeva','faddeeva_bridge','plasma_z_faddeeva','faddeeva_probe']:
            obj=self.directory/f'{m}.o'
            if obj not in objects:objects.append(obj)
        self._run(['gfortran','-shared',*map(str,objects),'-lm','-o',str(library)])
        self.native=ct.CDLL(str(library))
        self.native.probe_faddeeva.argtypes=[ct.c_double]*3+[ct.c_int,DARRAY,IARRAY]
        self.native.probe_faddeeva.restype=None
        self.native.native_tensor.argtypes=[DARRAY,DARRAY,ct.c_int,ct.c_int,DARRAY,ct.c_int,
                                            DARRAY,DARRAY,DARRAY,IARRAY,DARRAY,DARRAY]
        self.native.native_tensor.restype=None

    def accurate_z(self,z,k=1.,repeats=1):
        out=np.zeros(3);info=np.zeros(3,np.int32)
        self.native.probe_faddeeva(k,z.real,z.imag,repeats,out,info)
        return {'value':complex(*out[:2]),'seconds':float(out[2]),'status':int(info[0]),
                'c_complex_bits':int(info[1]),'plume_complex_bits':int(info[2])}

    def accurate_tensor(self,s,w,kpar=.31,kperp=3.,va=.02,beta=1.,alpha_scale=1.,n=1024,
                        collision=0.,repeats=1,gamma=None):
        species=np.asarray([*s,collision],dtype=np.float64)
        parameters=np.asarray([w.real,w.imag,kpar,kperp,va,beta,alpha_scale],dtype=np.float64)
        if gamma is None:gamma=ive(np.arange(n+2),kperp**2*s[3]**2*s[2]*s[0]/(2*s[1]))
        gamma=np.ascontiguousarray(gamma)
        output=np.empty(18);arguments=np.empty((2*n+1,2));values=np.empty_like(arguments)
        low=np.empty(54);info=np.zeros(6,np.int32);seconds=np.zeros(1)
        self.native.native_tensor(species,parameters,n,int(collision!=0),gamma,repeats,
                                 output,arguments,values,info,seconds,low)
        return {'chi':(output[::2]+1j*output[1::2]).reshape(3,3),
                'z':values[:,0]+1j*values[:,1], 'arguments':arguments[:,0]+1j*arguments[:,1],
                'seconds':float(seconds[0]),'flags':info[:4].tolist()}


def audit_arguments(h):
    rows=list(z_grid());sampled={}
    for _,system in stress_systems():
        beta,alpha=normalization_state(system)
        for sign in [-1,1]:
            for kp in [.07,3.,30.]:
                for w in [.73+.13j,.73-.08j]:
                    for s in system:
                        if s[0]==0:continue
                        result=h.tensor(s,w,kpar=sign*.31,kperp=kp,beta=beta,alpha_scale=alpha,n=1024)
                        for z in result['arguments']:sampled.setdefault((sign*.31,z),None)
    rows.extend(('recovered_stress',np.sign(k)*z) for k,z in sampled)
    return rows
