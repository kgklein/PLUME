from faddeeva_build import compile_c,c_objects
"""Isolated J1.B.1 harness; never edits or replaces a production implementation.

A temporary instrumented copy changes only the local stopping tolerance and
records branch/termination metadata. The ordinary method calls disprels::zet_in
from an unchanged source build. Species algebra includes the production .inc
files, with explicitly supplied identical Bessel coefficients for every Z method.
"""
from pathlib import Path
import ctypes as ct
import hashlib
import json
import subprocess
import tempfile
import time
import numpy as np
from scipy.special import ive,wofz
ROOT=Path(__file__).resolve().parents[1]
MODULES=('nrtype','nrutil_trim','vars','plasma_z_faddeeva','plasma_z_legacy','plasma_z_backend','cold_susceptibility','auxiliary_normalization',
         'kinetic_core','harmonic_selection','bessel_cache','adaptive_kinetic','cold_reference','electromagnetic_fields','scan_support',
         'functions','bessel','disprels')
STRICT={'electromagnetic_fields','scan_support','plasma_z_faddeeva','plasma_z_legacy','plasma_z_backend','cold_susceptibility','auxiliary_normalization','kinetic_core','harmonic_selection',
        'bessel_cache','adaptive_kinetic'}
TOLERANCES=(1e-7,1e-10,1e-13,1e-15)
CALLBACK=ct.CFUNCTYPE(None,ct.c_double,ct.c_double,ct.c_double,
                    ct.POINTER(ct.c_double),ct.POINTER(ct.c_double))
DARRAY=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
IARRAY=np.ctypeslib.ndpointer(dtype=np.int32,flags='C_CONTIGUOUS')


def variant_source():
    source=(ROOT/'src/plasma_z_legacy.f90').read_text().replace('pure complex function','complex function')
    start=source.index('   complex function zetout_legacy(zin)')
    end=source.index('   end function zetout_legacy',start)+len('   end function zetout_legacy')
    body=source[start:end]
    assert body.count('error=1.e-7')==1
    body=body.replace('zetout_legacy','variant_out')
    body=body.replace('     error=1.e-7','''     error=audit_tolerance
     audit_branch=0;audit_iterations=0;audit_stop=0
     variant_out=cmplx(ieee_value(0.,ieee_quiet_nan),ieee_value(0.,ieee_quiet_nan))''')
    body=body.replace('     fail=.true.','     audit_stop=1\n     fail=.true.')
    body=body.replace('     realte=-2.*x','     audit_branch=1\n     realte=-2.*x')
    body=body.replace('100  realse=realte','100  audit_iterations=audit_iterations+1\n     realse=realte')
    body=body.replace('200  fn=5.','200  audit_branch=2\n     fn=5.')
    body=body.replace('201  terme=term','201  audit_iterations=audit_iterations+1\n     terme=term')
    body=body.replace('     if (abs(real(zetaoz))+abs(aimag(zetaoz)) > 1.) go to 250',
                      '''     if (abs(real(zetaoz))+abs(aimag(zetaoz)) > 1.) then
        audit_stop=2
        go to 250
     endif''')
    body=body.replace('300  if (y < 0.)','300  audit_branch=3\n     if (y < 0.)')
    body=body.replace('301  temp1=bn1*term+an1*terme',
                      '301  audit_iterations=audit_iterations+1\n     temp1=bn1*term+an1*terme')
    body=body.replace('     if (n < 30) go to 301','     if (n < 30) go to 301\n     audit_stop=3')
    return '''! Generated only in a temporary audit build.
module audit_zet_variant
  use, intrinsic :: ieee_arithmetic
  implicit none
  real :: audit_tolerance=1.e-7
  integer :: audit_branch=0,audit_iterations=0,audit_stop=0
contains
  complex function variant_in(k,z) result(v)
    real :: k
    complex :: z
    if(k>0.) then
       v=variant_out(z)
    else if(k<0.) then
       v=-variant_out(-z)
    else
       v=cmplx(ieee_value(0.,ieee_quiet_nan),ieee_value(0.,ieee_quiet_nan))
       audit_branch=0;audit_iterations=0;audit_stop=4
    endif
  end function
'''+body+'\nend module\n'


def scipy_z(k,z):
    return np.sign(k)*1j*np.sqrt(np.pi)*wofz(np.sign(k)*z)


def metrics(actual,reference):
    actual=np.asarray(actual);reference=np.asarray(reference)
    absolute=abs(actual-reference)
    relative=np.divide(absolute,abs(reference),out=np.full(absolute.shape,np.nan),where=abs(reference)>1e-280)
    scaled=absolute/np.maximum(1.,abs(reference))
    return {'absolute':float(np.max(absolute)),'relative':float(np.nanmax(relative)),
            'scaled':float(np.max(scaled))}


class AuditHarness:
    def __init__(self,checked=False):
        self.temp=tempfile.TemporaryDirectory(prefix='plume-z-audit-')
        self.directory=Path(self.temp.name);self.checked=checked
        self.flags=['-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow'] if checked else [
            '-O3','-funroll-loops','-ffast-math']
        base=['gfortran','-fPIC','-fdefault-real-8','-I',str(self.directory),'-I',str(ROOT/'src'),
              '-J',str(self.directory)]
        for module in MODULES:
            flags=self.flags+(['-fno-fast-math'] if module in STRICT else [])
            self._run(base+flags+['-c',str(ROOT/'src'/f'{module}.f90')])
        (self.directory/'variant.f90').write_text(variant_source())
        self._run(base+self.flags+['-c',str(self.directory/'variant.f90')])
        self._run(base+self.flags+['-fno-fast-math','-c',str(ROOT/'tests/plasma_dispersion_probe.f90')])
        compile_c(self.directory,pic=True)
        library=self.directory/'audit.so'
        self._run(['gfortran','-shared',*[str(self.directory/f'{m}.o') for m in MODULES],
                   str(self.directory/'variant.o'),str(self.directory/'plasma_dispersion_probe.o'),*c_objects(self.directory),'-o',str(library)])
        self.library=ct.CDLL(str(library))
        self.library.audit_z.argtypes=[ct.c_double]*4+[ct.c_int,ct.c_int,DARRAY,IARRAY]
        self.library.audit_z.restype=None
        self.library.audit_tensor.argtypes=[DARRAY,DARRAY,ct.c_int,ct.c_int,ct.c_int,ct.c_double,
                                           DARRAY,ct.c_void_p,ct.c_int,DARRAY,DARRAY,DARRAY,IARRAY,DARRAY,DARRAY]
        self.library.audit_tensor.restype=None
        self.callback=CALLBACK(self._reference)

    def _run(self,command):
        result=subprocess.run(command,cwd=self.directory,text=True,capture_output=True)
        if result.returncode: raise RuntimeError(result.stderr)

    def close(self): self.temp.cleanup()

    @staticmethod
    def _reference(k,x,y,vr,vi):
        v=scipy_z(k,complex(x,y));vr[0]=v.real;vi[0]=v.imag

    def z(self,z,k=1.,tolerance=1e-7,method=0,repeats=1):
        value=np.zeros(3);info=np.zeros(8,np.int32)
        self.library.audit_z(k,z.real,z.imag,tolerance,method,repeats,value,info)
        return {'value':complex(*value[:2]),'seconds':value[2],'branch':int(info[0]),
                'iterations':int(info[1]),'stop':int(info[2]),'nonfinite':bool(info[3]),
                'flags':info[4:].tolist()}

    def tensor(self,s,w,kpar=.31,kperp=3.,va=.02,beta=1.,alpha_scale=1.,n=1024,
               tolerance=1e-7,method=0,collision=0.,repeats=1,gamma=None,collisional=None):
        species=np.asarray([*s,collision],dtype=np.float64)
        parameters=np.asarray([w.real,w.imag,kpar,kperp,va,beta,alpha_scale],dtype=np.float64)
        if gamma is None:
            argument=kperp**2*s[3]**2*s[2]*s[0]/(2*s[1])
            gamma=ive(np.arange(n+2),argument)
        gamma=np.ascontiguousarray(gamma,dtype=np.float64)
        output=np.empty(18);arguments=np.empty((2*n+1,2));values=np.empty_like(arguments)
        info=np.zeros(6,np.int32);seconds=np.zeros(1);low=np.empty(54)
        self.library.audit_tensor(species,parameters,n,int(collision!=0 if collisional is None else collisional),method,tolerance,gamma,
                                  ct.cast(self.callback,ct.c_void_p),repeats,output,arguments,values,info,seconds,low)
        return {'chi':(output[::2]+1j*output[1::2]).reshape(3,3),
                'low':(low[::2]+1j*low[1::2]).reshape(3,3,3),
                'arguments':arguments[:,0]+1j*arguments[:,1],
                'z':values[:,0]+1j*values[:,1],'flags':info[:4].tolist(),
                'calls':int(info[4]),'nonfinite':bool(info[5]),'seconds':float(seconds[0])}


def stix_terms(s,w,kpar,kperp,va,n,gamma=None):
    """Independent affine-in-Z Stix terms, C_n + B_n*Z_n, plus drift.

    Derived from the existing test-only mixed_reference formulation (whose NHDS
    notice is in NHDS_LICENSE.txt). Keeping C and B separate measures cancellation
    before and after each harmonic is evaluated. No production equations copied.
    """
    tau,mu,alpha,q,density,u=s
    orders=np.arange(-n,n+1);v=np.sqrt(mu*tau);g=mu/q;d=w-kpar*u
    lam=kperp**2*q*q*alpha*tau/(2*mu)
    if gamma is None: gamma=ive(np.arange(n+2),lam)
    gamma=np.asarray(gamma)
    gn=gamma[abs(orders)]
    hp=gn-.5*(gamma[abs(orders-1)]+gamma[abs(orders+1)])
    c_a=np.full(2*n+1,alpha-1.)
    b_a=(alpha*(d-orders*g)+orders*g)/(kpar*v)
    c_b=(alpha*(w-orders*g)-kpar*u+orders*g)/kpar
    b_b=(w-orders*g)*b_a/kpar
    c=np.zeros((3,3,2*n+1),complex);b=np.zeros_like(c)
    for matrix,a1,b1 in [(c,c_a,c_b),(b,b_a,b_b)]:
        matrix[0,0]=orders**2*gn/lam*a1
        matrix[0,1]=-1j*orders*hp*a1
        matrix[0,2]=kperp*orders*gn/(g*lam)*b1
        matrix[1,1]=(orders**2*gn/lam+2*lam*hp)*a1
        matrix[1,2]=1j*kperp*hp/g*b1
        matrix[2,2]=2*(w-orders*g)*gn/(kpar*v*v*alpha)*b1
        matrix[1,0]=-matrix[0,1];matrix[2,0]=matrix[0,2];matrix[2,1]=-matrix[1,2]
    prefactor=density*mu/(q*q*va*va*w*w)
    extra=np.zeros((3,3),complex)
    extra[2,2]=prefactor*2*w*u/(kpar*v*v*alpha)
    return c*prefactor,b*prefactor,extra


def high_precision_stix(s,w,kpar,kperp,va,n=128,dps=80,gamma=None,zvalues=None):
    """Optional mpmath finite-harmonic reference, independent C/B evaluation.

    gamma/zvalues, when supplied, are promoted binary64 values: this isolates
    arithmetic conditioning from special-function accuracy. Otherwise evaluate
    both functions at the requested arbitrary precision and physical arguments.
    """
    import mpmath as mp
    with mp.workdps(dps):
        tau,mu,alpha,q,density,u=map(mp.mpf,s)
        w=mp.mpc(w);kp=mp.mpf(kpar);kt=mp.mpf(kperp);a=mp.mpf(va)
        thermal=mp.sqrt(mu*tau);g=mu/q;d=w-kp*u
        lam=kt*kt*q*q*alpha*tau/(2*mu)
        coefficients=([mp.exp(-lam)*mp.besseli(i,lam) for i in range(n+2)] if gamma is None
                      else [mp.mpf(v) for v in gamma])
        total=mp.matrix(3,3)
        for j in range(-n,n+1):
            gn=coefficients[abs(j)];hp=gn-(coefficients[abs(j-1)]+coefficients[abs(j+1)])/2
            z=(d-j*g)/(kp*thermal)
            zz=(mp.sign(kp)*mp.j*mp.sqrt(mp.pi)*mp.exp(-z*z)*mp.erfc(-mp.j*mp.sign(kp)*z)
                if zvalues is None else mp.mpc(zvalues[j+n]))
            fa=alpha-1+(alpha*(d-j*g)+j*g)*zz/(kp*thermal)
            fb=(alpha*(w-j*g)-kp*u+j*g)/kp+(w-j*g)*(alpha*(d-j*g)+j*g)*zz/(kp*kp*thermal)
            total[0,0]+=j*j*gn/lam*fa
            total[0,1]+=-mp.j*j*hp*fa
            total[0,2]+=kt*j*gn/(g*lam)*fb
            total[1,1]+=(j*j*gn/lam+2*lam*hp)*fa
            total[1,2]+=mp.j*kt*hp/g*fb
            total[2,2]+=2*(w-j*g)*gn/(kp*thermal*thermal*alpha)*fb
        total[2,2]+=2*w*u/(kp*thermal*thermal*alpha)
        total*=density*mu/(q*q*a*a*w*w)
        total[1,0]=-total[0,1];total[2,0]=total[0,2];total[2,1]=-total[1,2]
        return np.array([[complex(total[i,j]) for j in range(3)] for i in range(3)])


def z_grid():
    from scipy.optimize import root
    rows=[]
    for x in [-20.,-5.,-3.,-2.,-.1,0.,.1,2.,3.,5.,20.]:
        for y in [-3.,-1.,-.1,-1e-12,0.,1e-12,.1,1.,3.]:
            z=complex(x,y)
            regime=('nearly_real' if abs(y)<1e-8 else 'small_series' if abs(z)<=3 and abs(y)<1
                    else 'continued_fraction' if abs(z)<=3 else 'asymptotic')
            rows.append((regime,z))
    for radius in [3.-1e-9,3.+1e-9]:
        for angle in np.linspace(-np.pi,np.pi,25): rows.append(('radius_3_boundary',radius*np.exp(1j*angle)))
    for x in [3.000000001,3.01,3.1,4.]:
        for y in [-1e-12,0.,1e-12]:rows.append(('nearly_real_asymptotic',complex(x,y)))
    for y in [1.-1e-10,1.+1e-10,-1.-1e-10,-1.+1e-10]:
        for x in [0.,1.,2.]: rows.append(('imaginary_1_boundary',complex(x,y)))
    for x in [-1e6,-1e3,-100.,-13.2,13.2,100.,1e3,1e6]:
        for y in [-1.,0.,1.]:rows.append(('large_arguments',complex(x,y)))
    for z in [-13j,-14j,-26j,-27j,-30j,14-14j,200-200j,20j,1000j]:
        rows.append(('domain_and_range',z))
    zero=root(lambda q:[wofz(complex(*q)).real,wofz(complex(*q)).imag],[2.,-1.3],tol=1e-12)
    if not zero.success: raise RuntimeError(zero.message)
    for offset in [0.,1e-12,1e-8,1e-5]: rows.append(('near_zero',complex(*zero.x)+offset))
    return rows


def stress_systems():
    p=(.4,1.,1.7,1.,1.,0.);e=(.1,1836.,.8,-1.,1.1,-.2);a=(2.,.25,1.3,.5,.05,.4)
    return [('warm_pe',[p,e]),('warm_pe_alpha',[p,e,a]),
            ('warm_p_cold_e',[p,(0.,1836.,0.,-1.,1.,-.1)]),
            ('cold_p_warm_e_alpha',[(0.,1.,0.,1.,1.,0.),e,a]),
            ('cold_pe_warm_alpha',[(0.,1.,0.,1.,1.,0.),(0.,1836.,0.,-1.,1.1,0.),a])]


def normalization_state(species):
    if species[0][0]>0:return species[0][0],species[0][2]
    return next(s[0] for s in species if s[0]>0),1.


def json_clean(obj):
    if isinstance(obj,complex):return [obj.real,obj.imag]
    if isinstance(obj,np.ndarray):return json_clean(obj.tolist())
    if isinstance(obj,dict):return {str(k):json_clean(v) for k,v in obj.items()}
    if isinstance(obj,(list,tuple)):return [json_clean(v) for v in obj]
    if isinstance(obj,(np.integer,np.bool_)):return obj.item()
    if isinstance(obj,(float,np.floating)):return float(obj) if np.isfinite(obj) else str(obj)
    return obj
