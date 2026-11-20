"""Audit-only energy, absorption and independent arbitrary-precision derivatives.

Conventions follow the explicit real-frequency formula supplied for the audit.
The warm tensor is independently evaluated from NHDS/Stix (see mixed_reference
and NHDS_LICENSE.txt). Production formulas are not imported or edited.
"""
import numpy as np
try:
    import mpmath as mp
except ImportError:
    mp=None
from scipy.special import wofz
from mixed_reference import species_tensor
from plasma_dispersion_audit import stix_terms


def hermitian(a):return (a+np.swapaxes(a.conj(),-1,-2))/2

def antihermitian(a):return (a-np.swapaxes(a.conj(),-1,-2))/(2j)

def contraction(e,a):return np.einsum('i,...ij,j->...',e.conj(),a,e).real

def absorption(e,chi):return contraction(e,antihermitian(chi))


def channels(e,low):
    """Historical seven-channel semantics, with explicit signed harmonic axis."""
    zero=low[1];y=e[1];z=e[2];perp=np.array([e[0],e[1],0j])
    return np.array([abs(y)**2*zero[1,1].imag,(y.conjugate()*z*zero[1,2]).imag,
                     (z.conjugate()*y*zero[2,1]).imag,abs(z)**2*zero[2,2].imag,
                     absorption(e,zero),absorption(perp,low[2]),absorption(perp,low[0])])


def harmonic_reference(s,w,kpar,kperp,a,n=64):
    c,b,extra=stix_terms(s,w,kpar,kperp,a,n)
    orders=np.arange(-n,n+1);z=(w-kpar*s[5]-orders*s[1]/s[3])/(kpar*np.sqrt(s[0]*s[1]))
    zz=np.sign(kpar)*1j*np.sqrt(np.pi)*wofz(np.sign(kpar)*z)
    terms=np.moveaxis(c+b*zz,-1,0)
    return terms,extra


def mp_warm(s,w,kpar,kperp,a,n=32,dps=65):
    """chi and d(w*chi)/dw using Z'=-2(1+z Z), all arithmetic at dps."""
    with mp.workdps(dps):
        tau,mu,alpha,q,density,u=map(mp.mpf,s)
        w=mp.mpc(w);kp=mp.mpf(kpar);kt=mp.mpf(kperp);a=mp.mpf(a)
        thermal=mp.sqrt(mu*tau);g=mu/q;delta=w-kp*u
        lam=kt*kt*q*q*alpha*tau/(2*mu)
        gamma=[mp.exp(-lam)*mp.besseli(i,lam) for i in range(n+2)]
        matrix=mp.matrix(3);derivative=mp.matrix(3)
        for j in range(-n,n+1):
            m=j*g;gn=gamma[abs(j)];hp=gn-(gamma[abs(j-1)]+gamma[abs(j+1)])/2
            z=(delta-m)/(kp*thermal)
            zz=mp.sign(kp)*mp.j*mp.sqrt(mp.pi)*mp.exp(-z*z)*mp.erfc(-mp.j*mp.sign(kp)*z)
            zp=-2*(1+z*zz)
            h=alpha*delta+(1-alpha)*m
            fa=alpha-1+h*zz/(kp*thermal)
            da=alpha*zz/(kp*thermal)+h*zp/(kp*thermal)**2
            fb=(alpha*(w-m)-kp*u+m)/kp+(w-m)*h*zz/(kp*kp*thermal)
            db=alpha/kp+((h+(w-m)*alpha)*zz+(w-m)*h*zp/(kp*thermal))/(kp*kp*thermal)
            transverse=(mp.mpf('.5') if abs(j)==1 else 0) if lam==0 else j*j*gn/lam
            for response,ax,bx in [(matrix,fa,fb),(derivative,da,db)]:
                response[0,0]+=transverse*ax
                response[0,1]+=-mp.j*j*hp*ax
                if lam!=0:response[0,2]+=kt*j*gn*bx/(g*lam)
                response[1,1]+=(transverse+2*lam*hp)*ax
                response[1,2]+=mp.j*kt*hp*bx/g
                response[2,2]+=2*gn*((w-m)*bx)/(kp*thermal*thermal*alpha) if response is matrix else 0
            derivative[2,2]+=2*gn*(fb+(w-m)*db)/(kp*thermal*thermal*alpha)
        matrix[2,2]+=2*w*u/(kp*thermal*thermal*alpha)
        derivative[2,2]+=2*u/(kp*thermal*thermal*alpha)
        for response in [matrix,derivative]:
            response[1,0]=-response[0,1];response[2,0]=response[0,2];response[2,1]=-response[1,2]
        prefactor=density*mu/(q*q*a*a*w*w)
        result=matrix*prefactor;energy=(w*derivative-matrix)*prefactor
        return (np.array([[complex(result[i,j]) for j in range(3)] for i in range(3)]),
                np.array([[complex(energy[i,j]) for j in range(3)] for i in range(3)]))


def mp_cold(s,w,kpar,kperp,a,dps=65):
    with mp.workdps(dps):
        _,mu,_,q,d,u=s;mu,q,d,u=map(mp.mpf,[mu,q,d,u]);kp,kt,a=map(mp.mpf,[kpar,kperp,a]);g=mu/q
        plasma=d*mu/(q*q*a*a)
        def chi(w):
            delta=w-kp*u;b=kt*u;den=(delta-g)*(delta+g)
            r=-plasma*delta/(w*w*(delta+g));l=-plasma*delta/(w*w*(delta-g))
            j=-plasma*b*delta/(w*w*den);m=mp.j*plasma*b*g/(w*w*den)
            return mp.matrix([[(r+l)/2,-mp.j*(r-l)/2,j],[mp.j*(r-l)/2,(r+l)/2,m],
                              [j,-m,-plasma/(delta*delta)-plasma*b*b/(w*w*den)]])
        point=mp.mpc(w);matrix=chi(point)
        derivative=mp.matrix(3)
        for i in range(3):
            for j in range(3):derivative[i,j]=mp.diff(lambda z:z*chi(z)[i,j],point)
        return (np.array([[complex(matrix[i,j]) for j in range(3)] for i in range(3)]),
                np.array([[complex(derivative[i,j]) for j in range(3)] for i in range(3)]))


def cold_linear_derivative(s,w,kpar,kperp,a):
    """Independent derivative of the fluid linear solve, without finite steps."""
    _,mu,_,q,d,u=s;g=mu/q;delta=w-kpar*u
    matrix=np.array([[-1j*delta,-g,0],[g,-1j*delta,0],[0,0,-1j*delta]],complex)
    chi=np.zeros((3,3),complex);derivative=chi.copy()
    for col in range(3):
        e=np.eye(3,dtype=complex)[col];b=np.cross([kperp,0,kpar],e)/(a*w)
        magnetic_force=a*u*np.cross([0,0,1],b)
        v=np.linalg.solve(matrix,g*(e+magnetic_force))
        dv=np.linalg.solve(matrix,-g*magnetic_force/w+1j*v)
        eta=np.dot([kperp,0,kpar],v)/(a*delta)
        deta=np.dot([kperp,0,kpar],dv)/(a*delta)-eta/delta
        chi[:,col]=1j*d/(q*a*a*w)*(v+np.array([0,0,a*u*eta]))
        derivative[:,col]=1j*d/(q*a*a)*(dv+np.array([0,0,a*u*deta]))
    return chi,derivative


def energy_reference(species,w,kpar,kperp,a,e,n=32,dps=65):
    if mp is None:
        # Independent fourth-order real-frequency derivative when the optional
        # audit-only arbitrary-precision package is not importable.
        results=[]
        for s in species:
            if s[0]==0:
                results.append(cold_linear_derivative(s,w,kpar,kperp,a));continue
            evaluate=lambda x:species_tensor(s,x,kpar,kperp,a,max(n,128))
            h=w*5e-4
            f=lambda x:x*evaluate(x)
            derivative=(f(w-2*h)-8*f(w-h)+8*f(w+h)-f(w+2*h))/(12*h)
            results.append((evaluate(w),derivative))
    else:
        results=[mp_warm(s,w,kpar,kperp,a,n,dps) if s[0]>0 else mp_cold(s,w,kpar,kperp,a,dps) for s in species]
    chi=np.array([r[0] for r in results]);derivative=hermitian(np.array([r[1] for r in results]))
    magnetic=np.cross([kperp,0,kpar],e)/(a*w)
    contributions=contraction(e,derivative)
    vacuum=float(np.vdot(e,e).real);mag=float(np.vdot(magnetic,magnetic).real)
    return dict(chi=chi,derivative=derivative,contributions=contributions,vacuum=vacuum,magnetic=mag,
                energy=vacuum+mag+float(sum(contributions)),q=absorption(e,chi))


def fd_derivative(evaluate,w,relative,centered=True):
    h=w*relative
    if centered:return hermitian(((w+h)*evaluate(w+h)-(w-h)*evaluate(w-h))/(2*h))
    return hermitian(((1+relative)*evaluate(w+h)-evaluate(w))/relative)
