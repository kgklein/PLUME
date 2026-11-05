"""Independent longitudinal response of a drifting bi-Maxwellian with C[f1]=-nu*f1.

Solve the longitudinal kinetic moment integral before pairing gyro harmonics.
The electric drive and v x B1 force use the physical omega; only the resolvent
contains omega+i*nu. No expanded polynomial from the production include is used.
"""
import numpy as np
from scipy.special import ive,roots_hermite
from plasma_dispersion_audit import scipy_z


def coefficients(species,kperp,n):
    tau,mu,alpha,q,_,_=species
    return ive(np.arange(n+2),kperp**2*q*q*alpha*tau/(2*mu))


def zz_moments(species,omega,kpar,kperp,va,nu,n,gamma=None,zvalues=None):
    """Return signed chi_zz,n, regrouping remainder, and the total response.

    M_j=<v_parallel**j/(omega+i*nu-kpar*v_parallel-n*Omega)>.
    The recurrences use only Gaussian normalization and mean drift. Optional
    supplied Gamma/Z values hold all special-function inputs fixed for audits.
    """
    tau,mu,alpha,q,density,u=species
    gamma=coefficients(species,kperp,n) if gamma is None else gamma
    orders=np.arange(-n,n+1)
    thermal=np.sqrt(mu*tau);gyro=mu/q;wp2=density*mu/(q*q*va*va)
    r=np.asarray(omega+1j*nu-orders*gyro,dtype=np.clongdouble)
    zeta=(omega+1j*nu-kpar*u-orders*gyro)/(kpar*thermal)
    zv=scipy_z(kpar,zeta) if zvalues is None else zvalues
    m0=-np.asarray(zv,dtype=np.clongdouble)/(kpar*thermal)
    m1=(r*m0-1.)/kpar
    m2=(r*m1-u)/kpar
    electric_magnetic_v2=1.+orders*gyro/omega*(1./alpha-1.)
    electric_magnetic_v1=(orders*gyro/omega-1.)*u
    response=-2.*wp2/(omega*thermal**2)*np.asarray(gamma)[abs(orders)]*(
        electric_magnetic_v2*m2+electric_magnetic_v1*m1)
    # The production bulk expression adds the complete density/drift term after
    # its truncated Bessel sum. Signed harmonics each contain their own part.
    gamma_sum=np.longdouble(gamma[0])+2.*np.sum(np.asarray(gamma[1:n+1],dtype=np.longdouble))
    drift=2.*wp2*u/(omega*kpar*thermal**2*alpha)
    remainder=drift*(1.-gamma_sum)
    return np.asarray(response,dtype=np.complex128),complex(remainder),complex(np.sum(response)+remainder)


def zz_quadrature(species,omega,kpar,kperp,va,nu,n,gamma=None,points=256):
    """Direct Gaussian quadrature, valid only for Im(omega+i*nu)>0.

    The numerator follows -grad_v f0 dot (E + v x B1/c), with E along z and
    B1=c*k x E/omega. The current includes the full equilibrium-frame velocity.
    Lower-half-plane continuation must instead use zz_moments.
    """
    if (omega+1j*nu).imag<=0:raise ValueError('Quadrature requires Im(omega+i*nu)>0')
    tau,mu,alpha,q,density,u=species
    gamma=coefficients(species,kperp,n) if gamma is None else gamma
    thermal=np.sqrt(mu*tau);gyro=mu/q;wp2=density*mu/(q*q*va*va)
    nodes,weights=roots_hermite(points)
    velocity=u+thermal*nodes
    orders=np.arange(-n,n+1)
    source=(velocity-u)[None,:]+orders[:,None]*gyro/omega*((1./alpha-1.)*velocity[None,:]+u)
    denominator=omega+1j*nu-kpar*velocity[None,:]-orders[:,None]*gyro
    response=-2.*wp2/(omega*thermal**2)*np.asarray(gamma)[abs(orders)]*np.sum(
        weights[None,:]/np.sqrt(np.pi)*velocity[None,:]*source/denominator,axis=1)
    gamma_sum=np.longdouble(gamma[0])+2.*np.sum(np.asarray(gamma[1:n+1],dtype=np.longdouble))
    drift=2.*wp2*u/(omega*kpar*thermal**2*alpha)
    return response,complex(np.sum(response)+drift*(1.-gamma_sum))
