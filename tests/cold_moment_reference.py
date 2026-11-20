"""Independent normalized fluid solve and direct bi-Maxwellian Vlasov quadrature.

No PLUME formulas/includes are used for these moment oracles. The gyroangle
collocation solves linearized Vlasov, then integrates velocity and density.
Its real-velocity contour is valid for Im(w)>0 only; damped cases separately
use the established Landau-continued independent NHDS/Stix susceptibility.
"""
import numpy as np
from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.laguerre import laggauss


def fluid_moments(s,w,kpar,kperp,a,e,b):
    _,mu,_,q,_,u=s;g=mu/q;delta=w-kpar*u
    matrix=np.array([[-1j*delta,-g,0],[g,-1j*delta,0],[0,0,-1j*delta]],complex)
    force=g*(e+a*u*np.cross([0.,0.,1.],b))
    v=np.linalg.solve(matrix,force)
    eta=np.dot([kperp,0.,kpar],v)/(a*delta)
    return v,eta


def susceptibility_moments(s,w,kpar,kperp,a,e,chi):
    # Independent coupled current+continuity 4x4 solve, without the production
    # algebraic elimination; chi comes from the independent NHDS/Stix oracle.
    _,_,_,q,d,u=s
    matrix=np.eye(4,dtype=complex);matrix[2,3]=a*u
    matrix[3,:3]=-np.array([kperp,0.,kpar])/a
    matrix[3,3]=w-kpar*u
    rhs=np.r_[-1j*w*a*a*q/d*(chi@e),0j]
    answer=np.linalg.solve(matrix,rhs)
    return answer[:3],answer[3]


def vlasov_quadrature(s,w,kpar,kperp,a,e,b,order=24,angles=32):
    if w.imag<=0:raise ValueError('Real-contour Vlasov quadrature requires Im(w)>0')
    tau,mu,alpha,q,_,u=s;tpar=np.sqrt(mu*tau);tperp=tpar*np.sqrt(alpha);g=mu/q
    x,wx=hermgauss(order);r,wr=laggauss(order)
    phi=2*np.pi*np.arange(angles)/angles
    freq=np.fft.fftfreq(angles,1/angles)
    derivative=np.fft.ifft(1j*freq[:,None]*np.fft.fft(np.eye(angles),axis=0),axis=0)
    velocity=np.empty((order,order,angles,3))
    velocity[:,:,:,0]=tperp*np.sqrt(r)[None,:,None]*np.cos(phi)
    velocity[:,:,:,1]=tperp*np.sqrt(r)[None,:,None]*np.sin(phi)
    velocity[:,:,:,2]=(u+tpar*x)[:,None,None]
    force=e+a*np.cross(velocity,b)
    grad=-2*velocity/np.array([tperp**2,tperp**2,tpar**2])
    grad[:,:,:,2]=-2*(velocity[:,:,:,2]-u)/tpar**2
    rhs=-(g/a)*np.sum(force*grad,axis=-1)
    doppler=w-kpar*velocity[:,:,:,2]-kperp*velocity[:,:,:,0]
    matrix=np.broadcast_to(-g*derivative,(*doppler.shape,angles)).copy()
    indices=np.arange(angles);matrix[...,indices,indices]+=-1j*doppler
    perturbation=np.linalg.solve(matrix,rhs[...,None])[...,0]
    weight=wx[:,None,None]*wr[None,:,None]/(np.sqrt(np.pi)*angles)
    eta=np.sum(weight*perturbation)
    flux=np.sum(weight[...,None]*velocity*perturbation[...,None],axis=(0,1,2))
    v=a*(flux-np.array([0,0,u*eta]))
    return v,eta
