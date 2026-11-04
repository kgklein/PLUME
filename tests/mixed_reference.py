"""Independent test-only NHDS/Stix bi-Maxwellian harmonic sum in d_ref units.

The Y-tensor formulas follow NHDS src/suscep.f90 (calc_ypsilon), implementing
Stix (1992), equation 10-57. This vectorized Python adaptation uses SciPy's
independent Bessel/Faddeeva implementations, not PLUME's kernels. NHDS's
BSD notice is reproduced in NHDS_LICENSE.txt. Cold response uses the independent
fluid linear solve from Stage 1. No production Fortran is invoked.
"""

import numpy as np
from scipy.special import ive, wofz
from test_cold_susceptibility import Case, fluid_tensor


def kinetic_tensor(s, omega, kpar, kperp, va, harmonics=60):
    tau, mu, alpha, charge, density, drift = s
    thermal = np.sqrt(mu*tau)
    gyro = mu/charge
    z = 0.5*(kperp*thermal/gyro)**2*alpha
    n = np.arange(-harmonics, harmonics+1)
    gamma = ive(abs(n), z)
    derivative = (ive(abs(n-1), z)+ive(abs(n+1), z))/2
    # Exact Bessel limits: n^2*Gamma_n/z is 1/2 for n=+/-1, zero otherwise.
    transverse_weight = np.where(abs(n) == 1, .5, 0.) if z == 0 else n*n*gamma/z
    resonance = omega-kpar*drift-n*gyro
    argument = resonance/(kpar*thermal)
    # Causal continuation for either sign of parallel wavevector.
    plasma_z = np.sign(kpar)*1j*np.sqrt(np.pi)*wofz(np.sign(kpar)*argument)
    a = alpha-1 + (alpha*resonance+n*gyro)*plasma_z/(kpar*thermal)
    b = (alpha*(omega-n*gyro)-(kpar*drift-n*gyro))/kpar
    b = b+(omega-n*gyro)*(alpha*resonance+n*gyro)*plasma_z/(kpar**2*thermal)
    response = np.zeros((3, 3), dtype=complex)
    response[0, 0] = np.sum(transverse_weight*a)
    response[0, 1] = np.sum(-1j*n*(gamma-derivative)*a)
    response[0, 2] = 0 if z == 0 else np.sum(kperp*n*gamma*b/(gyro*z))
    response[1, 1] = np.sum((transverse_weight+2*z*(gamma-derivative))*a)
    response[1, 2] = np.sum(1j*kperp*(gamma-derivative)*b/gyro)
    response[2, 2] = np.sum(2*(omega-n*gyro)*gamma*b/(kpar*thermal**2*alpha))
    response[2, 2] += 2*omega*drift/(kpar*thermal**2*alpha)
    response[1, 0] = -response[0, 1]
    response[2, 0] = response[0, 2]
    response[2, 1] = -response[1, 2]
    return response*density*mu/(charge**2*va**2*omega**2)


def species_tensor(s, omega, kpar, kperp, va, harmonics=60):
    tau, mu, alpha, charge, density, drift = s
    if tau > 0:
        return kinetic_tensor(s, omega, kpar, kperp, va, harmonics)
    return fluid_tensor(Case(omega=omega, kpar=kpar, kperp=kperp, va=va,
                            mu=mu, charge=charge, density=density, drift=drift))


def wave_matrix(species, omega, kpar, kperp, va, harmonics=60):
    n = np.array([kperp, 0., kpar])/(va*omega)
    return (np.eye(3, dtype=complex)*(1-n@n)+np.outer(n, n) +
            sum(species_tensor(s, omega, kpar, kperp, va, harmonics) for s in species))
