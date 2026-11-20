"""Direct cold-kernel tests: paper equations and an independent fluid solve.

No production dispersion routine is linked into either probe build.
See COLD_SUSCEPTIBILITY_DESIGN.md for normalization, signs and tolerances.
"""

from dataclasses import dataclass, replace
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class Case:
    omega: complex = 0.73 + 0.19j
    kpar: float = 0.31
    kperp: float = 0.47
    va: float = 0.02
    mu: float = 1.0
    charge: float = 1.0
    density: float = 1.0
    drift: float = 0.0
    tau: float = 0.0
    alpha: float = 0.0
    collision: float = 0.0

    def record(self):
        values = (self.omega.real, self.omega.imag, self.kpar, self.kperp,
                  self.va, self.mu, self.charge, self.density, self.drift,
                  self.tau, self.alpha, self.collision)
        return ' '.join(format(x, '.17g') for x in values)


def paper_tensor(case):
    """Literal equations (5)-(9), independent of the kernel's combined S/H."""
    w = case.omega
    shifted = w - case.kpar * case.drift
    gyro = case.mu / case.charge
    wp2 = case.density * case.mu / (case.charge**2 * case.va**2)
    right = -wp2 / w**2 * shifted / (shifted + gyro)
    left = -wp2 / w**2 * shifted / (shifted - gyro)
    j = -wp2 / w**2 * case.kperp * case.drift * shifted / (shifted**2 - gyro**2)
    m = 1j * wp2 / w**2 * case.kperp * case.drift * gyro / (shifted**2 - gyro**2)
    p = -wp2 / w**2 * (w**2 / shifted**2 +
                       case.kperp**2 * case.drift**2 / (shifted**2 - gyro**2))
    s, h = (right + left) / 2, (right - left) / 2
    return np.array([[s, -1j*h, j], [1j*h, s, m], [j, -m, p]])


def fluid_tensor(case):
    """Solve momentum, Faraday and continuity equations for each E basis.

    h is velocity divided by q_s/(m_s*Omega_R). The force includes U0 x B1;
    the current includes the perturbed density times U0. Neither R/L nor
    the kernel's rational tensor expressions are used here.
    """
    k = np.array([case.kperp, 0.0, case.kpar])
    u = np.array([0.0, 0.0, case.drift])
    shifted = case.omega - k @ u
    gyro = case.mu / case.charge
    wp2 = case.density * case.mu / (case.charge**2 * case.va**2)
    momentum = np.array([[-1j*shifted, -gyro, 0],
                         [gyro, -1j*shifted, 0],
                         [0, 0, -1j*shifted]], dtype=complex)
    response = np.empty((3, 3), dtype=complex)
    for column, electric in enumerate(np.eye(3)):
        magnetic = np.cross(k, electric) / (case.va * case.omega)
        force = electric + np.cross(case.va * u, magnetic)
        h = np.linalg.solve(momentum, force)
        density_current = u * (k @ h) / shifted
        response[:, column] = 1j * wp2 / case.omega * (h + density_current)
    return response


class ColdSusceptibilityTests(unittest.TestCase):
    max_scaled_error = {'paper': 0.0, 'fluid': 0.0, 'near_pole': 0.0}

    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.temporary.cleanup)
        build = Path(cls.temporary.name)
        subprocess.run(['make', '-j2'], cwd=ROOT, check=True, capture_output=True)
        cls.probes = {}
        production = build / 'production_probe'
        subprocess.run(['gfortran', '-O3', '-DDOUBLE', '-fdefault-real-8',
                        '-funroll-loops', '-ffast-math', '-I', str(ROOT / 'include'),
                        str(ROOT / 'tests/cold_susceptibility_probe.f90'),
                        str(ROOT / 'cold_susceptibility.o'), str(ROOT / 'vars.o'),
                        str(ROOT / 'nrtype.o'), '-o', str(production)],
                       cwd=build, check=True, capture_output=True)
        cls.probes['production'] = production
        checked = build / 'checked_probe'
        subprocess.run(['gfortran', '-O0', '-g', '-fdefault-real-8',
                        '-fcheck=all', '-ffpe-trap=invalid,zero,overflow',
                        '-I', str(build), '-J', str(build),
                        str(ROOT / 'src/nrtype.f90'), str(ROOT / 'src/vars.f90'),
                        str(ROOT / 'src/cold_susceptibility.f90'),
                        str(ROOT / 'tests/cold_susceptibility_probe.f90'),
                        '-o', str(checked)], cwd=build, check=True, capture_output=True)
        cls.probes['checked'] = checked

    @classmethod
    def tearDownClass(cls):
        print('Cold tensor max |actual-reference|/max(1,|reference|): ' +
              ', '.join(f'{key}={value:.3e}' for key, value in cls.max_scaled_error.items()))

    def evaluate(self, build, cases):
        result = subprocess.run([str(self.probes[build])],
                                input='\n'.join(c.record() for c in cases) + '\n',
                                text=True, capture_output=True, check=True, timeout=20)
        lines = result.stdout.splitlines()
        self.assertEqual(lines[0].split(), ['PRECISION', '15', '64'])
        self.assertEqual(len(lines), 1 + 10 * len(cases))
        responses = []
        for start in range(1, len(lines), 10):
            self.assertEqual(lines[start].split()[0], 'STATUS')
            status = int(lines[start].split()[1])
            values = [[float(x) for x in line.split()] for line in lines[start+1:start+10]]
            self.assertTrue(all(len(pair) == 2 for pair in values))
            chi = np.array([complex(*pair) for pair in values]).reshape(3, 3)
            responses.append((status, chi))
        return responses

    def compare(self, cases, near=False):
        for build in self.probes:
            actual = self.evaluate(build, cases)
            for case, (status, chi) in zip(cases, actual):
                with self.subTest(build=build, case=case):
                    self.assertEqual(status, 0)
                    self.assertTrue(np.isfinite(chi).all())
                    for name, oracle in [('paper', paper_tensor), ('fluid', fluid_tensor)]:
                        expected = oracle(case)
                        scaled = float(np.max(abs(chi-expected) / np.maximum(1, abs(expected))))
                        metric = 'near_pole' if near else name
                        type(self).max_scaled_error[metric] = max(self.max_scaled_error[metric], scaled)
                        np.testing.assert_allclose(chi, expected, rtol=1e-9 if near else 5e-13,
                                                   atol=5e-13, err_msg=f'{build}: {name}: {case}')
                    # Transpose identities remain exact for complex frequencies.
                    self.assertEqual(chi[0, 0], chi[1, 1])
                    self.assertEqual(chi[1, 0], -chi[0, 1])
                    self.assertEqual(chi[2, 0], chi[0, 2])
                    self.assertEqual(chi[2, 1], -chi[1, 2])

    def test_stationary(self):
        cases = [Case(omega=w, mu=mu, charge=q, kperp=kp)
                 for w in (0.73+0j, 0.73+0.19j, 0.73-0.19j)
                 for mu, q in ((1, 1), (1836, -1), (0.25, 0.5))
                 for kp in (0.0, 0.47)]
        self.compare(cases)
        for build in self.probes:
            for case, (_, chi) in zip(cases, self.evaluate(build, cases)):
                a = case.density*case.mu/(case.charge**2*case.va**2)
                g, w = case.mu/case.charge, case.omega
                np.testing.assert_allclose(chi[0, 0], -a/(w*w-g*g), rtol=5e-13)
                np.testing.assert_allclose(chi[0, 1], -1j*a*g/(w*(w*w-g*g)), rtol=5e-13)
                np.testing.assert_allclose(chi[2, 2], -a/(w*w), rtol=5e-13)
                np.testing.assert_array_equal(chi[[0, 1, 2, 2], [2, 2, 0, 1]], 0)

    def test_parallel_drift(self):
        cases = [Case(kperp=0, drift=u, kpar=k) for u in (-1.1, 0.3, 1.1) for k in (-0.4, 0.4)]
        self.compare(cases)
        for build in self.probes:
            for case, (_, chi) in zip(cases, self.evaluate(build, cases)):
                delta = case.omega - case.kpar*case.drift
                a = case.density*case.mu/(case.charge**2*case.va**2)
                g = case.mu/case.charge
                for handedness in (-1, 1):
                    vector = np.array([1, handedness*1j, 0])
                    eigenvalue = -a*delta/(case.omega**2*(delta+handedness*g))
                    np.testing.assert_allclose(chi @ vector, eigenvalue*vector, rtol=5e-13, atol=5e-13)
                np.testing.assert_array_equal(chi[[0, 1, 2, 2], [2, 2, 0, 1]], 0)

    def test_oblique_drift_and_angles(self):
        cases = [Case(kpar=0.8*np.cos(theta), kperp=sign*0.8*np.sin(theta), drift=u)
                 for theta in np.deg2rad([15, 45, 75, 135])
                 for sign in (-1, 1) for u in (-0.9, 0.9)]
        cases += [Case(kpar=0, kperp=0.8, drift=u) for u in (-0.9, 0.9)]
        self.compare(cases)
        for build in self.probes:
            for case, (_, chi) in zip(cases, self.evaluate(build, cases)):
                self.assertTrue(np.all(abs(chi[[0, 1, 2, 2], [2, 2, 0, 1]]) > 0))
                electrostatic = -case.density*case.mu/(case.charge**2*case.va**2) / (
                    case.omega-case.kpar*case.drift)**2
                self.assertGreater(abs(chi[2, 2]-electrostatic), 0)

    def test_complex_frequency(self):
        cases = [Case(omega=w, drift=0.7) for w in
                 (0.8+0.3j, 0.8-0.3j, -0.8+0.3j, -0.8-0.3j, 0.3j, -0.3j)]
        self.compare(cases)
        for build in self.probes:
            for _, chi in self.evaluate(build, cases):
                self.assertGreater(np.max(abs(chi-chi.conj().T)), 1)
            _, real_chi = self.evaluate(build, [Case(omega=0.8+0j, drift=0.7)])[0]
            np.testing.assert_array_equal(real_chi, real_chi.conj().T)

    def test_species_and_normalization(self):
        self.compare([Case(mu=mu, charge=q, density=d, va=a, drift=u)
                      for mu, q in ((1, 1), (100, -1), (1836, -1), (1836.15267343, -1), (0.25, 0.5))
                      for d in (0.03, 1.7) for a in (1e-4, 0.02, 0.1) for u in (-0.6, 0.6)])

    def test_no_thermal_dependency(self):
        # The probe also sets global beta and vtp to zero, with no species arrays.
        cases = [Case(drift=0.7, tau=t, alpha=a, collision=c)
                 for t, a, c in ((0, 0, 0), (1, 2, 3), (-4, -5, -6),
                                  (float('nan'), float('nan'), float('nan')))]
        for build in self.probes:
            values = self.evaluate(build, cases)
            for status, chi in values:
                self.assertEqual(status, 0)
                np.testing.assert_array_equal(chi, values[0][1])

    def test_charge_density_and_speed_scaling(self):
        case = Case(drift=0.6)
        cases = [case, replace(case, charge=-1), replace(case, density=2), replace(case, va=2*case.va)]
        self.compare(cases)
        for build in self.probes:
            original, reversed_charge, doubled_density, doubled_speed = [chi for _, chi in self.evaluate(build, cases)]
            expected = original.copy()
            expected[[0, 1, 1, 2], [1, 0, 2, 1]] *= -1
            np.testing.assert_array_equal(reversed_charge, expected)
            np.testing.assert_array_equal(doubled_density, 2*original)
            np.testing.assert_array_equal(doubled_speed, original/4)

    def test_seeded_nonsingular_points(self):
        rng = np.random.default_rng(73621)
        cases = [Case(omega=complex(rng.uniform(-3, 3), rng.choice([-1, 1])*rng.uniform(0.2, 0.8)),
                      kpar=rng.uniform(-2, 2), kperp=rng.uniform(-2, 2),
                      drift=rng.uniform(-2, 2), mu=10**rng.uniform(-1, 3),
                      charge=rng.choice([-1, 1])*rng.uniform(0.5, 2),
                      density=rng.uniform(0.01, 2), va=10**rng.uniform(-4, -1)) for _ in range(40)]
        self.compare(cases)

    def test_near_poles(self):
        case = Case(kpar=0.5, kperp=0.4, drift=0.5, va=0.1)
        # Distance 1e-5..1e-3 from poles; allow 1e-9 relative error in the
        # independent, increasingly ill-conditioned momentum solve.
        self.compare([replace(case, omega=complex(pole, sign*distance))
                      for pole in (0, 0.25, 1.25, -0.75)
                      for distance in (1e-3, 1e-4, 1e-5) for sign in (-1, 1)], near=True)

    def test_exact_poles(self):
        cases = [Case(omega=0j), Case(omega=0.25+0j, kpar=0.5, drift=0.5),
                 Case(omega=1.25+0j, kpar=0.5, drift=0.5),
                 Case(omega=-0.75+0j, kpar=0.5, drift=0.5),
                 Case(omega=1+0j, charge=-1), Case(omega=-1+0j, charge=-1)]
        for build in self.probes:
            for expected, (status, chi) in zip([2, 3, 4, 4, 4, 4], self.evaluate(build, cases)):
                self.assertEqual(status, expected)
                self.assertTrue(np.isnan(chi.real).all() and np.isnan(chi.imag).all())

    def test_zero_density(self):
        cases = [Case(density=0, omega=w, kpar=0.5, drift=0.5)
                 for w in (0j, 0.25+0j, 1.25+0j, -0.75+0j, 0.7+0.2j)]
        for build in self.probes:
            for status, chi in self.evaluate(build, cases):
                self.assertEqual(status, 0)
                np.testing.assert_array_equal(chi, np.zeros((3, 3), dtype=complex))

    def test_invalid_inputs(self):
        case = Case()
        cases = [replace(case, **{key: value}) for key, values in
                 [('mu', [0, -1]), ('charge', [0]), ('density', [-1]), ('va', [0, -1])]
                 for value in values]
        for key in ('kpar', 'kperp', 'va', 'mu', 'charge', 'density', 'drift'):
            cases += [replace(case, **{key: value}) for value in (float('nan'), float('inf'), -float('inf'))]
        cases += [replace(case, omega=complex(value, 0)) for value in (float('nan'), float('inf'))]
        cases += [replace(case, omega=complex(1, value)) for value in (float('nan'), float('inf'))]
        cases += [replace(case, density=0, charge=0)]  # validation precedes zero-density shortcut
        for build in self.probes:
            for status, chi in self.evaluate(build, cases):
                self.assertEqual(status, 1)
                self.assertTrue(np.isnan(chi.real).all() and np.isnan(chi.imag).all())


if __name__ == '__main__':
    unittest.main()
