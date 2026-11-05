from faddeeva_build import c_objects
"""J1.B.1 characterization tests; all probes and tolerance changes are isolated.

These assertions document detected defects; they do not demand that the legacy
Z implementation meets reference accuracy it does not currently provide.
"""
import hashlib
from pathlib import Path
import subprocess
import unittest
import numpy as np
from scipy.special import ive
from plasma_dispersion_audit import (AuditHarness,ROOT,MODULES,TOLERANCES,scipy_z,
                                     stix_terms,metrics,z_grid)
from mixed_reference import species_tensor

S=(2.,.25,1.3,.5,.05,.4)
W=.73+.13j

class PlasmaDispersionAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(['make','-j2'],cwd=ROOT,check=True,capture_output=True)
        cls.hashes={p:hashlib.sha256(p.read_bytes()).digest() for p in (ROOT/'src').glob('*') if p.is_file()}
        cls.harnesses={False:AuditHarness(),True:AuditHarness(checked=True)}
        for h in cls.harnesses.values():cls.addClassCleanup(h.close)
        cls.direct=cls.harnesses[False].directory/'direct_z'
        result=subprocess.run(['gfortran','-O3','-fdefault-real-8','-I',str(ROOT/'include'),
                              str(ROOT/'tests/plasma_z_direct_probe.f90'),
                              *[str(ROOT/f'{m}.o') for m in MODULES],*c_objects(ROOT),'-o',str(cls.direct)],capture_output=True,text=True)
        if result.returncode:raise RuntimeError(result.stderr)

    @classmethod
    def tearDownClass(cls):
        for p,digest in cls.hashes.items():
            if hashlib.sha256(p.read_bytes()).digest()!=digest:raise AssertionError(f'Production source changed: {p}')

    def test_original_linked_implementation_and_checked_build(self):
        points=[(k,z) for k in [-1.,1.] for z in [0j,1+.1j,1-1j,3.01+.1j,10-.2j,-2.9+1e-12j]]
        data=''.join(f'{k} {z.real} {z.imag}\n' for k,z in points)
        result=subprocess.run([self.direct],input=data,text=True,capture_output=True,check=True)
        values=[complex(*map(float,line.split())) for line in result.stdout.splitlines()]
        for (k,z),value in zip(points,values):
            for h in self.harnesses.values():
                np.testing.assert_allclose(h.z(z,k=k)['value'],value,rtol=2e-14,atol=2e-14)

    def test_series_tolerance_and_asymptotic_floor(self):
        h=self.harnesses[False]
        error=[]
        for tolerance in TOLERANCES:
            value=h.z(1+.1j,method=1,tolerance=tolerance)
            error.append(abs(value['value']-scipy_z(1.,1+.1j)))
            self.assertEqual(value['branch'],1)
        self.assertLess(error[-1],2e-14)
        self.assertGreater(error[0],error[-1]*100)
        for h in self.harnesses.values():
            for tolerance in TOLERANCES:
                value=h.z(3.01+.1j,method=1,tolerance=tolerance)
                self.assertEqual((value['branch'],value['stop']),(2,2))
                self.assertGreater(abs(value['value']-scipy_z(1.,3.01+.1j)),1e-4)

    def test_landau_symmetries_and_negative_wavevector(self):
        h=self.harnesses[False]
        for z in [.2+.1j,1.2+.8j,3.5+.2j,2.+1.2j]:
            a=h.z(z)['value'];reflected=h.z(-z)['value']
            conjugate=h.z(-z.conjugate())['value']
            reference=scipy_z(1.,z)
            self.assertLess(abs(a-reference),1e-5)
            np.testing.assert_allclose(h.z(z,k=-1.)['value'],-reflected,rtol=1e-14,atol=1e-14)
            np.testing.assert_allclose(conjugate,-a.conjugate(),rtol=1e-14,atol=1e-14)
            np.testing.assert_allclose(a+reflected,2j*np.sqrt(np.pi)*np.exp(-z*z),rtol=2e-5,atol=2e-7)
            # The independent reference also obeys the entire-function identity.
            np.testing.assert_allclose(scipy_z(1.,z)+scipy_z(1.,-z),
                                      2j*np.sqrt(np.pi)*np.exp(-z*z),rtol=2e-13,atol=2e-13)

    def test_failure_paths_and_hard_iteration_cap(self):
        for h in self.harnesses.values():
            value=h.z(-14j,method=1)
            self.assertEqual(value['stop'],1)
            self.assertTrue(value['nonfinite'])
            self.assertTrue(np.isfinite(scipy_z(1.,-14j)))
            zero=h.z(1+.1j,k=0,method=1)
            self.assertEqual(zero['stop'],4)
            self.assertTrue(zero['nonfinite'])
            # Tightened continued-fraction stopping can reach the historical cap.
            results=[h.z(z,method=1,tolerance=1e-15) for z in [1j,1+1j,2+1j]]
            self.assertTrue(any(r['stop']==3 for r in results))
            for r in results:self.assertLessEqual(r['iterations'],30)

    def test_reproduce_stress_and_substitute_only_z(self):
        for h in self.harnesses.values():
            kwargs=dict(kpar=-.31,kperp=3.,beta=.4,alpha_scale=1.7)
            old=h.tensor(S,W,**kwargs);new=h.tensor(S,W,method=2,**kwargs)
            expected=species_tensor(S,W,-.31,3.,.02,1024)
            self.assertAlmostEqual(metrics(old['chi'],expected)['scaled'],7.624471347e-4,delta=2e-10)
            self.assertLess(metrics(new['chi'],expected)['scaled'],2e-11)
            np.testing.assert_array_equal(old['arguments'],new['arguments'])
            for tol in TOLERANCES:
                value=h.tensor(S,W,method=1,tolerance=tol,**kwargs)
                self.assertGreater(metrics(value['chi'],expected)['scaled'],7.6e-4)
            smaller=h.tensor(S,W,method=2,**{**kwargs,'n':128})
            np.testing.assert_allclose(smaller['chi'],new['chi'],rtol=1e-13,atol=1e-12)

    def test_affine_error_budget_and_cancellation(self):
        h=self.harnesses[False]
        old=h.tensor(S,W,kpar=-.31,kperp=3.);new=h.tensor(S,W,kpar=-.31,kperp=3.,method=2)
        c,b,extra=stix_terms(S,W,-.31,3.,.02,1024)
        prediction=np.sum(b*(old['z']-new['z']),axis=-1)
        np.testing.assert_allclose(prediction,old['chi']-new['chi'],rtol=1e-8,atol=2e-11)
        component=(2,2)
        dominant=int(np.argmax(abs(b[component]*(old['z']-new['z']))))-1024
        self.assertEqual(dominant,3)
        raw=(np.sum(abs(c[component])+abs(b[component]*new['z']))+abs(extra[component]))/abs(new['chi'][component])
        self.assertGreater(raw,10.)

    def test_auxiliary_scale_and_collision_controls(self):
        h=self.harnesses[False]
        expected=species_tensor(S,W,-.31,3.,.02,1024)
        outputs=[]
        for beta in [.01,.1,1.,10.]:
            value=h.tensor(S,W,kpar=-.31,kperp=3.,beta=beta,method=2)
            outputs.append(value['chi'])
            self.assertLess(metrics(value['chi'],expected)['scaled'],2e-11)
        for value in outputs[1:]:np.testing.assert_allclose(value,outputs[0],rtol=2e-12,atol=2e-11)
        for nu in [1e-4,.01,.1]:
            old=h.tensor(S,W,kpar=-.31,kperp=3.,collision=nu)
            new=h.tensor(S,W,kpar=-.31,kperp=3.,collision=nu,method=2)
            self.assertFalse(old['nonfinite']);self.assertFalse(new['nonfinite'])
            self.assertGreater(metrics(old['chi'],new['chi'])['scaled'],1e-8)

    def test_nearly_real_branch_continuity_defect(self):
        for h in self.harnesses.values():
            upper=h.z(3.01+1e-12j)['value']
            real=h.z(3.01+0j)['value']
            lower=h.z(3.01-1e-12j)['value']
            self.assertGreater(abs(upper-real),1e-4)
            self.assertGreater(abs(lower-real),1e-4)
            reference=[scipy_z(1.,3.01+1e-12j),scipy_z(1.,3.01+0j),scipy_z(1.,3.01-1e-12j)]
            self.assertLess(max(abs(v-reference[1]) for v in reference),1e-11)

    def test_collisional_zero_rate_anisotropy_control(self):
        h=self.harnesses[False]
        for alpha in [1.,1.3]:
            species=tuple([*S[:2],alpha,*S[3:]])
            base=h.tensor(species,W,kpar=-.31,kperp=3.,method=2)
            col=h.tensor(species,W,kpar=-.31,kperp=3.,method=2,collisional=True)
            # The historical audit records the pre-correction defect. Current
            # collisional zz must recover the collisionless limit for both
            # isotropic and anisotropic distributions.
            np.testing.assert_allclose(base['chi'],col['chi'],rtol=2e-13,atol=2e-11)

    def test_grid_is_deterministic_and_spans_regimes(self):
        a=z_grid();b=z_grid()
        self.assertEqual(a,b)
        self.assertGreater(len(a),150)
        for name in ['small_series','nearly_real','continued_fraction','asymptotic',
                     'radius_3_boundary','imaginary_1_boundary','large_arguments','near_zero','domain_and_range']:
            self.assertTrue(any(label==name for label,z in a),name)
