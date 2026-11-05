"""Milestone 1: vendor source, C ABI, standalone accuracy and native intervention."""
import json
import struct
import subprocess
import unittest
import numpy as np
from scipy.special import wofz
from faddeeva_harness import FaddeevaHarness,audit_arguments,ROOT
from plasma_dispersion_audit import scipy_z,metrics
from mixed_reference import species_tensor


class FaddeevaStandaloneTests(unittest.TestCase):
    maxima={'grid':0.,'mpmath':0.,'stress':0.}
    @classmethod
    def setUpClass(cls):
        cls.h={False:FaddeevaHarness(),True:FaddeevaHarness(checked=True)}
        for h in cls.h.values():cls.addClassCleanup(h.close)
        cls.points=audit_arguments(cls.h[False])

    @classmethod
    def tearDownClass(cls):print('Faddeeva standalone scaled errors:',cls.maxima)

    def test_entire_audit_grid(self):
        h=self.h[False]
        for regime,z in self.points:
            out=h.accurate_z(z)
            with np.errstate(all='ignore'):ref=scipy_z(1.,z)
            if not np.isfinite(ref):
                self.assertNotEqual(out['status'],0,(regime,z,out))
                continue
            self.assertEqual(out['status'],0,(regime,z,out))
            error=abs(out['value']-ref)/max(1.,abs(ref))
            type(self).maxima['grid']=max(type(self).maxima['grid'],error)
            self.assertLess(error,1e-13,(regime,z,out,ref,error))
            self.assertEqual(out['c_complex_bits'],128)
            self.assertEqual(out['plume_complex_bits'],128)

    def test_signed_convention_checked_and_frozen_high_precision(self):
        data=json.loads((ROOT/'tests/PLASMA_DISPERSION_AUDIT_RESULTS.json').read_text())
        for h in self.h.values():
            for row in data['mpmath_z']:
                z=complex(*row['z']);ref=complex(*row['value'])
                out=h.accurate_z(z)
                self.assertEqual(out['status'],0,(z,out))
                error=abs(out['value']-ref)/max(1.,abs(ref))
                type(self).maxima['mpmath']=max(type(self).maxima['mpmath'],error)
                self.assertLess(error,1e-13)
            for z in [0j,3.01+1e-12j,3.01-1e-12j,-14j,200-200j,1.+1.j]:
                for k in [-.31,.31]:
                    out=h.accurate_z(z,k=k)
                    self.assertEqual(out['status'],0)
                    np.testing.assert_allclose(out['value'],scipy_z(k,z),rtol=1e-13,atol=1e-13)

    def test_invalid_range_and_recovery(self):
        for h in self.h.values():
            for z,k,expected in [(1+.1j,0.,2),(complex(np.nan,0),1.,1),
                                 (complex(np.inf,0),1.,1),(1+.1j,np.inf,1),(-30j,1.,3)]:
                bad=h.accurate_z(z,k=k)
                self.assertEqual(bad['status'],expected,(z,k,bad))
                self.assertTrue(np.isnan(bad['value'].real) and np.isnan(bad['value'].imag))
                self.assertEqual(h.accurate_z(1+.1j)['status'],0)
            self.assertEqual(h.accurate_z(-26j)['status'],0)

    def test_native_stress_intervention(self):
        s=(2.,.25,1.3,.5,.05,.4);w=.73+.13j
        for h in self.h.values():
            result=h.accurate_tensor(s,w,kpar=-.31,kperp=3.,beta=.4,alpha_scale=1.7)
            ref=species_tensor(s,w,-.31,3.,.02,1024)
            error=metrics(result['chi'],ref)['scaled']
            type(self).maxima['stress']=max(type(self).maxima['stress'],error)
            self.assertLess(error,2e-13)

    def test_fast_link_subnormal_result(self):
        h=self.h[False];exe=h.directory/'direct_fast'
        subprocess.run(['gfortran','-O3','-ffast-math','-fdefault-real-8','-I',str(h.directory),
                        str(ROOT/'tests/faddeeva_direct_probe.f90'),str(h.directory/'plasma_z_faddeeva.o'),
                        str(h.directory/'faddeeva.o'),str(h.directory/'faddeeva_bridge.o'),'-lm','-o',str(exe)],check=True,capture_output=True)
        result=subprocess.run([exe],input='1 1e308 0\n1 0 -14\n',text=True,capture_output=True,check=True)
        rows=[line.split() for line in result.stdout.splitlines()]
        self.assertEqual(rows[0][0],'0')
        re=struct.unpack('d',struct.pack('Q',int(rows[0][-2],16)))[0]
        self.assertNotEqual(re,0.)
        np.testing.assert_allclose(re,-1e-308,rtol=1e-14,atol=0.)
