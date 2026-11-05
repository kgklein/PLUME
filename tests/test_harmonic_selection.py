"""Checkpoint A: independent positive-tail verification; no production kernel."""
from pathlib import Path
import subprocess
import tempfile
import unittest
import numpy as np
from scipy.special import ive
ROOT = Path(__file__).resolve().parents[1]

class HarmonicSelectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.tmp.cleanup)
        cls.exe = Path(cls.tmp.name)/'select'
        subprocess.run(['gfortran','-O2','-fdefault-real-8','-fcheck=all',
                        '-ffpe-trap=invalid,zero,overflow', '-J',cls.tmp.name,
                        str(ROOT/'src/harmonic_selection.f90'),
                        str(ROOT/'tests/harmonic_selection_probe.f90'),'-o',str(cls.exe)],check=True)

    def select(self,x,t=1e-10,cap=512):
        out=subprocess.check_output([self.exe],input=f'{x} {t} {cap}\n',text=True).split()
        return int(out[0]),int(out[1]),int(out[2]),float(out[3]),float(out[4])

    def test_positive_tails_and_adjacent_order(self):
        for x in [1e-12,1e-6,.01,.5,1,10,100,1000,10000]:
            for tol in [1e-4,1e-10,1e-14,1e-17]:
                with self.subTest(x=x,tol=tol):
                    status,n,ok,b0,b2=self.select(x,tol,10000)
                    self.assertEqual((status,ok),(0,1))
                    orders=np.arange(n+1, max(n+200,int(15*np.sqrt(x))+200))
                    coefficients=ive(orders,x)
                    tails=(2*coefficients.sum(),2*np.sum(orders**2*coefficients)/x)
                    self.assertLess(max(tails),tol)
                    self.assertLessEqual(tails[0],b0*(1+1e-12))
                    self.assertLessEqual(tails[1],b2*(1+1e-12))
                    self.assertTrue(np.isfinite(ive(n+1,x)))

    def test_zero_tiny_and_cap(self):
        self.assertEqual(self.select(0.)[:3],(0,0,1))
        self.assertEqual(self.select(1e-300)[:3],(0,1,1))
        self.assertEqual(self.select(1000.,1e-10,2)[:3],(0,2,0))

    def test_monotonic_and_invalid(self):
        counts=[self.select(x)[1] for x in [1e-8,.01,1,10,100,1000]]
        self.assertEqual(counts,sorted(counts))
        counts=[self.select(30.,t)[1] for t in [1e-3,1e-6,1e-10,1e-14]]
        self.assertEqual(counts,sorted(counts))
        for args in [(-1.,1e-10,10),(1.,0.,10),(1.,1.,10),(1.,1e-10,-1)]:
            self.assertEqual(self.select(*args)[0],1)
