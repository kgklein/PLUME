"""Checkpoint B: isolated storage, invalidation, precision and cap behavior."""
from pathlib import Path
import subprocess
import tempfile
import unittest
import numpy as np
from scipy.special import ive
ROOT=Path(__file__).resolve().parents[1]

class BesselCacheTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp=tempfile.TemporaryDirectory(); cls.addClassCleanup(cls.tmp.cleanup)
        cls.exe=Path(cls.tmp.name)/'cache'
        subprocess.run(['gfortran','-O2','-fdefault-real-8','-fcheck=all',
                        '-ffpe-trap=invalid,zero,overflow','-J',cls.tmp.name,
                        str(ROOT/'src/harmonic_selection.f90'),str(ROOT/'src/bessel_cache.f90'),
                        str(ROOT/'tests/bessel_cache_probe.f90'),'-o',str(cls.exe)],check=True)

    def run_cache(self,records):
        output=subprocess.check_output([self.exe],input='\n'.join(records)+'\n',text=True)
        snapshots=[]
        for block in output.split('STATUS ')[1:]:
            lines=block.splitlines(); self.assertEqual(lines[0],'0')
            rows=[line.split() for line in lines[2:] if line and line[0]!='#' and not line.startswith(('VALUES','WARNING'))]
            values=next(line[7:] for line in lines if line.startswith('VALUES '))
            snapshots.append((rows,np.fromstring(values,sep=' ')))
        return output,snapshots

    def test_coefficients_and_adjacent_orders(self):
        xs=[0.,1e-300,1e-12,.01,1.,10.,100.,1000.,10000.,1000000.]
        _,results=self.run_cache([f'1 1 {x} adaptive 1e-12 2000 128' for x in xs])
        for x,(_,v) in zip(xs,results):
            np.testing.assert_allclose(v,ive(np.arange(len(v)),x),atol=2e-16,rtol=2e-12)

    def test_reuse_resize_and_tolerance(self):
        _,rows=self.run_cache(['1 2 1 adaptive 1e-10 512 128',
                              '1 2 1 adaptive 1e-10 512 128',
                              '1 2 100 adaptive 1e-10 512 128',
                              '1 2 .01 adaptive 1e-10 512 128',
                              '1 2 .01 adaptive 1e-5 512 128',
                              '1 2 .01 adaptive 1e-14 512 128',
                              '1 3 .01 adaptive 1e-14 512 128',
                              '2 3 100 adaptive 1e-10 512 128',
                              '1 1 100 adaptive 1e-10 512 128'])
        stats=[r[0][0] for r in rows]
        self.assertEqual(stats[1][-3:],['1','1','1'])
        self.assertGreater(int(stats[2][3]),int(stats[1][3]))
        self.assertEqual(stats[2][3],stats[5][3])
        self.assertEqual(stats[3][-1],stats[4][-1])
        self.assertEqual(int(stats[5][-1]),int(stats[4][-1])+1)
        self.assertEqual(stats[6][-1],stats[5][-1])
        np.testing.assert_array_equal(rows[0][1],rows[1][1])

    def test_cap_warning_once_and_recovery(self):
        out,rows=self.run_cache(['1 1 1000 adaptive 1e-10 3 128']*3+
                               ['1 1 1000 adaptive 1e-10 512 128'])
        self.assertEqual(out.count('WARNING:'),1)
        self.assertEqual(rows[0][0][0][2],'3')
        self.assertEqual(rows[0][0][0][4],'F')
        self.assertEqual(rows[-1][0][0][4],'T')

    def test_fixed_range(self):
        _,rows=self.run_cache(['1 1 10 fixed 1e-10 512 256'])
        self.assertEqual(rows[0][0][0][2],'256')
        self.assertEqual(len(rows[0][1]),258)
