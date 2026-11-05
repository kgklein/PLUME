from faddeeva_build import compile_c,c_objects
"""Stage 3A: native cold-reference assembly and existing production root finder.

Oracles use the Stage 1 momentum/continuity solve and literal circular response,
never the production kernel. No legacy numerical expectations are changed.
"""
from pathlib import Path
import subprocess
import tempfile
import sys
import unittest

import numpy as np
from scipy.optimize import brentq, root as solve_root
from test_cold_susceptibility import Case, paper_tensor
from mixed_reference import species_tensor, wave_matrix
from test_mixed_dispersion import ROOT, input_text, run, parse, replace_species

P = (0., 1., 1., 1., 1., 0.)
E = (0., 1836.15267343, 1., -1., 1., 0.)
ALPHA = (0., .25, 0., .5, .05, 2.)


def circular(species, w, va):
    """Stationary R, L, P including identity, literal finite-mass formulas."""
    r = l = p = 1.+0j
    for _, mu, _, q, density, _ in species:
        a = density*mu/(q*q*va*va)
        g = mu/q
        r -= a/(w*(w+g))
        l -= a/(w*(w-g))
        p -= a/(w*w)
    return r, l, p


class FullyColdDispersionTests(unittest.TestCase):
    tensor_error = 0.
    root_error = 0.
    benchmark_roots = {}

    @classmethod
    def setUpClass(cls):
        subprocess.run(['make','-j2'], cwd=ROOT, check=True, capture_output=True)
        cls.temp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.temp.cleanup)
        cls.probes = {}
        modules = ('nrtype','nrutil_trim','vars','plasma_z_faddeeva','plasma_z_legacy','plasma_z_backend','cold_susceptibility','auxiliary_normalization','kinetic_core','harmonic_selection','bessel_cache','adaptive_kinetic','cold_reference','electromagnetic_fields','scan_support','functions','bessel','disprels')
        for checked in (False, True):
            directory = Path(cls.temp.name)/str(checked)
            directory.mkdir()
            exe = directory/'probe'
            if checked:
                compile_c(directory)
                command = ['gfortran','-O0','-g','-fdefault-real-8','-fcheck=all',
                           '-ffpe-trap=invalid,zero,overflow','-J',str(directory),'-I',str(directory),
                           *[str(ROOT/'src'/f'{name}.f90') for name in modules],
                           str(ROOT/'tests/fully_cold_probe.f90')]
            else:
                command = ['gfortran','-O3','-fdefault-real-8','-I',str(ROOT/'include'),
                           str(ROOT/'tests/fully_cold_probe.f90'),
                           *[str(ROOT/f'{name}.o') for name in modules]]
            subprocess.run(command+c_objects(directory if checked else ROOT)+['-o',str(exe)], cwd=directory, check=True, capture_output=True)
            cls.probes[checked] = exe

    @classmethod
    def tearDownClass(cls):
        print(f'Fully cold: maximum element scaled error={cls.tensor_error:.3e}; '
              f'maximum root difference={cls.root_error:.3e}; roots={cls.benchmark_roots}')

    def evaluate(self, species, frequencies, checked=False, roots=False, **kwargs):
        stdin = ''.join(f'{w.real:.17g} {w.imag:.17g} {int(roots)}\n' for w in frequencies)
        result, _ = run(self.probes[checked], input_text(species,**kwargs), stdin)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        blocks, found = parse(result.stdout, len(species))
        self.assertEqual(len(found) if roots else len(blocks), len(frequencies), result.stdout)
        return blocks, found, result.stdout

    def compare_tensor(self, actual, expected):
        error = np.max(abs(actual-expected)/np.maximum(1.,abs(expected)))
        type(self).tensor_error = max(type(self).tensor_error,float(error))
        np.testing.assert_allclose(actual, expected, rtol=5e-13, atol=5e-13)

    def test_tensor_assembly(self):
        systems = [[P,E], [P,replace_species(E,density=1.1,drift=.2/1.1),ALPHA],
                   [replace_species(P,drift=-.2),replace_species(E,density=.85,drift=-.1),
                    replace_species(ALPHA,drift=-.8), (0.,.13,-8.,-.4,.07,.3)]]
        for checked in self.probes:
            for species in systems:
                for va in (.001,.03,.1):
                    for kpar,kperp in ((.2,0.),(0.,.3),(.17,-.21),(-.17,.21),(0.,0.)):
                        frequencies = [.73+.19j,.73-.19j]
                        blocks,_,_ = self.evaluate(species, frequencies, checked=checked,
                                                   va=va,kpar=kpar,kperp=kperp)
                        for w,b in zip(frequencies,blocks):
                            self.assertEqual(b['status'],0)
                            chi = np.array([species_tensor(s,w,kpar,kperp,va) for s in species])
                            self.compare_tensor(b['chi'],chi)
                            # A second oracle uses literal paper equations, not a fluid solve.
                            paper = np.array([paper_tensor(Case(omega=w,kpar=kpar,kperp=kperp,va=va,
                                                mu=s[1],charge=s[3],density=s[4],drift=s[5])) for s in species])
                            self.compare_tensor(b['chi'],paper)
                            epsilon = np.eye(3)+chi.sum(axis=0)
                            n = np.array([kperp,0.,kpar])/(va*w)
                            self.compare_tensor(np.eye(3)+b['chi'].sum(axis=0),epsilon)
                            self.compare_tensor(b['lam'],epsilon+np.outer(n,n)-np.eye(3)*(n@n))
                            np.testing.assert_allclose(b['det'],np.linalg.det(b['lam']),rtol=5e-12)

    def check_root(self, species, expected, guess, *, kpar,kperp,va=.01,label='',checked=False):
        blocks,roots,_ = self.evaluate(species,[guess],roots=True,checked=checked,
                                       kpar=kpar,kperp=kperp,va=va)
        flag,w = roots[0]
        self.assertGreater(flag,0)
        self.assertLess(flag,75)
        self.assertEqual(blocks[0]['status'],0)
        np.testing.assert_allclose([w.real,w.imag],[expected.real,expected.imag],rtol=5e-10,atol=2e-10)
        type(self).root_error = max(type(self).root_error,abs(w-expected))
        if label:
            type(self).benchmark_roots[label] = w
        self.assertLess(abs(blocks[0]['det'])/np.max(abs(blocks[0]['lam']))**3,1e-10)
        # Independent residual, scaled by the tensor norm, not a frequency floor.
        matrix = wave_matrix(species,w,kpar,kperp,va)
        self.assertLess(abs(np.linalg.det(matrix/np.max(abs(matrix)))),1e-10)
        return w

    def independent_root(self,species,guess,kpar,kperp,va):
        def residual(x):
            w = complex(*x)
            d = np.linalg.det(wave_matrix(species,w,kpar,kperp,va)*(va*w)**2)
            return [d.real,d.imag]
        solution = solve_root(residual,[guess.real,guess.imag],tol=1e-11)
        self.assertTrue(solution.success,solution.message)
        w = complex(*solution.x)
        matrix = wave_matrix(species,w,kpar,kperp,va)
        self.assertLess(abs(np.linalg.det(matrix/np.max(abs(matrix)))),1e-12)
        return w

    def test_parallel_circular_roots(self):
        species=[P,E]
        for k in (.2,.6):
            for index,label in ((0,'R'),(1,'L')):
                va=.01
                f=lambda w: (circular(species,w,va)[index]-(k/(va*w))**2).real
                expected=complex(brentq(f,.01,.99,xtol=1e-14))
                for sign in (-1,1):
                    w=self.check_root(species,expected,.98*expected+.001j,kpar=sign*k,kperp=0.,
                                      label=f'{label} k={sign*k}',checked=sign<0)
                    r,l,p=circular(species,w,va)
                    self.assertLess(abs((r,l)[index]-(k/(va*w))**2)/abs((k/(va*w))**2),1e-9)

    def test_perpendicular_ordinary_extraordinary_roots(self):
        species=[P,E]; va=.1; k=.2
        ordinary=complex(np.sqrt(sum(s[4]*s[1]/s[3]**2 for s in species)/va**2+(k/va)**2))
        def extraordinary(w):
            r,l,p=circular(species,w,va)
            return (r*l/((r+l)/2)-(k/(va*w))**2).real
        x=complex(brentq(extraordinary,.01,.4,xtol=1e-14))
        for checked in self.probes:
            for expected,label in ((ordinary,'O'),(x,'X')):
                self.check_root(species,expected,.999*expected+.0001j,kpar=0.,kperp=k,va=va,
                                label=label,checked=checked)
        self.assertLess(abs(extraordinary(x)),1e-8)

    def test_oblique_roots(self):
        for angle in (15,45,75):
            for magnitude in (.12,.4):
                for sign in (-1,1):
                    kpar=sign*magnitude*np.cos(np.deg2rad(angle))
                    kperp=magnitude*np.sin(np.deg2rad(angle))
                    guess=abs(kpar)*.9+.001j
                    expected=self.independent_root([P,E],guess,kpar,kperp,.01)
                    self.check_root([P,E],expected,guess,kpar=kpar,kperp=kperp,
                                    label=f'{angle}deg k={sign*magnitude}')

    def test_unstable_drifting_alpha_roots(self):
        # Charge neutral and current balanced, but a parallel alpha beam is unstable.
        species=[P,replace_species(E,density=1.1,drift=.2/1.1),ALPHA]
        for sign in (-1,1):
            directed=[replace_species(s,drift=sign*s[5]) for s in species]
            for growth in (-1,1):
                guess=.3+growth*.08j
                expected=self.independent_root(directed,guess,sign*.3,.1,.01)
                self.assertGreater(growth*expected.imag,.07)
                w=self.check_root(directed,expected,guess,kpar=sign*.3,kperp=.1,
                                  label=f'beam drift={sign}, growth={growth}',checked=growth>0)
                self.assertGreater(growth*w.imag,.07)

    def test_no_thermal_dependency_checked(self):
        # Probe poisons all solver_* scalars with NaN and asserts solver_spec stays
        # unallocated. Checked executable traps any accidental thermal arithmetic.
        for kpar,kperp in ((0.,.2),(.2,0.),(.15,.2)):
            baseline=None
            for alpha in (0.,1.,-10.,1e100):
                species=[replace_species(P,alpha=alpha),replace_species(E,alpha=-alpha)]
                blocks,_,_=self.evaluate(species,[.3+.04j],checked=True,kpar=kpar,kperp=kperp)
                self.assertEqual(blocks[0]['status'],0)
                if baseline is not None:
                    np.testing.assert_array_equal(blocks[0]['chi'],baseline['chi'])
                    np.testing.assert_array_equal(blocks[0]['lam'],baseline['lam'])
                    self.assertEqual(blocks[0]['det'],baseline['det'])
                baseline=blocks[0]

    def test_poles_recovery_and_failed_root(self):
        species=[P,E,replace_species(ALPHA,drift=.5)]
        ws=[.7+.01j,0j,.125+0j,.625+0j,-.375+0j,1+0j,.7+.01j]
        for checked in self.probes:
            blocks,_,_=self.evaluate(species,ws,checked=checked,kpar=.25)
            self.assertEqual([b['status'] for b in blocks],[0,2,3,4,4,4,0])
            for b in blocks[1:-1]:
                self.assertTrue(np.isnan(b['chi']).all())
                self.assertTrue(np.isnan(b['lam']).all())
                self.assertTrue(np.isnan(b['det']))
            np.testing.assert_array_equal(blocks[0]['lam'],blocks[-1]['lam'])
            _,roots,_=self.evaluate([P,E],[0j,.1+.001j],checked=checked,roots=True,kpar=.1,kperp=0.)
            self.assertEqual(roots[0][0],-1)
            self.assertGreater(roots[1][0],0)

    def test_warning_only_equilibria(self):
        species=[replace_species(P,drift=.1),replace_species(E,density=.8,drift=-.1)]
        _,_,stdout=self.evaluate(species,[.5+.05j,.6-.01j])
        self.assertEqual(stdout.count('WARNING: Plasma not charge neutral'),1)
        self.assertEqual(stdout.count('WARNING: Nonzero equilibrium current'),1)

    def test_public_boundaries(self):
        source=input_text([P,E])
        bad=[(source.replace("'d_ref'","'rho_ref'"),'require normalization=d_ref'),
             (source.replace('betap=0','betap=0.2'),'betap=0'),
             (input_text([P,replace_species(E,tau=.1)],kpar=0.),'nonzero kpar'),
             (input_text([P,replace_species(E,tau=-.1)]),'negative tauS'),
             (input_text([P,E],extra='collision_type=1'),'collision_type=0'),
             (input_text([P,E],option=6),'FPC'),(input_text([P,E],option=7),'FPC'),
             (input_text([P,E],option=5),'scans'),(source.replace('nscan=0','nscan=1'),'scans'),
             (input_text([P,E],option=2),'option=-1, 0, or 1'),
             (input_text([P,E],va=0),'positive vtp')]
        for text,message in bad:
            result,_=run(ROOT/'plume.e',text)
            self.assertNotEqual(result.returncode,0)
            self.assertIn(message,result.stdout+result.stderr)
        result,files=run(ROOT/'plume.e',source)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        self.assertIn('raw det(Lambda)',result.stdout)
        self.assertIn('determinant only',result.stdout)
        self.assertIn('no thermal reference scale',files['normalization_test.txt'])

    def test_map_with_poles(self):
        # Exact perpendicular grid includes w=0 and both ion cyclotron poles.
        source=input_text([P,E],option=0,kpar=0.,kperp=.2)
        source+=('&maps\nloggridw=.false.\nloggridg=.false.\nomi=-1.25\nomf=1.25\n'
                 'gami=-.1\ngamf=.1\nnr=10\nni=2\npositive_roots=.false.\n/\n')
        result,files=run(ROOT/'plume.e',source)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        self.assertIn('singular/nonfinite mixed map points',result.stdout)
        rows=np.array([list(map(float,l.split())) for l in files['dispersion_test.map'].splitlines() if l.strip()])
        self.assertTrue(np.isnan(rows[:,2:]).all(axis=1).any())
        self.assertTrue(np.isfinite(rows).all(axis=1).any())
        for line in files['dispersion_test.roots'].splitlines():
            self.assertTrue(np.isfinite(list(map(float,line.split()))).all())

    def test_map_recovers_perpendicular_root(self):
        source=input_text([P,E],option=0,kpar=0.,kperp=.2)
        source+=('&maps\nloggridw=.false.\nloggridg=.false.\nomi=.12\nomf=.3\n'
                 'gami=-.01\ngamf=.01\nnr=36\nni=10\npositive_roots=.true.\n/\n')
        result,files=run(ROOT/'plume.e',source)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        rows=[list(map(float,l.split())) for l in files['dispersion_test.roots'].splitlines() if l.strip()]
        self.assertGreater(len(rows),0)
        expected=self.independent_root([P,E],.199+.001j,0.,.2,.001)
        # Root files print only five significant digits.
        self.assertTrue(any(abs(complex(*row[4:6])-expected)<5e-6 for row in rows))

    def test_example_and_wrapper(self):
        sys.path.insert(0,str(ROOT))
        from linfpclib.linfpc import plume_input
        for name in ('cold_proton_electron.in','cold_proton_electron_perpendicular.in'):
            source=(ROOT/'inputs/example'/name).read_text().replace("dataName='example'","dataName='mixed'")
            result,_=run(ROOT/'plume.e',source)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            self.assertIn('Dispersion Solutions',result.stdout)
        wrapped=plume_input('mixed')
        wrapped.set_params(0.,.2,0.,.01,2,0,1,1,False,True,normalization='d_ref')
        wrapped.make_species(*P); wrapped.make_species(*E)
        wrapped.make_guess(.19,.001)
        with tempfile.TemporaryDirectory() as directory:
            name=str(Path(directory)/'cold')
            wrapped.write_input(name,'test')
            result,_=run(ROOT/'plume.e',Path(name+'.in').read_text())
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        wrapped.species[1]['tauS']=.1
        # A warm perpendicular species remains unsupported after Stage 3B.2.
        with self.assertRaisesRegex(ValueError,'nonzero kpar'):
            wrapped.check_input_for_errors()


if __name__ == '__main__':
    unittest.main()
