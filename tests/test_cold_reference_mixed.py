from faddeeva_build import compile_c,c_objects
"""Milestone B: native cold-reference dispatch, independent tensors and roots."""
from pathlib import Path
import subprocess
import tempfile
import unittest
import sys
import numpy as np
from scipy.optimize import root as solve_root
from mixed_reference import species_tensor,wave_matrix
import test_mixed_dispersion as mixed

P=(0.,1.,0.,1.,1.,0.)
E=(.2,1836.,.8,-1.,1.,0.)
ALPHA=(.03,.25,1.2,.5,.05,.3)


class ColdReferenceMixedTests(unittest.TestCase):
    tensor_error=0.
    root_error=0.
    order_error=0.
    maximum_residual=0.
    roots={}

    @classmethod
    def setUpClass(cls):
        mixed.MixedDispersionTests.setUpClass()
        cls.addClassCleanup(mixed.MixedDispersionTests.temp.cleanup)
        cls.probes={False:mixed.MixedDispersionTests.probe,True:mixed.MixedDispersionTests.checked_probe}
        # Poison thermal globals and assert that solver_spec is never allocated.
        cls.native_probe=cls.probes[True].parent/'native_probe'
        modules=('nrtype','nrutil_trim','vars','plasma_z_faddeeva','plasma_z_legacy','plasma_z_backend','cold_susceptibility','auxiliary_normalization',
                 'kinetic_core','harmonic_selection','bessel_cache','adaptive_kinetic','cold_reference','functions','bessel','disprels')
        directory=cls.probes[True].parent
        subprocess.run(['gfortran','-O0','-g','-fdefault-real-8','-fcheck=all','-ffpe-trap=invalid,zero,overflow',
                        '-I',str(directory),'-J',str(directory),
                        *[str(mixed.ROOT/'src'/f'{m}.f90') for m in modules],
                        str(mixed.ROOT/'tests/fully_cold_probe.f90'),*c_objects(directory),'-o',str(cls.native_probe)],
                       cwd=directory,check=True,capture_output=True)

    @classmethod
    def tearDownClass(cls):
        print(f'Cold-reference mixed: tensor error={cls.tensor_error:.3e}; root error={cls.root_error:.3e}; '
              f'order error={cls.order_error:.3e}; residual={cls.maximum_residual:.3e}; roots={cls.roots}')

    def evaluate(self,species,frequencies,*,checked=False,root=False,native=False,**kwargs):
        result,files=mixed.run(self.native_probe if native else self.probes[checked],
                               mixed.input_text(species,**kwargs),
                               ''.join(f'{w.real:.17g} {w.imag:.17g} {int(root)}\n' for w in frequencies))
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        blocks,roots=mixed.parse(result.stdout,len(species))
        self.assertEqual(len(roots) if root else len(blocks),len(frequencies),result.stdout)
        return blocks,roots,result.stdout,files

    def compare_tensor(self,a,b,order=False):
        error=float(np.max(abs(a-b)/np.maximum(1.,abs(b))))
        if order:
            type(self).order_error=max(self.order_error,error)
        else:
            type(self).tensor_error=max(self.tensor_error,error)
        np.testing.assert_allclose(a,b,rtol=5e-9 if order else 2e-7,atol=1e-5)

    def test_species_and_wave_assembly(self):
        cases=[[P,E], [P,mixed.replace_species(E,drift=.15)],
               [P,mixed.replace_species(E,alpha=.65,drift=-.15)],
               [P,mixed.replace_species(E,density=1.1),ALPHA],
               [P,mixed.replace_species(E,tau=0,alpha=0),ALPHA],
               [P,E,mixed.replace_species(ALPHA,tau=0,drift=-.3)],
               [P,mixed.replace_species(E,tau=2.),mixed.replace_species(ALPHA,tau=.002),(.2,.1,.7,-.5,.03,-.1)],
               [mixed.replace_species(P,drift=.1),mixed.replace_species(E,density=.8,drift=-.1)]]
        for checked in self.probes:
            for species in cases:
                for kpar in (.06,-.06):
                    for kperp in (0.,.05):
                        for va in (.001,.05):
                            ws=[.13+.03j,.13-.03j]
                            blocks,_,_,_=self.evaluate(species,ws,checked=checked,kpar=kpar,kperp=kperp,va=va)
                            for w,b in zip(ws,blocks):
                                self.assertEqual(b['status'],0)
                                chi=np.array([species_tensor(s,w,kpar,kperp,va) for s in species])
                                self.compare_tensor(b['chi'],chi)
                                self.compare_tensor(np.eye(3)+b['chi'].sum(axis=0),np.eye(3)+chi.sum(axis=0))
                                self.compare_tensor(b['lam'],wave_matrix(species,w,kpar,kperp,va))
                                np.testing.assert_allclose(b['det'],np.linalg.det(b['lam']),rtol=5e-12)

    def check_root(self,species,guess,*,kpar=.06,kperp=.05,va=.001,label='',checked=False):
        def residual(x):
            w=complex(*x)
            d=np.linalg.det(wave_matrix(species,w,kpar,kperp,va)*(va*w)**2)
            return [d.real,d.imag]
        solution=solve_root(residual,[guess.real,guess.imag],tol=1e-10)
        self.assertTrue(solution.success,solution.message)
        expected=complex(*solution.x)
        # Select the branch independently before giving the secant method a
        # nearby (not solved) guess: different methods can select different
        # legitimate branches from a broad initial guess.
        branch_guess=.98*expected+1e-4j
        blocks,roots,_,_=self.evaluate(species,[branch_guess],root=True,checked=checked,kpar=kpar,kperp=kperp,va=va)
        flag,w=roots[0]
        self.assertGreater(flag,0)
        self.assertLess(flag,75)
        np.testing.assert_allclose([w.real,w.imag],[expected.real,expected.imag],rtol=1e-8,atol=1e-10)
        type(self).root_error=max(self.root_error,abs(w-expected))
        if label:
            type(self).roots[label]=w
        self.assertEqual(blocks[0]['status'],0)
        residual_norm=abs(blocks[0]['det'])/np.max(abs(blocks[0]['lam']))**3
        type(self).maximum_residual=max(self.maximum_residual,float(residual_norm))
        # Match PLUME's existing scaled convergence criterion. The frequency
        # agreement above is an independent, much tighter branch check.
        self.assertLess(residual_norm,1e-8)
        for value in (w,expected):
            matrix=wave_matrix(species,value,kpar,kperp,va)
            self.assertLess(abs(np.linalg.det(matrix/np.max(abs(matrix)))),1e-8)
        return w

    def test_exact_parallel_roots(self):
        for checked in self.probes:
            for drift in (0.,-.15,.15):
                for sign in (-1,1):
                    self.check_root([P,mixed.replace_species(E,drift=drift)],.05-.001j,kpar=sign*.06,kperp=0.,
                                    checked=checked,label=f'parallel u={drift},k_sign={sign}')

    def test_oblique_roots(self):
        for species,label in (([P,E],'electrons'),([P,mixed.replace_species(E,drift=-.15)],'electron drift'),
                              ([P,E,ALPHA],'warm electrons/alphas'),
                              ([P,mixed.replace_species(E,tau=0,alpha=0),ALPHA],'warm alpha drift')):
            self.check_root(species,.05-.001j,label=label)

    def test_growing_and_damped_beam(self):
        species=[P,mixed.replace_species(E,tau=0,alpha=0,density=1.1,drift=.2/1.1),
                 mixed.replace_species(ALPHA,alpha=1.,drift=2.)]
        growing=self.check_root(species,.3+.08j,kpar=.3,kperp=.1,va=.01,label='growing warm alpha beam',checked=True)
        self.assertGreater(growing.imag,.01)
        # Oblique warm-electron Landau damping, independently predicted.
        damped=self.check_root([P,E],.05-.001j,label='damped oblique')
        self.assertLess(damped.imag,-1e-8)

    def test_species_order_invariance(self):
        base=[P,mixed.replace_species(E,density=1.1),ALPHA]
        w=.13+.03j
        for checked in self.probes:
            for kperp in (0.,.05):
                a=self.evaluate(base,[w],checked=checked,kperp=kperp)[0][0]
                reordered=[base[0],base[2],base[1]]
                b=self.evaluate(reordered,[w],checked=checked,kperp=kperp)[0][0]
                self.compare_tensor(a['chi'],b['chi'][[0,2,1]],order=True)
                self.compare_tensor(a['lam'],b['lam'],order=True)
                r1=self.check_root(base,.05-.001j,kperp=kperp,checked=checked)
                r2=self.check_root(reordered,.98*r1+1e-4j,kperp=kperp,checked=checked)
                self.assertLess(abs(r1-r2),5e-10)

    def test_poles_numerical_failure_and_recovery(self):
        for checked in self.probes:
            blocks,_,_,_=self.evaluate([P,E],[.13+.03j,0j,1+0j,-1+0j,.13+.03j],checked=checked)
            self.assertEqual([b['status'] for b in blocks],[0,2,4,4,0])
            np.testing.assert_array_equal(blocks[0]['chi'],blocks[-1]['chi'])
            for b in blocks[1:-1]:
                self.assertTrue(np.isnan(b['chi']).all())
                self.assertTrue(np.isnan(b['lam']).all())
            bad=[P,mixed.replace_species(E,tau=1e160)]
            block=self.evaluate(bad,[.13+.03j],checked=checked)[0][0]
            self.assertEqual(block['status'],5)
            self.assertTrue(np.isnan(block['chi']).all())
            parallel_bad=[P,mixed.replace_species(E,tau=1e308)]
            block=self.evaluate(parallel_bad,[.13+.03j],checked=checked,kperp=0.)[0][0]
            self.assertEqual(block['status'],5)
            self.assertTrue(np.isnan(block['lam']).all())

    def test_native_state_and_warning_only_equilibria(self):
        base=None
        for alpha in (0.,-3.,1.,1e100):
            species=[mixed.replace_species(P,alpha=alpha,drift=.1),mixed.replace_species(E,density=.8,drift=-.1)]
            blocks,_,stdout,files=self.evaluate(species,[.13+.03j,.13-.03j],native=True)
            self.assertEqual(stdout.count('WARNING: Plasma not charge neutral'),1)
            self.assertEqual(stdout.count('WARNING: Nonzero equilibrium current'),1)
            if base is not None:
                np.testing.assert_array_equal(blocks[0]['chi'],base)
            base=blocks[0]['chi']
            self.assertIn('auxiliary_first_warm_species=2',files['normalization_test.txt'])

    def test_unsupported_requests(self):
        source=mixed.input_text([P,E])
        cases=[(source.replace("'d_ref'","'rho_ref'"),'normalization=d_ref'),
               (source.replace('betap=0','betap=.1'),'betap=0'),
               (mixed.input_text([P,E],kpar=0),'nonzero kpar'),
               (mixed.input_text([P,E],extra='collision_type=1'),'collision_type=0'),
               (mixed.input_text([P,E],option=6),'FPC'),(mixed.input_text([P,E],option=7),'FPC'),
               (mixed.input_text([P,E],option=5),'scans'),(source.replace('nscan=0','nscan=1'),'scans'),
               (mixed.input_text([P,mixed.replace_species(E,tau=-.1)]),'negative tauS')]
        for text,message in cases:
            result,_=mixed.run(mixed.ROOT/'plume.e',text)
            self.assertNotEqual(result.returncode,0)
            self.assertIn(message,result.stdout+result.stderr)

    def test_map_recovers_parallel_mode(self):
        source=mixed.input_text([P,E],option=0,kperp=0.)
        source+=('&maps\nloggridw=.false.\nloggridg=.false.\nomi=.04\nomf=.09\n'
                 'gami=-.005\ngamf=.005\nnr=40\nni=10\npositive_roots=.true.\n/\n')
        result,files=mixed.run(mixed.ROOT/'plume.e',source)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        roots=[list(map(float,l.split())) for l in files['dispersion_test.roots'].splitlines() if l.strip()]
        expected=self.check_root([P,E],.05-.001j,kperp=0.)
        self.assertTrue(any(abs(complex(*r[4:6])-expected)<5e-6 for r in roots))

    def test_example_and_wrapper(self):
        sys.path.insert(0,str(mixed.ROOT))
        from linfpclib.linfpc import plume_input
        source=(mixed.ROOT/'inputs/example/mixed_cold_protons_warm_electrons.in').read_text()
        result,_=mixed.run(mixed.ROOT/'plume.e',source.replace("dataName='example'","dataName='mixed'"))
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        wrapped=plume_input('mixed')
        wrapped.set_params(0.,.05,.06,.001,2,0,1,1,False,True,normalization='d_ref')
        wrapped.make_species(*P);wrapped.make_species(*E);wrapped.make_guess(.05,-.001)
        with tempfile.TemporaryDirectory() as directory:
            path=str(Path(directory)/'input')
            wrapped.write_input(path,'test')
            result,_=mixed.run(mixed.ROOT/'plume.e',Path(path+'.in').read_text())
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
