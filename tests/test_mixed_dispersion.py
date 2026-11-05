from faddeeva_build import compile_c,c_objects
"""Stage 2 production assembly, roots, boundaries and kinetic-to-cold limits."""
from pathlib import Path
import subprocess
import tempfile
import unittest
import sys

import numpy as np
from scipy.optimize import root as solve_root
from mixed_reference import species_tensor, wave_matrix

ROOT = Path(__file__).resolve().parents[1]
PROTON = (.4, 1., 1.3, 1., 1., 0.)
ELECTRON = (.2, 1836., .8, -1., 1., 0.)
COLD_E = (0., 1836., 0., -1., 1., 0.)
COLD_ALPHA = (0., .25, 0., .5, .05, .3)


def input_text(species, kpar=.06, kperp=.05, va=1e-3, option=-1, extra='', guess=.06-.0001j):
    lines = ['&params', "heating_convention='legacy'", "z_backend='legacy'", "bessel_mode='fixed'", "normalization='d_ref'", f'betap={species[0][0]:.17g}',
             f'kpar={kpar:.17g}', f'kperp={kperp:.17g}', f'vtp={va:.17g}',
             f'nspec={len(species)}', 'nscan=0', f'option={option}',
             'nroot_max=1', 'use_map=.false.', 'writeOut=.true.',
             "dataName='mixed'", "outputName='test'", extra, '/']
    for i, s in enumerate(species, 1):
        lines += [f'&species_{i}'] + [f'{name}={value:.17g}' for name, value in
                  zip(('tauS', 'muS', 'alphS', 'Qs', 'Ds', 'vvS'), s)] + ['/']
    lines += ['&guess_1', f'g_om={guess.real:.17g}', f'g_gam={guess.imag:.17g}', '/']
    return '\n'.join(lines)+'\n'


def run(executable, source, stdin=''):
    with tempfile.TemporaryDirectory() as temporary:
        cwd = Path(temporary)
        (cwd/'data/mixed').mkdir(parents=True)
        (cwd/'case.in').write_text(source)
        result = subprocess.run([str(executable), 'case.in'], input=stdin, cwd=cwd,
                                text=True, capture_output=True, timeout=30)
        files = {p.name: p.read_text() for p in (cwd/'data/mixed').iterdir()}
    return result, files


def parse(text, nspec):
    blocks, roots = [], []
    for line in text.splitlines():
        f = line.split()
        if not f:
            continue
        if f[0] == 'ROOT':
            roots.append((int(f[1]), complex(float(f[2]), float(f[3]))))
        if f[0] == 'DET':
            blocks.append({'status': int(f[1]), 'det': complex(float(f[2]), float(f[3])),
                           'chi': np.empty((nspec, 3, 3), complex), 'lam': np.empty((3, 3), complex)})
        elif f[0] == 'LAM':
            blocks[-1]['lam'][int(f[1])-1, int(f[2])-1] = complex(float(f[3]), float(f[4]))
        elif f[0] == 'CHI':
            blocks[-1]['chi'][int(f[1])-1, int(f[2])-1, int(f[3])-1] = complex(float(f[4]), float(f[5]))
    return blocks, roots


class MixedDispersionTests(unittest.TestCase):
    tensor_error = 0.
    root_difference = 0.
    limit_errors = {}

    @classmethod
    def setUpClass(cls):
        subprocess.run(['make', '-j2'], cwd=ROOT, check=True, capture_output=True)
        cls.temp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.temp.cleanup)
        cls.probe = Path(cls.temp.name)/'mixed_probe'
        subprocess.run(['gfortran', '-O3', '-fdefault-real-8', '-Iinclude',
                        'tests/mixed_probe.f90', 'nrtype.o', 'nrutil_trim.o', 'vars.o', 'plasma_z_faddeeva.o','plasma_z_legacy.o','plasma_z_backend.o',*c_objects(ROOT),
                        'functions.o', 'cold_susceptibility.o', 'auxiliary_normalization.o', 'kinetic_core.o', 'harmonic_selection.o', 'bessel_cache.o', 'adaptive_kinetic.o', 'cold_reference.o','species_moments.o','electromagnetic_fields.o','heating_response.o','heating_energy.o','heating_diagnostics.o','scan_support.o', 'bessel.o', 'disprels.o',
                        '-o', str(cls.probe)], cwd=ROOT, check=True, capture_output=True)
        checked = Path(cls.temp.name)/'checked'
        checked.mkdir()
        compile_c(checked)
        cls.checked_probe = checked/'mixed_probe'
        subprocess.run(['gfortran', '-O0', '-g', '-fdefault-real-8', '-fcheck=all',
                        '-ffpe-trap=invalid,zero,overflow', '-I', str(checked), '-J', str(checked),
                        *[str(ROOT/'src'/name) for name in ('nrtype.f90', 'nrutil_trim.f90',
                          'vars.f90', 'plasma_z_faddeeva.f90', 'plasma_z_legacy.f90', 'plasma_z_backend.f90', 'cold_susceptibility.f90', 'auxiliary_normalization.f90', 'kinetic_core.f90', 'harmonic_selection.f90', 'bessel_cache.f90', 'adaptive_kinetic.f90', 'cold_reference.f90', 'species_moments.f90','electromagnetic_fields.f90','heating_response.f90','heating_energy.f90','heating_diagnostics.f90', 'scan_support.f90', 'functions.f90', 'bessel.f90', 'disprels.f90')],
                        str(ROOT/'tests/mixed_probe.f90'), *c_objects(checked), '-o', str(cls.checked_probe)],
                       cwd=checked, check=True, capture_output=True)

    @classmethod
    def tearDownClass(cls):
        print(f'Mixed maximum tensor scaled error={cls.tensor_error:.3e}; '
              f'maximum root difference={cls.root_difference:.3e}; cold-limit errors={cls.limit_errors}')

    def evaluate(self, species, frequencies, root=False, **kwargs):
        source = input_text(species, **kwargs)
        stdin = ''.join(f'{w.real:.17g} {w.imag:.17g} {int(root)}\n' for w in frequencies)
        result, _ = run(self.probe, source, stdin)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        blocks, roots = parse(result.stdout, len(species))
        self.assertEqual(len(roots) if root else len(blocks), len(frequencies), result.stdout)
        return blocks, roots, result.stdout

    def test_assembly(self):
        systems = [[PROTON, COLD_E],
                   [PROTON, ELECTRON, COLD_ALPHA],
                   [PROTON, replace_species(COLD_E, drift=.07)],
                   [PROTON, replace_species(COLD_E, drift=-.07, density=.9)],
                   [PROTON, ELECTRON, replace_species(COLD_ALPHA, drift=-.3)],
                   [PROTON, replace_species(ELECTRON, density=1.12), COLD_ALPHA,
                    (0., .1, -7., -.5, .03, -.2)]]
        frequencies = [.07+.005j, .07-.005j]
        for species in systems:
            for kpar in (.06, -.06):
                blocks, _, stdout = self.evaluate(species, frequencies, kpar=kpar)
                for w, block in zip(frequencies, blocks):
                    with self.subTest(species=species, kpar=kpar, frequency=w):
                        self.assertEqual(block['status'], 0)
                        expected = np.array([species_tensor(s, w, kpar, .05, 1e-3) for s in species])
                        scale = np.maximum(1., abs(expected))
                        type(self).tensor_error = max(self.tensor_error, float(np.max(abs(block['chi']-expected)/scale)))
                        np.testing.assert_allclose(block['chi'], expected, rtol=2e-7, atol=1e-5)
                        # Independent sum of each contribution + one vacuum/wave term.
                        np.testing.assert_allclose(block['lam'], wave_matrix(species,w,kpar,.05,1e-3),
                                                   rtol=2e-7, atol=1e-5)
                        np.testing.assert_allclose(block['det'], np.linalg.det(block['lam']), rtol=2e-12)
                self.assertLessEqual(stdout.count('WARNING: Plasma not charge neutral'), 1)
                self.assertLessEqual(stdout.count('WARNING: Nonzero equilibrium current'), 1)

    def test_warm_contributions_unchanged(self):
        # Changing only electron/alpha temperature must not affect any warm chi.
        mixed = [PROTON, ELECTRON, COLD_ALPHA]
        warm = [PROTON, ELECTRON, replace_species(COLD_ALPHA, tau=.03, alpha=1.)]
        a = self.evaluate(mixed, [.07+.005j])[0][0]
        b = self.evaluate(warm, [.07+.005j])[0][0]
        np.testing.assert_array_equal(a['chi'][:2], b['chi'][:2])
        altered = [PROTON, ELECTRON, replace_species(COLD_ALPHA, alpha=-100.)]
        c = self.evaluate(altered, [.07+.005j])[0][0]
        np.testing.assert_array_equal(a['chi'], c['chi'])

    def test_roots(self):
        systems = [([PROTON, COLD_E], 0.),
                   ([PROTON, COLD_E], .01),
                   ([PROTON, ELECTRON, COLD_ALPHA], .05),
                   ([PROTON, replace_species(COLD_E, drift=-.07)], .09)]
        for species, kperp in systems:
            guess = .05-.001j
            blocks, roots, _ = self.evaluate(species, [guess], root=True, kperp=kperp)
            self.assertGreater(roots[0][0], 0)
            self.assertLess(roots[0][0], 75)
            plume_root = roots[0][1]
            def residual(x):
                w = complex(*x)
                determinant = np.linalg.det(wave_matrix(species, w, .06, kperp, 1e-3)*(1e-3*w)**2)
                return [determinant.real, determinant.imag]
            independent = solve_root(residual, [guess.real, guess.imag], tol=1e-10)
            self.assertTrue(independent.success, independent.message)
            reference = complex(*independent.x)
            self.assertLess(np.linalg.norm(residual(independent.x)), 1e-12)
            delta = abs(plume_root-reference)
            type(self).root_difference = max(self.root_difference, delta)
            np.testing.assert_allclose([plume_root.real, plume_root.imag],
                                       [reference.real, reference.imag], rtol=1e-8, atol=1e-10)
            block = blocks[0]
            self.assertEqual(block['status'], 0)
            self.assertLess(abs(block['det'])/np.max(abs(block['lam']))**3, 1e-10)
            # Oracle harmonic truncation independently checked.
            np.testing.assert_allclose(wave_matrix(species,reference,.06,kperp,1e-3,60),
                                       wave_matrix(species,reference,.06,kperp,1e-3,90), rtol=1e-12)

    def test_kinetic_to_cold_limit(self):
        # A finite-mass non-reference ion is well-conditioned across this range.
        # Keep cold electrons throughout to test the mixed exact-parallel path.
        for drift in (0., .25):
            for kperp in (0., .001, .05):
                target = (0., .25, 1., .5, .05, drift)
                cold = self.evaluate([PROTON, COLD_E, target], [.12+.03j], kperp=kperp)[0][0]['chi'][2]
                errors = []
                for tau in (.1, .01, .001, .0001, .00001):
                    warm = replace_species(target, tau=tau)
                    response = self.evaluate([PROTON, COLD_E, warm], [.12+.03j], kperp=kperp)[0][0]['chi'][2]
                    errors.append(float(np.max(abs(response-cold))/np.max(abs(cold))))
                self.assertTrue(all(a > b for a, b in zip(errors, errors[1:])), errors)
                self.assertLess(errors[-1], errors[0]/100, errors)
                self.assertLess(errors[-1], 1e-4)
                type(self).limit_errors[f'u={drift},kperp={kperp}'] = errors

    def test_exact_parallel_assembly(self):
        species = [PROTON, replace_species(COLD_E, drift=.03), replace_species(ELECTRON, drift=-.05)]
        for kpar in (-.06, .06):
            blocks, _, _ = self.evaluate(species, [.07+.005j, .07-.005j], kpar=kpar, kperp=0.)
            for w, block in zip([.07+.005j, .07-.005j], blocks):
                self.assertEqual(block['status'], 0)
                expected = np.array([species_tensor(s,w,kpar,0.,1e-3) for s in species])
                np.testing.assert_allclose(block['chi'], expected, rtol=2e-6, atol=1e-5)
                np.testing.assert_allclose(block['lam'], wave_matrix(species,w,kpar,0.,1e-3),
                                           rtol=2e-6, atol=1e-5)
                np.testing.assert_array_equal(block['chi'][:,[0,1,2,2],[2,2,0,1]],0.)
        # The oblique expressions converge to this same parallel response.
        parallel = self.evaluate(species,[.07+.005j],kperp=0.)[0][0]['chi']
        oblique = self.evaluate(species,[.07+.005j],kperp=1e-5)[0][0]['chi']
        np.testing.assert_allclose(oblique[:,[0,1,2],[0,1,2]], parallel[:,[0,1,2],[0,1,2]],rtol=2e-6)

    def test_checked_mixed_kernel(self):
        for kperp in (0., .05):
            species = [PROTON, ELECTRON, COLD_ALPHA]
            source = input_text(species, kperp=kperp)
            # A singular point between valid samples must not cause a thermal
            # divide-by-zero or leave stale data, even with floating-point traps.
            result, _ = run(self.checked_probe, source, '.07 .005 0\n0 0 0\n.07 -.005 0\n')
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
            blocks, _ = parse(result.stdout, len(species))
            self.assertEqual([b['status'] for b in blocks], [0, 2, 0])
            for w, block in zip([.07+.005j, .07-.005j], [blocks[0], blocks[2]]):
                expected = np.array([species_tensor(s,w,.06,kperp,1e-3) for s in species])
                np.testing.assert_allclose(block['chi'], expected, rtol=2e-7, atol=1e-5)

    def test_equilibrium_warnings_once(self):
        species = [PROTON, replace_species(COLD_E, density=.9, drift=.1)]
        _, _, stdout = self.evaluate(species, [.07+.005j, .08+.005j, .09+.005j])
        self.assertEqual(stdout.count('WARNING: Plasma not charge neutral'), 1)
        self.assertEqual(stdout.count('WARNING: Nonzero equilibrium current'), 1)

    def test_example_and_python_wrapper(self):
        sys.path.insert(0,str(ROOT))
        from linfpclib.linfpc import plume_input
        example = (ROOT/'inputs/example/mixed_warm_protons_cold_electrons.in').read_text()
        source = example.replace("dataName='example'", "dataName='mixed'")
        result, _ = run(ROOT/'plume.e', source)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        self.assertIn('Dispersion Solutions', result.stdout)
        wrapped = plume_input('mixed')
        wrapped.set_params(.4,.05,.06,1e-3,2,0,-1,1,False,True,normalization='d_ref',Kn=10.)
        wrapped.make_species(*PROTON)
        wrapped.make_species(*COLD_E)
        wrapped.make_guess(.07, .005)
        self.assertEqual(wrapped.get_nuns(wrapped.species[1]), 0.)
        with tempfile.TemporaryDirectory() as directory:
            name = str(Path(directory)/'mixed')
            wrapped.write_input(name,'test')
            result, _ = run(ROOT/'plume.e', Path(name+'.in').read_text())
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        wrapped.params['collision_type'] = 1
        with self.assertRaisesRegex(ValueError, 'collision_type=0'):
            wrapped.check_input_for_errors()
        wrapped.params['collision_type'] = 0
        wrapped.params['nscan'] = 1
        with self.assertRaisesRegex(ValueError, 'scans'):
            wrapped.check_input_for_errors()

    def test_poles_and_recovery(self):
        species = [PROTON, replace_species(COLD_ALPHA, drift=.5)]
        frequencies = [.07+.005j, 0j, .03+0j, .53+0j, -.47+0j, .07+.005j]
        blocks, _, _ = self.evaluate(species, frequencies)
        self.assertEqual([b['status'] for b in blocks], [0,2,3,4,4,0])
        for block in blocks[1:-1]:
            self.assertTrue(np.isnan(block['det']))
            self.assertTrue(np.isnan(block['lam']).all())
            self.assertTrue(np.isnan(block['chi']).all())
        np.testing.assert_array_equal(blocks[0]['chi'], blocks[-1]['chi'])
        # No cold pole is returned as a root.
        blocks, roots, _ = self.evaluate(species, [0j], root=True)
        self.assertEqual(roots[0][0], -1)
        self.assertTrue(np.isnan(roots[0][1]))

    def test_public_paths_and_boundaries(self):
        species = [PROTON, COLD_E]
        for option in (-1, 1):
            result, _ = run(ROOT/'plume.e', input_text(species, option=option))
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
            self.assertNotIn('Ps =', result.stdout)
            self.assertIn('determinant only' if option == -1 else 'Dispersion Solutions', result.stdout)
        source = input_text(species)
        bad_cases = [(source.replace("'d_ref'", "'rho_ref'"), 'require normalization=d_ref'),
                     (input_text([replace_species(PROTON,tau=0),ELECTRON],kpar=0), 'nonzero kpar'),
                     (input_text([PROTON,replace_species(COLD_E,tau=-1)]), 'negative tauS'),
                     (input_text(species,extra='collision_type=1'), 'collision_type=0'),
                     (input_text(species,option=6), 'FPC'), (input_text(species,option=7), 'FPC'),
                     (source.replace('nscan=0','nscan=1'), 'scans'),
                     (input_text(species,option=5), 'scans'),
                     (input_text(species,option=2), 'option=-1, 0, or 1'),
                     (input_text(species,option=4), 'option=-1, 0, or 1'),
                     (input_text(species,kpar=0), 'nonzero kpar'),
                     (input_text([PROTON,replace_species(COLD_E,density=-1)]), 'nonnegative density'),
                     (input_text(species,option=1,guess=0j), 'zero-frequency pole')]
        for bad, message in bad_cases:
            result, _ = run(ROOT/'plume.e', bad)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(message, result.stdout+result.stderr)

    def test_map_with_singular_points(self):
        # Zero frequency, zero Doppler shift and ion cyclotron poles lie on grid.
        source = input_text([PROTON, (0.,1.,0.,1.,.1,.5)], option=0, kpar=.5)
        source += ('&maps\nloggridw=.false.\nloggridg=.false.\nomi=-.75\nomf=1.25\n'
                   'gami=-.05\ngamf=.05\nnr=8\nni=2\npositive_roots=.false.\n/\n')
        result, files = run(ROOT/'plume.e', source)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        self.assertIn('singular/nonfinite mixed map points', result.stdout)
        rows = [list(map(float,l.split())) for l in files['dispersion_test.map'].splitlines() if l.strip()]
        self.assertTrue(any(np.isnan(r[2:]).all() for r in rows))
        self.assertTrue(any(np.isfinite(r).all() for r in rows))
        for line in files['dispersion_test.roots'].splitlines():
            self.assertTrue(np.isfinite(list(map(float,line.split()))).all())


def replace_species(species, **updates):
    names = ('tau', 'mu', 'alpha', 'charge', 'density', 'drift')
    return tuple(updates.get(name, value) for name, value in zip(names, species))


if __name__ == '__main__':
    unittest.main()
