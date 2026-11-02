"""Finite-temperature equivalence tests for PLUME's two input conventions."""

import math
import re
import subprocess
import tempfile
import unittest
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
PROBE = ROOT / 'tests' / 'normalization_probe.f90'


def input_text(mode, beta, alpha_ref, kperp_d, kpar_d, species, scan=None):
    scale = math.sqrt(beta * alpha_ref)
    if mode == 'rho_ref':
        kperp, kpar, vtp = kperp_d * scale, kpar_d * scale, 1e-4 * math.sqrt(beta)
    else:
        kperp, kpar, vtp = kperp_d, kpar_d, 1e-4
    option = 1 if scan else -1
    lines = [
        '&params', f"normalization='{mode}'", f'betap={beta:.17g}',
        f'kperp={kperp:.17g}', f'kpar={kpar:.17g}', f'vtp={vtp:.17g}',
        f'nspec={len(species)}', f'nscan={1 if scan else 0}', f'option={option}',
        'nroot_max=1', 'use_map=.false.', 'writeOut=.false.',
        "dataName='normalization'", "outputName='equivalence'", '/',
    ]
    for idx, (tau, mu, alpha, charge, density, drift) in enumerate(species, 1):
        tau_input = 1 / tau if mode == 'rho_ref' else beta * tau
        lines += [f'&species_{idx}', f'tauS={tau_input:.17g}', f'muS={mu:.17g}',
                  f'alphS={alpha:.17g}', f'Qs={charge:.17g}',
                  f'Ds={density:.17g}', f'vvS={drift:.17g}', '/']
    if scan:
        axis, end_d = scan
        end = end_d * (scale if mode == 'rho_ref' else 1)
        lines += ['&scan_input_1', f'scan_type={1 if axis == "kpar" else 0}',
                  'scan_style=0', f'swi={kpar if axis == "kpar" else kperp:.17g}',
                  f'swf={end:.17g}', 'swlog=.false.', 'ns=2', 'nres=1',
                  'heating=.true.', 'eigen=.true.', 'tensor=.true.', '/']
        lines += ['&guess_1', 'g_om=0.01', 'g_gam=-0.001', '/']
    return '\n'.join(lines) + '\n'


def run_case(executable, source):
    with tempfile.TemporaryDirectory() as directory:
        cwd = Path(directory)
        (cwd / 'data' / 'normalization').mkdir(parents=True)
        (cwd / 'case.in').write_text(source)
        result = subprocess.run([str(executable), 'case.in'], cwd=cwd,
                                text=True, capture_output=True, timeout=90)
        files = {p.name: p.read_text() for p in (cwd / 'data' / 'normalization').iterdir()}
        return result, files


def numeric_records(text):
    records = {}
    for line in text.splitlines():
        fields = line.split()
        if fields and fields[0] in {'DET', 'CHI', 'FIELD', 'HEAT', 'SPLIT', 'DENS', 'FLOW', 'ROOT'}:
            values = [float(x.replace('D', 'E')) for x in fields[1:]]
            indices = {'CHI': 3, 'FIELD': 1, 'HEAT': 1, 'SPLIT': 2, 'DENS': 1, 'FLOW': 2, 'DET': 0, 'ROOT': 0}[fields[0]]
            key = (fields[0], *[int(x) for x in values[:indices]])
            records[key] = values[indices:]
    return records


class NormalizationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(['make', '-j2'], cwd=ROOT, check=True, capture_output=True)
        cls.tempdir = tempfile.TemporaryDirectory()
        cls.probe = Path(cls.tempdir.name) / 'probe'
        subprocess.run(['gfortran', '-O3', '-fdefault-real-8', '-Iinclude',
                        '-o', str(cls.probe), str(PROBE), 'nrtype.o', 'nrutil_trim.o',
                        'vars.o', 'functions.o', 'bessel.o', 'disprels.o', 'fpc.o'],
                       cwd=ROOT, check=True, capture_output=True)

    @classmethod
    def tearDownClass(cls):
        cls.tempdir.cleanup()

    def compare_records(self, left, right):
        self.assertEqual(left.keys(), right.keys())
        for key in left:
            for a, b in zip(left[key], right[key]):
                self.assertTrue(math.isfinite(a) and math.isfinite(b), key)
                self.assertAlmostEqual(a, b, delta=1e-10 + 5e-8 * max(abs(a), abs(b)), msg=str(key))

    def test_equivalent_plasmas(self):
        proton_electron = [(1, 1, 1, 1, 1, 0), (0.5, 1836, 1, -1, 1, 0)]
        cases = [
            ('isotropic', 1, 1, proton_electron),
            ('anisotropic', 1, 1.7, [(1, 1, 1.7, 1, 1, 0), (0.5, 1836, 0.8, -1, 1, 0)]),
            ('electron drift', 1, 1, [(1, 1, 1, 1, 1, 0), (0.5, 1836, 1, -1, 1, 0.03)]),
            ('alpha drift', 0.3, 1, [(1, 1, 1, 1, 1, 0), (0.7, 0.25, 1.2, 0.5, 0.05, 0.08),
                                    (0.5, 1836, 1, -1, 1.1, 0.003636363636)]),
            ('low beta', 0.02, 1.3, [(1, 1, 1.3, 1, 1, 0), (0.5, 1836, 1, -1, 1, 0.02)]),
        ]
        for name, beta, alpha, species in cases:
            with self.subTest(name=name):
                sources = [input_text(mode, beta, alpha, .06, .08, species) for mode in ('rho_ref', 'd_ref')]
                left = run_case(self.probe, sources[0])[0]
                right = run_case(self.probe, sources[1])[0]
                self.assertEqual(left.returncode, 0, left.stderr + left.stdout)
                self.assertEqual(right.returncode, 0, right.stderr + right.stdout)
                self.compare_records(numeric_records(left.stdout), numeric_records(right.stdout))

    def test_wavevector_scans(self):
        species = [(1, 1, 1.4, 1, 1, 0), (0.5, 1836, 1, -1, 1, 0.02)]
        for axis in ('kpar', 'kperp'):
            with self.subTest(axis=axis):
                outputs = []
                for mode in ('rho_ref', 'd_ref'):
                    result, files = run_case(ROOT / 'plume.e', input_text(mode, .4, 1.4, .05, .06, species, (axis, .09)))
                    self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
                    mode_name = next(name for name in files if name.endswith('.mode1') and '_tensor_' not in name)
                    rows = [[float(v.replace('D', 'E')) for v in line.split()] for line in files[mode_name].splitlines() if line.strip()]
                    self.assertEqual(len(rows), 3)
                    outputs.append(rows)
                scale = math.sqrt(.4 * 1.4)
                for rho, inertial in zip(*outputs):
                    self.assertAlmostEqual(rho[0] / scale, inertial[0], delta=1e-4)
                    self.assertAlmostEqual(rho[1] / scale, inertial[1], delta=1e-4)
                    for a, b in zip(rho[4:-6 * len(species)], inertial[4:-6 * len(species)]):
                        self.assertAlmostEqual(a, b, delta=1e-4 * max(1, abs(a), abs(b)))
                    for index in range(len(species)):
                        start = -6 * len(species) + 6 * index
                        self.assertAlmostEqual(.4 / rho[start], inertial[start], delta=1e-4)
                        for a, b in zip(rho[start + 1:start + 6], inertial[start + 1:start + 6]):
                            self.assertAlmostEqual(a, b, delta=1e-4 * max(1, abs(a), abs(b)))

    def test_validation(self):
        species = [(1, 1, 1, 1, 1, 0), (0.5, 1836, 1, -1, 1, 0)]
        valid = input_text('d_ref', 1, 1, .05, .06, species)
        for bad, message in [(valid.replace("'d_ref'", "'wrong'"), 'normalization'),
                             (valid.replace('tauS=1', 'tauS=0', 1), 'cold-fluid'),
                             (valid.replace('tauS=1', 'tauS=0.9', 1), 'tauS = betap')]:
            result, _ = run_case(self.probe, bad)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(message, result.stdout + result.stderr)

    def test_reference_beta_and_temperature_scans(self):
        species = [(1, 1, 1.4, 1, 1, 0), (0.5, 1836, 1, -1, 1, 0)]
        base = input_text('d_ref', .4, 1.4, .05, .06, species, ('kpar', .09))
        scan_block = re.compile(r'&scan_input_1\n.*?\n/\n', re.DOTALL)
        outputs = []
        for style, kind in ((0, 2), (1, 0)):
            block = (f'&scan_input_1\nscan_type={kind}\nscan_style={style}\n'
                     'swi=0.4\nswf=0.6\nswlog=.false.\nns=2\nnres=1\n'
                     'heating=.true.\neigen=.true.\ntensor=.true.\n/\n')
            source = scan_block.sub(block, base)
            result, files = run_case(ROOT / 'plume.e', source)
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            filename = next(name for name in files if name.endswith('.mode1') and '_tensor_' not in name)
            rows = [[float(value) for value in line.split()] for line in files[filename].splitlines() if line.strip()]
            for row in rows:
                self.assertAlmostEqual(row[2], row[-6 * len(species)], delta=1e-6)
            outputs.append(rows)
        for beta_row, tau_row in zip(*outputs):
            for a, b in zip(beta_row, tau_row):
                self.assertAlmostEqual(a, b, delta=1e-4 * max(1, abs(a), abs(b)))

        first = scan_block.sub('&scan_input_1\nscan_type=2\nscan_style=0\nswi=0.4\nswf=0.6\n'
                               'swlog=.false.\nns=1\nnres=1\nheating=.false.\neigen=.false.\ntensor=.false.\n/\n', base)
        second = '&scan_input_2\nscan_type=0\nscan_style=1\nswi=0.4\nswf=0.6\n'
        second += 'swlog=.false.\nns=1\nnres=1\nheating=.false.\neigen=.false.\ntensor=.false.\n/\n'
        source = first.replace('nscan=1', 'nscan=2').replace('option=1', 'option=2') + second
        result, _ = run_case(ROOT / 'plume.e', source)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('cannot scan betap and species_1 tauS', result.stdout + result.stderr)

    def test_two_dimensional_wavevector_scan(self):
        species = [(1, 1, 1.2, 1, 1, 0), (0.5, 1836, 1, -1, 1, 0.01)]
        outputs = []
        for mode in ('rho_ref', 'd_ref'):
            source = input_text(mode, .5, 1.2, .05, .06, species, ('kpar', .08))
            scale = math.sqrt(.5 * 1.2) if mode == 'rho_ref' else 1
            source = source.replace('nscan=1', 'nscan=2').replace('option=1', 'option=2')
            source += ('&scan_input_2\nscan_type=0\nscan_style=0\n'
                       f'swi={.05 * scale:.17g}\nswf={.07 * scale:.17g}\n'
                       'swlog=.false.\nns=1\nnres=1\nheating=.true.\n'
                       'eigen=.true.\ntensor=.true.\n/\n')
            result, files = run_case(ROOT / 'plume.e', source)
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            name = next(name for name in files if name.endswith('.mode1') and '_tensor_' not in name)
            rows = [[float(value) for value in line.split()] for line in files[name].splitlines() if line.strip()]
            outputs.append(rows)
        rho_rows, d_rows = outputs
        self.assertEqual(len(rho_rows), len(d_rows))
        for rho, inertial in zip(rho_rows, d_rows):
            self.assertAlmostEqual(rho[0] / math.sqrt(.6), inertial[0], delta=1e-4)
            self.assertAlmostEqual(rho[1] / math.sqrt(.6), inertial[1], delta=1e-4)
            for a, b in zip(rho[4:-12], inertial[4:-12]):
                self.assertAlmostEqual(a, b, delta=1e-4 * max(1, abs(a), abs(b)))

    def test_matched_example_files(self):
        rows = []
        for mode in ('rho_ref', 'd_ref'):
            source = (ROOT / 'inputs' / 'example' / f'normalization_{mode}.in').read_text()
            with tempfile.TemporaryDirectory() as directory:
                cwd = Path(directory)
                (cwd / 'data' / 'example').mkdir(parents=True)
                (cwd / 'case.in').write_text(source)
                result = subprocess.run([str(ROOT / 'plume.e'), 'case.in'], cwd=cwd,
                                        text=True, capture_output=True, timeout=90)
                self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
                if mode == 'd_ref':
                    self.assertTrue((cwd / 'data' / 'example' / 'normalization_normalization_d_ref.txt').exists())
                name = next((cwd / 'data' / 'example').glob('*.mode1'))
                rows.append([[float(value) for value in line.split()]
                             for line in name.read_text().splitlines() if line.strip()])
        self.assertEqual(len(rows[0]), len(rows[1]))
        for rho, inertial in zip(*rows):
            self.assertAlmostEqual(rho[0] / math.sqrt(.6), inertial[0], delta=1e-4)
            self.assertAlmostEqual(rho[1] / math.sqrt(.6), inertial[1], delta=1e-4)
            for a, b in zip(rho[4:-12], inertial[4:-12]):
                self.assertAlmostEqual(a, b, delta=1e-4 * max(1, abs(a), abs(b)))

    def test_existing_inputs_match_unmodified_solver(self):
        with tempfile.TemporaryDirectory() as directory:
            baseline = Path(directory)
            tracked = subprocess.run(['git', 'ls-files', 'src', 'Makefile'], cwd=ROOT,
                                     text=True, capture_output=True, check=True).stdout.splitlines()
            for relative in tracked:
                content = subprocess.run(['git', 'show', f'HEAD:{relative}'], cwd=ROOT,
                                         capture_output=True, check=True).stdout
                path = baseline / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(content)
            (baseline / 'include').mkdir(exist_ok=True)
            subprocess.run(['make', '-j2'], cwd=baseline, check=True, capture_output=True)
            old_probe = baseline / 'probe'
            subprocess.run(['gfortran', '-O3', '-fdefault-real-8', '-Iinclude',
                            '-o', str(old_probe), str(PROBE), 'nrtype.o', 'nrutil_trim.o',
                            'vars.o', 'functions.o', 'bessel.o', 'disprels.o', 'fpc.o'],
                           cwd=baseline, check=True, capture_output=True)
            for relative in ('inputs/example/example_map_par.in', 'tests/drift_cbe_test.in'):
                source = (ROOT / relative).read_text()
                old = run_case(old_probe, source)[0]
                new = run_case(self.probe, source)[0]
                self.assertEqual(old.returncode, 0, old.stderr + old.stdout)
                self.assertEqual(new.returncode, 0, new.stderr + new.stdout)
                self.assertEqual(numeric_records(old.stdout), numeric_records(new.stdout))

    def test_weak_collision_conversion(self):
        species = [(1, 1, 1.3, 1, 1, 0), (0.5, 1836, 1, -1, 1, 0.01)]
        beta = .4
        scale = math.sqrt(beta * 1.3)
        sources = []
        for mode, kn in (('rho_ref', 100.0), ('d_ref', 100.0 * scale)):
            source = input_text(mode, beta, 1.3, .05, .06, species)
            source = source.replace('nroot_max=1', f'nroot_max=1\ncollision_type=1\nKn={kn:.17g}')
            sources.append(source)
        results = [run_case(self.probe, source)[0] for source in sources]
        for result in results:
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        self.compare_records(numeric_records(results[0].stdout), numeric_records(results[1].stdout))

    def test_velocity_and_species_temperature_scans(self):
        beta, alpha = .4, 1.2
        species = [(1, 1, alpha, 1, 1, 0), (0.5, 1836, 1, -1, 1, 0.01)]
        for kind in ('vtp', 'electron_tau'):
            with self.subTest(kind=kind):
                outputs = []
                for mode in ('rho_ref', 'd_ref'):
                    source = input_text(mode, beta, alpha, .05, .06, species, ('kpar', .08))
                    if kind == 'vtp':
                        style, scan_type = 0, 3
                        initial = 1e-4 * (math.sqrt(beta) if mode == 'rho_ref' else 1)
                        final = 1.5 * initial
                    else:
                        style, scan_type = 2, 0
                        initial = 2 if mode == 'rho_ref' else beta * .5
                        final = 1.25 if mode == 'rho_ref' else beta * .8
                    block = (f'&scan_input_1\nscan_type={scan_type}\nscan_style={style}\n'
                             f'swi={initial:.17g}\nswf={final:.17g}\nswlog=.false.\n'
                             'ns=1\nnres=1\nheating=.true.\neigen=.true.\ntensor=.true.\n/\n')
                    source = re.sub(r'&scan_input_1\n.*?\n/\n', block, source, flags=re.DOTALL)
                    result, files = run_case(ROOT / 'plume.e', source)
                    self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
                    name = next(name for name in files if name.endswith('.mode1') and '_tensor_' not in name)
                    outputs.append([[float(value) for value in line.split()]
                                    for line in files[name].splitlines() if line.strip()])
                for rho, inertial in zip(*outputs):
                    for a, b in zip(rho[4:-12], inertial[4:-12]):
                        self.assertAlmostEqual(a, b, delta=1e-4 * max(1, abs(a), abs(b)))

    def test_python_wrapper_round_trip(self):
        from linfpclib.linfpc import plume_input

        wrapped = plume_input('normalization')
        wrapped.set_params(.4, .05, .06, 1e-4, 2, 0, -1, 1, False, False,
                           normalization='d_ref')
        wrapped.make_species(.4, 1, 1.5, 1, 1, 0)
        wrapped.make_species(.2, 1836, 1, -1, 1, .02)
        with tempfile.TemporaryDirectory() as directory:
            name = str(Path(directory) / 'wrapped')
            wrapped.write_input(name, 'wrapped')
            self.assertIn("normalization='d_ref'", Path(name + '.in').read_text())
            loaded = plume_input('normalization')
            loaded.load_from_file(name + '.in')
            self.assertEqual(loaded.params['normalization'], 'd_ref')
            self.assertEqual(loaded.species[0]['tauS'], .4)

    def test_frequency_map_scan(self):
        species = [(1, 1, 1.3, 1, 1, 0), (0.5, 1836, 1, -1, 1, 0)]
        maps = []
        for mode in ('rho_ref', 'd_ref'):
            source = input_text(mode, .4, 1.3, .05, .06, species, ('kpar', .08))
            source = source.replace('option=1', 'option=4').replace('use_map=.false.', 'use_map=.true.')
            source += ('&maps\nloggridw=.false.\nloggridg=.false.\nomi=-0.2\nomf=0.2\n'
                       'gami=-0.1\ngamf=0.02\nnr=8\nni=8\npositive_roots=.false.\n/\n')
            result, files = run_case(ROOT / 'plume.e', source)
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
            maps.append(sorted((name, content) for name, content in files.items() if name.endswith('.map')))
        self.assertEqual(len(maps[0]), len(maps[1]))
        for (_, rho), (_, inertial) in zip(*maps):
            a = [[float(value) for value in line.split()] for line in rho.splitlines() if line.strip()]
            b = [[float(value) for value in line.split()] for line in inertial.splitlines() if line.strip()]
            self.assertEqual(len(a), len(b))
            for left, right in zip(a, b):
                for x, y in zip(left, right):
                    if math.isfinite(x) and math.isfinite(y):
                        self.assertAlmostEqual(x, y, delta=1e-3 * max(1, abs(x), abs(y)))


if __name__ == '__main__':
    unittest.main()
