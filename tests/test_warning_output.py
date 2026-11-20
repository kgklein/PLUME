"""Warning routing, lifecycle, fallback, and numerical output invariance."""
from pathlib import Path
import subprocess
import tempfile
import unittest

import test_mixed_dispersion as mixed

ROOT = mixed.ROOT


class WarningOutputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(['make', '-j2'], cwd=ROOT, check=True, capture_output=True)
        cls.temp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.temp.cleanup)
        directory = Path(cls.temp.name)
        probe = directory/'warnings.f90'
        probe.write_text('''program probe
  use vars, only: dataName, outputName, initialize_warnings, emit_warning, close_warnings
  implicit none
  integer :: old_unit, replacement_unit
  character(16) :: mode
  call get_command_argument(1,mode)
  call emit_warning('WARNING: before initialization')
  dataName='warnings'; outputName='probe'
  call initialize_warnings()
  if (mode=='readonly') then
     inquire(file='data/warnings/probe.warn',number=old_unit)
     close(old_unit)
     open(newunit=replacement_unit,file='data/warnings/probe.warn',status='old',action='read')
     if (replacement_unit/=old_unit) error stop 'fault injection requires recycled unit'
  endif
  call emit_warning('WARNING: first')
  call emit_warning('WARNING: second')
  call close_warnings()
end program
''')
        cls.probe = directory/'probe'
        subprocess.run(['gfortran', '-fdefault-real-8', '-I'+str(ROOT/'include'),
                        str(probe), str(ROOT/'vars.o'), '-o', str(cls.probe)], check=True)

    def test_policy_and_numerical_invariance(self):
        species = [mixed.PROTON, mixed.replace_species(mixed.COLD_E, density=.9, drift=.1)]
        baseline = None
        for verbose in (False, True):
            for echo in (False, True):
                text = mixed.input_text(species, extra=f'warnings_to_stdout=.{str(echo).lower()}.')
                text = text.replace('writeOut=.true.', f'writeOut=.{str(verbose).lower()}.')
                result, files = mixed.run(ROOT/'plume.e', text)
                self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
                warning = files['test.warn']
                self.assertEqual(warning.count('WARNING: Plasma not charge neutral'), 1)
                self.assertEqual(warning.count('WARNING: Nonzero equilibrium current'), 1)
                self.assertEqual('WARNING:' in result.stdout, echo)
                if echo:
                    self.assertEqual([line for line in result.stdout.splitlines() if 'WARNING:' in line],
                                     warning.splitlines())
                numerical = {k:v for k,v in files.items() if not k.endswith('.warn')}
                if baseline is not None:
                    self.assertEqual(numerical, baseline)
                baseline = numerical

    def test_heating_status_default_and_echo(self):
        text = (ROOT/'inputs/example/heating_corrected_cold.in').read_text()
        text = text.replace("dataName='example'", "dataName='mixed'")
        text = text.replace('option=1', 'option=-1').replace('nscan=1', 'nscan=0').replace('g_gam=-0.0001', 'g_gam=0.01')
        baseline = None
        for echo in (False, True):
            result, files = mixed.run(ROOT/'plume.e', text.replace('warnings_to_stdout=.false.',
                f'warnings_to_stdout=.{str(echo).lower()}.', 1))
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
            warnings = next(v for k,v in files.items() if k.endswith('.warn'))
            self.assertIn('WARNING: corrected heating/energy status, applicability, point, mode:', warnings)
            self.assertEqual('WARNING:' in result.stdout, echo)
            numerical = {k:v for k,v in files.items() if not k.endswith('.warn')}
            if baseline is not None:
                self.assertEqual(numerical, baseline)
            baseline = numerical

    def test_rerun_replaces_log_and_warning_free_run_is_empty(self):
        with tempfile.TemporaryDirectory() as temporary:
            cwd = Path(temporary)
            (cwd/'data/mixed').mkdir(parents=True)
            for density in (.9, .9, 1.):
                text = mixed.input_text([mixed.PROTON, mixed.replace_species(mixed.COLD_E, density=density)], option=0)
                text += '&maps loggridw=.false.,loggridg=.false.,omi=.04,omf=.08,gami=.001,gamf=.002,nr=1,ni=1 /\n'
                (cwd/'case.in').write_text(text)
                result = subprocess.run([str(ROOT/'plume.e'), 'case.in'], cwd=cwd, text=True, capture_output=True)
                self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
                warning = (cwd/'data/mixed/test.warn').read_text()
                self.assertEqual(warning.count('WARNING:'), int(density!=1.))

    def test_early_warning_and_open_failure_fallback(self):
        with tempfile.TemporaryDirectory() as temporary:
            result = subprocess.run([str(self.probe)], cwd=temporary, text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, '')
        self.assertEqual(result.stderr.count('unable to write warning log'), 1)
        for message in ('before initialization', 'first', 'second'):
            self.assertEqual(result.stderr.count('WARNING: '+message), 1)

    def test_write_failure_fallback(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)/'data/warnings'
            directory.mkdir(parents=True)
            result = subprocess.run([str(self.probe), 'readonly'], cwd=temporary, text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stderr.count('unable to write warning log'), 1)
        for message in ('first', 'second'):
            self.assertEqual(result.stderr.count('WARNING: '+message), 1)

    def test_fpc_grid_warnings(self):
        from test_normalization import input_text, run_case
        species = [(1,1,1.3,1,1,0), (.5,1836,1,-1,1,.02)]
        for option, warning in ((6, 'vparmin'), (7, 'vxmin')):
            text = input_text('rho_ref', .4, 1.3, .05, .06, species)
            text = text.replace('option=-1', f'option={option}')
            text += ('&fpc vperpmin=0, vperpmax=.5, vparmin=-.6, vparmax=.5, '
                     'vxmin=-.6, vxmax=.5, vymin=-.5, vymax=.5, vzmin=-.5, vzmax=.5, delv=.5 /\n'
                     '&guess_1\ng_om=.08\ng_gam=-.001\n/\n')
            result, files = run_case(ROOT/'plume.e', text)
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
            log = next(v for k,v in files.items() if k.endswith('.warn'))
            self.assertIn('WARNING: '+warning+' not integer multiple of delv', log)
            self.assertNotIn('WARNING:', result.stdout)

    def test_python_wrapper_control(self):
        from linfpclib.linfpc import plume_input
        obj = plume_input('warnings')
        for echo in (False, True):
            obj.set_params(1., .05, .06, .001, 2, 0, -1, 1, False, False,
                           warnings_to_stdout=echo)
            self.assertEqual(obj.params['warnings_to_stdout'], '.true.' if echo else '.false.')

    def test_fatal_errors_remain_visible(self):
        text = mixed.input_text([mixed.PROTON, mixed.COLD_E], extra="z_backend='invalid'")
        result, _ = mixed.run(ROOT/'plume.e', text)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('ERROR: z_backend', result.stdout+result.stderr)


if __name__ == '__main__':
    unittest.main()
