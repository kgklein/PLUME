from faddeeva_build import compile_c,c_objects
"""Stage 3B.1 independent prototype, with current production isolation checks."""
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np
from auxiliary_kernel import experimental_kernel
from mixed_reference import species_tensor
from test_mixed_dispersion import ROOT, input_text, run, replace_species

P = (.4,1.,1.,1.,1.,0.)
E = (.2,1836.,1.,-1.,1.,0.)
ALPHA = (.03,.25,1.2,.5,.05,.3)
COLD_P = (0.,1.,float('nan'),1.,1.,0.)
COLD_E = (0.,1836.,float('nan'),-1.,1.,0.)
SCALES = (.01,.1,1.,10.)


class AuxiliaryNormalizationTests(unittest.TestCase):
    errors = {'scale':0., 'production':0., 'oracle':0., 'extreme':0.}
    conditioning = {}

    @classmethod
    def setUpClass(cls):
        subprocess.run(['make','-j2'],cwd=ROOT,check=True,capture_output=True)
        cls.temp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.temp.cleanup)
        cls.probes = {}
        modules = ('nrtype','nrutil_trim','vars','plasma_z_faddeeva','plasma_z_legacy','plasma_z_backend','cold_susceptibility','auxiliary_normalization','kinetic_core','harmonic_selection','bessel_cache','adaptive_kinetic','cold_reference','species_moments','electromagnetic_fields','heating_response','heating_energy','heating_diagnostics','scan_support','functions','bessel','disprels')
        for checked in (False,True):
            directory = Path(cls.temp.name)/str(checked)
            directory.mkdir()
            (directory/'kernel.f90').write_text(experimental_kernel())
            base = ['gfortran','-fdefault-real-8','-I',str(directory),'-J',str(directory)]
            if checked:
                compile_c(directory)
                flags = ['-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow']
                command = base+flags+['-c',*[str(ROOT/'src'/f'{m}.f90') for m in modules]]
                subprocess.run(command,cwd=directory,check=True,capture_output=True)
                objects = [str(directory/f'{m}.o') for m in modules]
            else:
                flags = ['-O3','-funroll-loops','-ffast-math']
                base += ['-I',str(ROOT/'include')]
                objects = [str(ROOT/f'{m}.o') for m in modules]
            # The adapter's IEEE checks remain meaningful under production flags.
            subprocess.run(base+flags+['-fno-fast-math','-c',str(ROOT/'src/auxiliary_normalization.f90')],
                           cwd=directory,check=True,capture_output=True)
            executable = directory/'probe'
            subprocess.run(base+flags+[str(directory/'kernel.f90'),str(ROOT/'tests/auxiliary_probe.f90'),
                           *objects,*c_objects(directory if checked else ROOT),'-o',str(executable)],
                           cwd=directory,check=True,capture_output=True)
            cls.probes[checked] = executable

    @classmethod
    def tearDownClass(cls):
        print('Auxiliary maximum scaled element errors:',cls.errors)
        print('Auxiliary conditioning:',cls.conditioning)

    def run_probe(self,species,requests,*,checked=False,kpar=.06,kperp=.05,va=.001):
        # A process has one fixed physical state; requests vary only auxiliary
        # scale/frequency, or select the unmodified production baseline.
        lines=[f'{len(species)} {kpar:.17g} {kperp:.17g} {va:.17g}']
        lines += [' '.join(f'{v:.17g}' for v in s) for s in species]
        lines += [f'{w.real:.17g} {w.imag:.17g} {b:.17g} {mode}' for w,b,mode in requests]
        return subprocess.run([str(self.probes[checked])],input='\n'.join(lines)+'\n',
                              text=True,capture_output=True,timeout=30)

    def evaluate(self,species,requests,**kwargs):
        result=self.run_probe(species,requests,**kwargs)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        blocks=[]
        for line in result.stdout.splitlines():
            fields=line.split()
            if fields[0]=='STATUS':
                blocks.append({'status':int(fields[1]),'chi':{}})
            elif fields[0]=='CHI':
                i,j,k=map(int,fields[1:4])
                blocks[-1]['chi'].setdefault(i,np.empty((3,3),complex))[j-1,k-1]=complex(*map(float,fields[4:6]))
        self.assertEqual(len(blocks),len(requests),result.stdout)
        for block in blocks:
            if block['status']==0:
                self.assertEqual(set(block['chi']),{i for i,s in enumerate(species,1) if s[0]>0})
        return blocks

    def compare(self,actual,expected,metric,rtol=5e-9,atol=1e-10):
        scaled=float(np.max(abs(actual-expected)/np.maximum(1.,abs(expected))))
        type(self).errors[metric]=max(self.errors[metric],scaled)
        np.testing.assert_allclose(actual,expected,rtol=rtol,atol=atol)
        return scaled

    def test_warm_scale_invariance_and_production(self):
        systems=[ [P,E],
                  [replace_species(P,alpha=1.7),replace_species(E,alpha=.7)],
                  [P,replace_species(E,drift=.2)],
                  [P,replace_species(E,density=1.1),ALPHA],
                  [replace_species(P,tau=.02,alpha=1.3),replace_species(E,tau=2.,drift=-.1),
                   replace_species(ALPHA,tau=.002,drift=-.2),(.2,.1,.8,-.5,.03,.1)] ]
        for checked in self.probes:
            for species in systems:
                for kpar in (.06,-.06):
                    for va in (.001,.05):
                        frequencies=[.13+.03j,.13-.03j]
                        requests=[(w,b,0) for w in frequencies for b in SCALES]
                        requests += [(w,1.,1) for w in frequencies]
                        blocks=self.evaluate(species,requests,checked=checked,kpar=kpar,va=va)
                        for block in blocks:
                            self.assertEqual(block['status'],0)
                        for f,w in enumerate(frequencies):
                            base=blocks[4*f]['chi']
                            production=blocks[8+f]['chi']
                            for i in base:
                                self.compare(base[i],production[i],'production')
                                for j in range(1,4):
                                    self.compare(blocks[4*f+j]['chi'][i],base[i],'scale')

    def test_original_thermal_normalization(self):
        # Explicit normalization anisotropy equals physical alpha_R only here.
        # This tests the extraction against the actual routine at the exact
        # original internal state, plus the independent isotropic auxiliary scale.
        species=[replace_species(P,alpha=1.7),replace_species(E,alpha=.7,drift=.1),ALPHA]
        w=.13+.03j
        for checked in self.probes:
            requests=[(w,.4,1),(w,.4,2),(w,.4,0),(w,.4,1)]
            blocks=self.evaluate(species,requests,checked=checked)
            for i in blocks[0]['chi']:
                # Separate compilation can change roundoff; no expressions change.
                self.compare(blocks[1]['chi'][i],blocks[0]['chi'][i],'production')
                self.compare(blocks[2]['chi'][i],blocks[0]['chi'][i],'production')
                # Experimental calls do not contaminate the production cache.
                np.testing.assert_array_equal(blocks[3]['chi'][i],blocks[0]['chi'][i])

    def test_cold_reference_against_independent_stix(self):
        systems=[[COLD_P,E],
                 [COLD_P,replace_species(E,drift=.15)],
                 [COLD_P,replace_species(E,alpha=.65,drift=-.15)],
                 [COLD_P,replace_species(E,density=1.1),ALPHA],
                 [COLD_P,COLD_E,replace_species(ALPHA,drift=-.3)]]
        for checked in self.probes:
            for species in systems:
                for kpar in (.06,-.06):
                    for va in (.001,.05):
                        for w in (.13+.03j,.13-.03j):
                            scales=SCALES+(next(s[0] for s in species if s[0]>0),)
                            blocks=self.evaluate(species,[(w,b,0) for b in scales],checked=checked,kpar=kpar,va=va)
                            for block in blocks:
                                self.assertEqual(block['status'],0)
                                for i,actual in block['chi'].items():
                                    reference=species_tensor(species[i-1],w,kpar,.05,va)
                                    self.compare(actual,reference,'oracle',rtol=2e-7,atol=1e-5)
                                    self.compare(actual,blocks[0]['chi'][i],'scale')

    def test_cold_reference_anisotropy_is_unused(self):
        w=.13+.03j
        baseline=None
        for placeholder in (float('nan'),0.,-3.,1e100):
            species=[replace_species(COLD_P,alpha=placeholder),replace_species(E,alpha=.6),ALPHA]
            block=self.evaluate(species,[(w,.1,0)],checked=True)[0]
            self.assertEqual(block['status'],0)
            if baseline is not None:
                for i in block['chi']:
                    np.testing.assert_array_equal(block['chi'][i],baseline[i])
            baseline=block['chi']

    def test_widely_separated_species_temperatures(self):
        # Five decades of native temperature, without a cold-species threshold.
        species=[replace_species(P,tau=.001,alpha=1.2),replace_species(E,tau=100.,alpha=.8),
                 replace_species(ALPHA,tau=.01),(.1,.1,1.1,1.,.02,-.1)]
        w=.3+.05j
        for checked in self.probes:
            blocks=self.evaluate(species,[(w,b,0) for b in SCALES]+[(w,1.,1)],checked=checked,
                                 kpar=.03,kperp=.01)
            for block in blocks:
                self.assertEqual(block['status'],0)
                for i in block['chi']:
                    self.compare(block['chi'][i],blocks[-1]['chi'][i],'production')

    def test_large_auxiliary_range(self):
        species=[COLD_P,replace_species(E,alpha=.8,drift=.1),ALPHA]
        w=.13+.03j
        scales=[1e-70,1e-40,1e-12,1e-6,.01,.1,1.,10.,1e6,1e12,1e40,1e70]
        for checked in self.probes:
            blocks=self.evaluate(species,[(w,b,0) for b in scales],checked=checked)
            for b,block in zip(scales,blocks):
                self.assertEqual(block['status'],0,(checked,b))
                for i in block['chi']:
                    error=self.compare(block['chi'][i],blocks[6]['chi'][i],'extreme')
                    key=f'b={b:g}'
                    type(self).conditioning[key]=max(self.conditioning.get(key,0.),error)

    def test_extreme_conditioning_and_recovery(self):
        # These scales are deliberately beyond the safe window for some
        # compiler arithmetic paths. Do not hard-code a universal threshold.
        species=[COLD_P,replace_species(E,alpha=.8,drift=.1),ALPHA]
        w=.13+.03j
        for checked in self.probes:
            baseline=self.evaluate(species,[(w,1.,0)],checked=checked)[0]
            for b in (1e-80,1e80,1e-160,1e160):
                result=self.run_probe(species,[(w,b,0),(w,1.,0)],checked=checked)
                key=f'extreme checked={checked},b={b:g}'
                if result.returncode:
                    # Trapping builds may stop on unrepresentable intermediate
                    # products. This is measured, not regularized or hidden.
                    self.assertTrue(checked,result.stderr)
                    self.assertIn('SIGFPE',result.stderr)
                    type(self).conditioning[key]='floating-point trap'
                    continue
                blocks=self.evaluate(species,[(w,b,0),(w,1.,0)],checked=checked)
                if blocks[0]['status']==0:
                    error=max(self.compare(t,baseline['chi'][i],'extreme')
                              for i,t in blocks[0]['chi'].items())
                    type(self).conditioning[key]=error
                else:
                    self.assertEqual(blocks[0]['status'],4)
                    type(self).conditioning[key]='nonfinite kinetic intermediates'
                self.assertEqual(blocks[1]['status'],0)
                for i,t in blocks[1]['chi'].items():
                    np.testing.assert_array_equal(t,baseline['chi'][i])

    def test_invalid_scale_and_species(self):
        species=[COLD_P,E]
        w=.13+.03j
        for checked in self.probes:
            blocks=self.evaluate(species,[(w,b,0) for b in (0.,-1.,float('nan'),float('inf'),1.)],checked=checked)
            self.assertEqual([b['status'] for b in blocks],[1,1,1,1,0])
            for va in (0.,-1.,float('nan')):
                self.assertEqual(self.evaluate(species,[(w,1.,0)],checked=checked,va=va)[0]['status'],1)
            for warm in (replace_species(E,tau=-1.),replace_species(E,alpha=0.),replace_species(E,mu=0.)):
                self.assertEqual(self.evaluate([COLD_P,warm],[(w,1.,0)],checked=checked)[0]['status'],1)

    def test_production_uses_core_not_extraction(self):
        # This is the actual input path; no experimental build is involved.
        species=[replace_species(COLD_P,alpha=1.),E]
        result,_=run(ROOT/'plume.e',input_text(species))
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        # The validated adapter/core are production code; extraction remains test-only.
        symbols=subprocess.run(['nm',str(ROOT/'plume.e')],check=True,capture_output=True,text=True).stdout
        self.assertIn('__auxiliary_normalization',symbols)
        self.assertIn('__kinetic_core',symbols)
        self.assertNotIn('__auxiliary_kinetic_kernel',symbols)


if __name__=='__main__':
    unittest.main()
