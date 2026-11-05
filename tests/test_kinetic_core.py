from faddeeva_build import compile_c,c_objects
"""Milestone A checkpoint: explicit production core versus unchanged prototype."""
import subprocess
import unittest
from pathlib import Path
import numpy as np
import test_auxiliary_normalization as aux


class KineticCoreTests(unittest.TestCase):
    maximum_error=0.

    @classmethod
    def setUpClass(cls):
        aux.AuxiliaryNormalizationTests.setUpClass()
        cls.addClassCleanup(aux.AuxiliaryNormalizationTests.temp.cleanup)
        cls.probes={}
        modules=('nrtype','nrutil_trim','vars','plasma_z_faddeeva','plasma_z_legacy','plasma_z_backend','cold_susceptibility','auxiliary_normalization','kinetic_core','harmonic_selection','bessel_cache','adaptive_kinetic','cold_reference','species_moments','electromagnetic_fields','heating_response','heating_energy','heating_diagnostics','scan_support','functions','bessel','disprels')
        for checked,prototype in aux.AuxiliaryNormalizationTests.probes.items():
            directory=prototype.parent
            base=['gfortran','-fdefault-real-8','-I',str(directory),'-J',str(directory)]
            flags=['-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow'] if checked else ['-O3']
            if not checked:
                base+=['-I',str(aux.ROOT/'include')]
            subprocess.run(base+flags+['-c',str(aux.ROOT/'src/kinetic_core.f90')],cwd=directory,
                           check=True,capture_output=True)
            objects=[str((directory if checked else aux.ROOT)/f'{m}.o') for m in modules]
            exe=directory/'core_probe'
            subprocess.run(base+flags+[str(directory/'kernel.f90'),str(aux.ROOT/'tests/kinetic_core_probe.f90'),
                           *objects,*c_objects(directory if checked else aux.ROOT),'-o',str(exe)],cwd=directory,check=True,capture_output=True)
            cls.probes[checked]=exe

    run_probe=aux.AuxiliaryNormalizationTests.run_probe
    evaluate=aux.AuxiliaryNormalizationTests.evaluate

    @classmethod
    def tearDownClass(cls):
        print(f'Kinetic core maximum scaled prototype/invariance error={cls.maximum_error:.3e}')

    def compare(self,a,b):
        type(self).maximum_error=max(self.maximum_error,float(np.max(abs(a-b)/np.maximum(1.,abs(b)))))
        np.testing.assert_allclose(a,b,rtol=5e-9,atol=1e-10)

    def test_core_matches_prototype_and_scale_invariance(self):
        systems=[[aux.P,aux.E],
                 [aux.replace_species(aux.P,alpha=1.7),aux.replace_species(aux.E,alpha=.8,drift=-.1),aux.ALPHA],
                 [aux.COLD_P,aux.replace_species(aux.E,alpha=.65,drift=.15),aux.ALPHA]]
        for checked in self.probes:
            for species in systems:
                for kpar in (.06,-.06):
                    for va in (.001,.05):
                        for w in (.13+.03j,.13-.03j):
                            requests=[(w,b,mode) for b in aux.SCALES for mode in (0,3)]
                            blocks=self.evaluate(species,requests,checked=checked,kpar=kpar,va=va)
                            self.assertTrue(all(b['status']==0 for b in blocks))
                            for index in range(0,len(blocks),2):
                                for i in blocks[index]['chi']:
                                    self.compare(blocks[index+1]['chi'][i],blocks[index]['chi'][i])
                                    self.compare(blocks[index+1]['chi'][i],blocks[1]['chi'][i])

    def test_cold_placeholder_and_state_preservation(self):
        for checked in self.probes:
            baseline=None
            for alpha in (float('nan'),0.,-10.,1e100):
                species=[aux.replace_species(aux.COLD_P,alpha=alpha),aux.E,aux.ALPHA]
                block=self.evaluate(species,[(.13+.03j,.2,3)],checked=checked)[0]
                self.assertEqual(block['status'],0)
                if baseline:
                    for i in baseline:
                        np.testing.assert_array_equal(block['chi'][i],baseline[i])
                baseline=block['chi']

    def test_arithmetic_failure_and_recovery(self):
        for checked in self.probes:
            requests=[(.13+.03j,b,3) for b in (1.,1e160,1.)]
            blocks=self.evaluate([aux.COLD_P,aux.E,aux.ALPHA],requests,checked=checked)
            self.assertEqual([b['status'] for b in blocks],[0,5,0])
            for i in blocks[0]['chi']:
                np.testing.assert_array_equal(blocks[0]['chi'][i],blocks[2]['chi'][i])

    def test_invalid_core_wavevectors(self):
        for checked in self.probes:
            for kpar,kperp in ((0.,.05),(.06,0.)):
                block=self.evaluate([aux.P,aux.E],[(.13+.03j,1.,3)],checked=checked,kpar=kpar,kperp=kperp)[0]
                self.assertEqual(block['status'],1)
