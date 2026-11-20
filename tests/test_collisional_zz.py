from faddeeva_build import c_objects
"""Focused correction tests using moment and direct velocity-integral oracles."""
import json
import subprocess
import unittest
import numpy as np
from plasma_dispersion_audit import AuditHarness,ROOT,MODULES,metrics
from collisional_zz_reference import coefficients,zz_moments,zz_quadrature
import test_mixed_dispersion as mixed

S=(2.,.25,1.3,.5,.05,.4)


class CollisionalZZTests(unittest.TestCase):
    maxima={'moments':0.,'quadrature':0.,'harmonics':0.,'zero_rate':0.,'production':0.,'unchanged':0.,
            'heating_relative':0.}

    @classmethod
    def setUpClass(cls):
        subprocess.run(['make','-j2'],cwd=ROOT,check=True,capture_output=True)
        cls.harnesses={False:AuditHarness(),True:AuditHarness(checked=True)}
        cls.probes={};cls.diagnostic_probes={}
        for checked,h in cls.harnesses.items():
            cls.addClassCleanup(h.close)
            exe=h.directory/'production_zz'
            flags=['-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow'] if checked else ['-O3']
            objects=h.directory if checked else ROOT
            result=subprocess.run(['gfortran','-fdefault-real-8',*flags,'-I',str(h.directory if checked else ROOT/'include'),
                 str(ROOT/'tests/collisional_zz_probe.f90'),*[str(objects/f'{m}.o') for m in MODULES],*c_objects(objects),
                 '-o',str(exe)],capture_output=True,text=True)
            if result.returncode:raise RuntimeError(result.stderr)
            cls.probes[checked]=exe
            diagnostic=h.directory/'normalization_probe'
            result=subprocess.run(['gfortran','-fdefault-real-8',*flags,'-I',str(h.directory if checked else ROOT/'include'),
                 str(ROOT/'tests/normalization_probe.f90'),*[str(objects/f'{m}.o') for m in MODULES],*c_objects(objects),
                 '-o',str(diagnostic)],capture_output=True,text=True)
            if result.returncode:raise RuntimeError(result.stderr)
            cls.diagnostic_probes[checked]=diagnostic

    @classmethod
    def tearDownClass(cls):print('Collisional zz maximum scaled errors:',cls.maxima)

    def compare(self,actual,expected,category,tol=2e-10,atol=2e-11):
        error=metrics(actual,expected)['scaled']
        type(self).maxima[category]=max(type(self).maxima[category],error)
        np.testing.assert_allclose(actual,expected,rtol=tol,atol=atol)

    def test_velocity_integral_and_independent_moments(self):
        for s in [S,(.4,1.,1.,1.,1.,0.),(.2,1836.,.8,-1.,1.1,-.2),(.03,.25,2.,.5,.04,-.3)]:
            for sign in [-1,1]:
                k=sign*.31;w=.73+max(.5,abs(k)*np.sqrt(s[0]*s[1]))*1j
                for nu in [0.,.01,.1]:
                    with self.subTest(s=s,k=k,nu=nu):
                        q1,total1=zz_quadrature(s,w,k,3.,.02,nu,64,points=128)
                        q2,total2=zz_quadrature(s,w,k,3.,.02,nu,64,points=256)
                        self.compare(total1,total2,'quadrature')
                        terms,_,moment=zz_moments(s,w,k,3.,.02,nu,64)
                        self.compare(moment,total2,'moments')
                        for checked,h in self.harnesses.items():
                            d=h.tensor(s,w,kpar=k,kperp=3.,n=64,collision=nu,collisional=True,method=2)
                            self.assertFalse(d['nonfinite'])
                            self.compare(d['chi'][2,2],total2,'quadrature')
                            for index in [-1,0,1]:
                                self.compare(d['low'][index+1,2,2],q2[64+index],'harmonics')

    def test_fixed_special_functions_and_landau_continuation(self):
        for checked,h in self.harnesses.items():
            for s in [S,(.4,1.,.7,-1.,.7,-.15),(.2,25.,1.8,-1.,1.,.03)]:
                for k in [-.31,.31]:
                    for w in [.73+.13j,.73-.08j]:
                        for nu in [0.,1e-4,.01,.1]:
                            d=h.tensor(s,w,kpar=k,kperp=3.,n=64,collision=nu,collisional=True)
                            terms,_,total=zz_moments(s,w,k,3.,.02,nu,64,zvalues=d['z'])
                            self.compare(d['chi'][2,2],total,'moments')
                            self.assertTrue(np.isfinite(d['low']).all())
                            for index in [-1,0,1]:
                                self.compare(d['low'][index+1,2,2],terms[64+index],'harmonics')

    def test_zero_rate_limit_and_decreasing_rates(self):
        for h in self.harnesses.values():
            for s in [S,(2.,.25,1.,.5,.05,.4),(.4,1.,.7,-1.,.7,-.15)]:
                for k in [-.31,.31]:
                    base=h.tensor(s,.73+.5j,kpar=k,kperp=3.,n=64,method=2)['chi']
                    zero=h.tensor(s,.73+.5j,kpar=k,kperp=3.,n=64,method=2,collisional=True)['chi']
                    self.compare(zero,base,'zero_rate',tol=2e-11)
                    errors=[]
                    for nu in [1e-2,1e-3,1e-4,1e-5]:
                        d=h.tensor(s,.73+.5j,kpar=k,kperp=3.,n=64,method=2,collision=nu,collisional=True)
                        errors.append(abs(d['chi'][2,2]-base[2,2]))
                    self.assertTrue(all(b<a*.2 for a,b in zip(errors,errors[1:])),errors)

    def test_diagnostic_reconstruction_and_drift_allocation(self):
        for h in self.harnesses.values():
            for nu in [0.,.1]:
                gamma=coefficients(S,3.,1)
                d=h.tensor(S,.73+.5j,kpar=-.31,kperp=3.,n=1,collision=nu,collisional=True,method=2)
                terms,remainder,total=zz_moments(S,.73+.5j,-.31,3.,.02,nu,1,gamma=gamma,zvalues=d['z'])
                self.assertGreater(abs(remainder),1.)
                self.compare(d['low'][:,2,2],terms,'harmonics')
                self.compare(d['low'][:,2,2].sum()+remainder,d['chi'][2,2],'harmonics')
                self.compare(d['chi'][2,2],total,'moments')

    def test_real_frequency_dissipative_limit_and_repeat_slots(self):
        for h in self.harnesses.values():
            baseline=h.tensor(S,.73+0j,kpar=-.31,kperp=3.,n=64,method=2)
            first=h.tensor(S,.73+0j,kpar=-.31,kperp=3.,n=64,method=2,collisional=True)
            # Nonresonant allocation differs from legacy collisionless diagnostics;
            # its real-frequency zz anti-Hermitian part must have the same limit.
            self.compare(first['low'][:,2,2].imag,baseline['low'][:,2,2].imag,'harmonics')
            h.tensor(S,.81-.03j,kpar=.31,kperp=.5,n=3,collision=.1,method=2,collisional=True)
            repeat=h.tensor(S,.73+0j,kpar=-.31,kperp=3.,n=64,method=2,collisional=True)
            np.testing.assert_array_equal(first['low'],repeat['low'])
            self.assertNotEqual(first['low'][0,2,2],first['low'][2,2,2])

    def test_other_tensor_elements_preserve_frozen_audit_values(self):
        frozen=json.loads((ROOT/'tests/PLASMA_DISPERSION_AUDIT_RESULTS.json').read_text())
        expected=np.asarray(frozen['production_zero_collision_confirmation']['production_collisional_chi'])
        expected=expected[...,0]+1j*expected[...,1]
        current=self.harnesses[False].tensor(S,.73+.13j,kpar=-.31,kperp=3.,beta=.4,alpha_scale=1.7,
                      n=1024,collisional=True)['chi']
        mask=np.ones((3,3),bool);mask[2,2]=False
        self.compare(current[mask],expected[mask],'unchanged',tol=2e-13)

    def production(self,mode,normalization,checked,nu):
        species=[(.4,1.,1.7,1.,1.,0.),S]
        source=mixed.input_text(species,kpar=-.31,kperp=3.,va=.02,
                    extra=f"collision_type=1\nKn=1e30\nbessel_mode='{mode}'\nbessel_nfixed=128")
        if normalization=='rho_ref':
            import re
            length=np.sqrt(.4*1.7)
            source=source.replace("normalization='d_ref'","normalization='rho_ref'")
            for key,value in [('kpar',-.31*length),('kperp',3.*length),('vtp',.02*np.sqrt(.4))]:
                source=re.sub(rf'^{key}=.*$',f'{key}={value:.17g}',source,flags=re.M)
            source=source.replace('tauS=0.40000000000000002','tauS=1').replace('tauS=2\n','tauS=0.2\n')
        result,_=mixed.run(self.probes[checked],source,f'.73 .5 {nu}\n.81 -.03 {nu}\n.73 .5 {nu}\n')
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        blocks=[]
        for line in result.stdout.splitlines():
            fields=line.split()
            if not fields:continue
            if fields[0]=='DET':
                self.assertEqual(fields[1],'0')
                blocks.append({'chi':np.empty((2,3,3),complex),'low':np.empty((2,3,3,3),complex),
                               'gamma':[[],[]],'n':[0,0]})
            elif fields[0]=='CHI':
                _,s,i,j,vr,vi=fields;blocks[-1]['chi'][int(s)-1,int(i)-1,int(j)-1]=complex(float(vr),float(vi))
            elif fields[0]=='LOW':
                _,s,i,j,n,vr,vi=fields;blocks[-1]['low'][int(s)-1,int(n)+1,int(i)-1,int(j)-1]=complex(float(vr),float(vi))
            elif fields[0]=='HARMONICS':blocks[-1]['n'][int(fields[1])-1]=int(fields[2])
            elif fields[0]=='GAMMA':blocks[-1]['gamma'][int(fields[1])-1].append(float(fields[3]))
        self.assertEqual(len(blocks),3,result.stdout)
        np.testing.assert_array_equal(blocks[0]['low'],blocks[2]['low'])
        for b in blocks:self.assertTrue(np.isfinite(b['low']).all())
        return species,blocks[0]

    def test_zero_rate_root_and_supplementary_diagnostics(self):
        import test_normalization as norm
        species=[(1.,1.,1.7,1.,1.,0.),(.5,1836.,.8,-1.,1.,.02)]
        for exe in self.diagnostic_probes.values():
            for mode in ['fixed','adaptive']:
                source=norm.input_text('d_ref',.4,1.7,.05,.06,species)
                source=source.replace("bessel_mode='fixed'",f"bessel_mode='{mode}'")
                records=[]
                for collision in [0,1]:
                    text=source.replace('&params',f'&params\ncollision_type={collision}\nKn=1e30',1)
                    out,_=norm.run_case(exe,text)
                    self.assertEqual(out.returncode,0,out.stdout+out.stderr)
                    data=norm.numeric_records(out.stdout)
                    self.assertGreater(data[('ROOT',)][2],0)
                    self.assertLess(data[('ROOT',)][2],75)
                    # Equivalent equations may take a different last secant
                    # iteration; compare the converged frequency and diagnostics.
                    data[('ROOT',)]=data[('ROOT',)][:2]
                    records.append(data)
                # Wave-energy normalization uses a 1e-6 finite difference.
                # Last-bit changes between equivalent zz polynomials are
                # amplified there. Keep the original comparison for all other
                # records and explicitly bound this separate conditioning.
                stable=[{k:v for k,v in data.items() if k[0] not in ('HEAT','SPLIT')} for data in records]
                norm.NormalizationTests().compare_records(*stable)
                for key in records[0]:
                    if key[0] not in ('HEAT','SPLIT'):continue
                    a=np.array(records[0][key]);b=np.array(records[1][key])
                    relative=np.max(abs(a-b)/np.maximum(1e-8,abs(a)))
                    type(self).maxima['heating_relative']=max(type(self).maxima['heating_relative'],float(relative))
                    np.testing.assert_allclose(a,b,rtol=1e-5,atol=1e-10)

    def test_production_fixed_adaptive_and_normalizations(self):
        for checked,h in self.harnesses.items():
            for mode in ['fixed','adaptive']:
                for nu in [0.,.01,.1]:
                    normalized=[]
                    for convention in ['rho_ref','d_ref']:
                        species,data=self.production(mode,convention,checked,nu)
                        normalized.append(data['chi'])
                        for index,s in enumerate(species):
                            n=data['n'][index];gamma=np.array(data['gamma'][index])
                            # Actual production Gamma values; the Z inputs are
                            # evaluated by the unchanged production implementation.
                            diagnostic=h.tensor(s,.73+.5j,kpar=-.31,kperp=3.,va=.02,beta=.4,alpha_scale=1.7,
                                       n=n,gamma=gamma,collision=nu,collisional=True)
                            terms,_,total=zz_moments(s,.73+.5j,-.31,3.,.02,nu,n,gamma=gamma,zvalues=diagnostic['z'])
                            self.compare(data['chi'][index,2,2],total,'production')
                            self.compare(data['low'][index,:,2,2],terms[n-1:n+2],'production')
                    self.compare(normalized[0],normalized[1],'production')
