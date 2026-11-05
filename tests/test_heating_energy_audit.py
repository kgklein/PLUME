"""Read-only heating audit. Defect assertions characterize current behavior,
not requirements that future corrected implementations retain those defects.
"""
import hashlib
import json
import subprocess
import unittest
from pathlib import Path
import numpy as np
from heating_energy_audit import HeatingHarness,audit_case,P,E,CP,CE,A,CA,ROOT,scans,mixed
from heating_energy_reference import (mp,mp_warm,energy_reference,hermitian,antihermitian,contraction,absorption,
                                     channels,harmonic_reference,fd_derivative)
from cold_moment_reference import fluid_moments,vlasov_quadrature
from mixed_reference import species_tensor

BENCHMARKS={
 'kinetic':([P,E],{}),
 'warm_reference_mixed':([P,CE],{}),
 'cold_reference_mixed':([CP,E],{}),
 'fully_cold':([CP,CE],{}),
 'interleaved':([CP,E,A,CA],dict(guess=.043-.00002j)),
 'weak_growth':([(.8,1.,2.,1.,1.,0.),(.2,1836.,1.,-1.,1.,0.)],dict(kpar=.3,guess=.322+.0017j)),
 'reactive_mixed':([CP,mixed.replace_species(CE,density=1.1,drift=.2/1.1),(.03,.25,1.,.5,.05,2.)],
                   dict(kpar=.3,kperp=.1,guess=.3+.08j)),
 'reactive_cold':([CP,mixed.replace_species(CE,density=1.1,drift=.2/1.1),(0.,.25,0.,.5,.05,2.)],
                  dict(kpar=.3,kperp=.1,guess=.3+.08j)),
 'high_a_cold':([CP,(0.,4.,0.,-1.,1.,0.)],dict(kpar=.2,kperp=.1,a=.3,guess=8.+.01j)),
}


def error(actual,reference):return float(np.max(abs(np.asarray(actual)-reference)/np.maximum(1.,abs(np.asarray(reference)))))


class HeatingEnergyAuditTests(unittest.TestCase):
    measurements={};cases={}
    @classmethod
    def setUpClass(cls):
        cls.h=HeatingHarness();cls.addClassCleanup(cls.h.close)
    @classmethod
    def tearDownClass(cls):print('Heating audit measurements:',json.dumps(cls.measurements))

    def case(self,name,checked=False):
        key=f'{name}/{checked}'
        if key not in type(self).cases:
            species,kw=BENCHMARKS[name];type(self).cases[key]=audit_case(self.h,species,checked=checked,**kw)
        return type(self).cases[key]

    @unittest.skipIf(mp is None,'optional arbitrary precision: put mpmath on PYTHONPATH')
    def test_independent_high_precision_derivative(self):
        maxima=0.
        for s in [P,E,A]:
            c,d=mp_warm(s,.73,.31,.7,.02,n=24,dps=60)
            c2,d2=mp_warm(s,.73,.31,.7,.02,n=32,dps=90)
            maxima=max(maxima,error(c,c2),error(d,d2))
            np.testing.assert_allclose(d,d2,rtol=1e-13,atol=1e-10)
            np.testing.assert_allclose(c,species_tensor(s,.73,.31,.7,.02,128),rtol=1e-7,atol=1e-6)
        type(self).measurements['arbitrary_precision_stability']=maxima

    def test_legacy_species_heating_and_seven_channels(self):
        maxima=0.;closure=0.;shift=0.
        species=[P,E,A]
        for checked in [False,True]:
            for mode in ['adaptive','fixed128']:
                p=self.h.evaluate(species,[.06-.0001j],modes=1,checked=checked,mode=mode)[0]
                wr=p['w'].real;e=p['e'];atreal=self.h.evaluate(species,[complex(wr)],checked=checked,mode=mode)[0]
                derivative=hermitian((1.000001*p['after']-atreal['chi'])/1e-6)
                energy=sum(contraction(e,derivative))+np.vdot(p['b'],p['b']).real
                total=absorption(e,atreal['chi'])/energy
                splits=np.array([channels(e,c)/energy for c in p['after_low']])
                maxima=max(maxima,error(total,p['power']),error(splits,p['split']))
                np.testing.assert_allclose(p['power'],total,rtol=2e-9,atol=2e-13)
                np.testing.assert_allclose(p['split'],splits,rtol=2e-9,atol=2e-13)
                closure=max(closure,error(np.sum(splits[:,:4],axis=1),splits[:,4]))
                np.testing.assert_allclose(np.sum(splits[:,:4],axis=1),splits[:,4],rtol=5e-13,atol=5e-15)
                same=np.array([channels(e,c)/energy for c in atreal['low']])
                shift=max(shift,float(np.linalg.norm(splits-same)/max(1e-30,np.linalg.norm(same))))
        self.assertGreater(shift,1e-8)
        type(self).measurements.update(legacy_reconstruction_scaled=maxima,channel_n0_closure_scaled=closure,
                                      shifted_channel_relative_bias=shift)

    def test_signed_harmonics_and_channel_scope(self):
        e=np.array([1.,.3+.7j,.25-.1j]);maxima=0.;records=[]
        species=[P,E,A]
        for checked in [False,True]:
            for kz in [-.7,.7]:
                p=self.h.evaluate(species,[1.3+0j],checked=checked,kpar=kz,kperp=3.)[0]
                for s,chi,low in zip(species,p['chi'],p['low']):
                    terms,extra=harmonic_reference(s,1.3,kz,3.,.01,64)
                    selected=terms[[63,64,65]]
                    observed=antihermitian(low);expected=antihermitian(selected)
                    maxima=max(maxima,error(observed,expected))
                    np.testing.assert_allclose(observed,expected,rtol=2e-8,atol=1e-8)
                    np.testing.assert_allclose(absorption(e,terms).sum(),absorption(e,chi),rtol=2e-8,atol=1e-8)
                    c=channels(e,low);self.assertLess(abs(c[:4].sum()-c[4])/max(1,abs(c[4])),1e-12)
                    if s==P and not checked and kz>0:
                        q=float(absorption(e,chi));truncated=float(sum(c[4:]));all_low=float(sum(absorption(e,selected)))
                        records.append(dict(total=q,n0=float(c[4]),transverse_plus=float(c[5]),transverse_minus=float(c[6]),
                                            missing_fraction=(q-truncated)/q,
                                            higher_harmonic_fraction=(q-all_low)/q,
                                            excluded_pm1_longitudinal_fraction=(all_low-truncated)/q))
        self.assertGreater(abs(records[0]['missing_fraction']),1e-3)
        type(self).measurements.update(harmonic_absorptive_scaled=maxima,truncation=records)

    def test_fixed_scratch_initialization(self):
        outcomes=[]
        e=np.array([1.,.3+.7j,.25-.1j])
        for checked in [False,True]:
            for mode in ['fixed128','adaptive']:
                p=self.h.evaluate([P,E],[1.3+0j],modes=2,checked=checked,kpar=.7,kperp=3.,mode=mode)[0]
                unassigned=p['low'][:,1,:,:][:,np.array([0,0,0,1,2]),np.array([0,1,2,0,0])]
                if mode=='fixed128':np.testing.assert_array_equal(unassigned,2j*np.ones_like(unassigned))
                else:np.testing.assert_array_equal(unassigned,np.zeros_like(unassigned))
                c=channels(e,p['low'][0]);mismatch=float(abs(c[:4].sum()-c[4]))
                if mode=='fixed128':self.assertGreater(mismatch,1.)
                else:self.assertLess(mismatch,1e-10)
                outcomes.append(dict(checked=checked,mode=mode,n0_unassigned_entries=int(np.count_nonzero(unassigned)),channel_closure_error=mismatch))
        type(self).measurements['scratch_initialization']=outcomes

    def test_mixed_ingredients_availability_and_recovery(self):
        records=[];e=np.array([1,.2+.3j,.01-.02j])
        for checked in [False,True]:
            for species in [[P,CE],[CP,E],[CP,E,A,CA]]:
                p=self.h.evaluate(species,[.73+0j],modes=3,checked=checked,kpar=.31,kperp=.7)[0]
                for index,s in enumerate(species):
                    available=bool(np.isfinite(p['low'][index]).all())
                    if species[0][0]==0:self.assertFalse(available)
                    if s[0]>0:
                        self.assertTrue(np.isfinite(p['recover'][index]).all())
                        terms,_=harmonic_reference(s,.73,.31,.7,.01,64)
                        np.testing.assert_allclose(channels(e,p['recover'][index]),channels(e,terms[63:66]),rtol=3e-8,atol=1e-8)
                    else:self.assertTrue(np.isnan(p['recover'][index]).all())
                    records.append(dict(checked=checked,reference_cold=species[0][0]==0,species=index+1,kinetic=s[0]>0,production_low_available=available))
            parallel=self.h.evaluate([CP,E],[.06-.0001j],modes=1,checked=checked,kperp=0.)[0]
            self.assertEqual(parallel['status'],0);self.assertTrue(np.isnan(parallel['low']).all())
        type(self).measurements['mixed_availability']=records

    def test_stationary_cold_energy_and_vacuum(self):
        e=np.array([1.,.2+.3j,.1-.07j]);records=[]
        for a in [.01,.3]:
            species=[CP,CE];w=.73;kz=.31;kp=.7
            reference=energy_reference(species,w,kz,kp,a,e)
            b=np.cross([kp,0,kz],e)/(a*w)
            particle=sum(s[4]/(s[1]*a*a)*np.vdot(fluid_moments(s,w,kz,kp,a,e,b)[0],fluid_moments(s,w,kz,kp,a,e,b)[0]).real for s in species)
            direct=np.vdot(e,e).real+np.vdot(b,b).real+particle
            self.assertLess(abs(reference['energy']-direct)/direct,3e-12)
            # A binary64 fluid solve can leave roundoff-level antisymmetry;
            # the real-frequency cold absorption is zero within that bound.
            np.testing.assert_allclose(reference['q'],0.,rtol=0.,atol=16*np.finfo(float).eps*
                max(1.,np.linalg.norm(reference['chi']))*np.vdot(e,e).real)
            records.append(dict(a=a,dielectric_energy=reference['energy'],field_plus_cold_kinetic_energy=float(direct)))
        high=self.case('high_a_cold');self.assertGreater(high['vacuum_fraction'],.49);self.assertLess(high['vacuum_fraction'],.51)
        type(self).measurements['closed_cold_energy']=records

    def test_derivative_sweep_and_weak_energy_balance(self):
        results={}
        for name in ['kinetic','warm_reference_mixed','cold_reference_mixed','fully_cold','interleaved']:
            for checked in [False,True]:
                d=self.case(name,checked)
                self.assertLess(d['null_residual'],3e-8)
                self.assertLess(d['independent_root_error'],1e-9)
                if d['closure_relative'] is not None:self.assertLess(d['closure_relative'],2e-4)
                else:self.assertLess(abs(d['predicted_gamma']-d['w'].imag),2e-11)
                best=min(r['centered_matrix'] for r in d['derivative_errors'])
                self.assertLess(best,3e-7)
                self.assertLess(best,d['derivative_errors'][3]['forward_matrix'])
                results[f'{name}/{checked}']=dict(w=[d['w'].real,d['w'].imag],energy=d['energy_reference'],
                    closure_relative=d['closure_relative'],best_centered_matrix=best,
                    historical_forward_matrix=d['derivative_errors'][3]['forward_matrix'])
        type(self).measurements['weak_balance']=results

    def test_instability_interpretation(self):
        weak=self.case('weak_growth')
        self.assertGreater(weak['w'].imag,0.);self.assertLess(weak['weak_ratio'],.01)
        self.assertLess(weak['closure_relative'],.01)
        results={}
        for name in ['weak_growth','reactive_mixed','reactive_cold']:
            d=self.case(name)
            if name!='weak_growth':
                self.assertGreater(d['w'].imag,.07);self.assertLess(abs(d['predicted_gamma']),1e-5)
            results[name]=dict(w=[d['w'].real,d['w'].imag],energy=d['energy_reference'],
                              real_frequency_q=d['raw'].tolist(),predicted_gamma=d['predicted_gamma'],closure_relative=d['closure_relative'])
        type(self).measurements['instability']=results

    def test_rho_d_equivalence_and_historical_backend(self):
        maxima=0.
        for checked in [False,True]:
            for backend in ['legacy','faddeeva']:
                for mode in ['fixed128','adaptive']:
                    d=self.h.evaluate([P,E],[.06-.0001j],modes=1,checked=checked,mode=mode,backend=backend)[0]
                    r=self.h.evaluate([P,E],[d['w']],modes=0,checked=checked,mode=mode,backend=backend,normalization='rho_ref')[0]
                    maxima=max(maxima,error(d['power'],r['power']),error(d['split'],r['split']))
                    np.testing.assert_allclose(d['power'],r['power'],rtol=2e-7,atol=2e-12)
                    np.testing.assert_allclose(d['split'],r['split'],rtol=2e-7,atol=2e-12)
        type(self).measurements['normalization_scaled']=maxima

    def test_actual_output_and_python_column_mapping(self):
        import sys
        sys.path.insert(0,str(ROOT))
        from linfpclib.linfpc import load_plume_sweep
        records=[];species=[P,E,A];ns=3
        for option in [1,2]:
            for eigen in [False,True]:
                axes=[scans.scan(1,.06,.07,heat=True,eigen=eigen)] if option==1 else [scans.scan(0,.05,.06,steps=1),scans.scan(1,.06,.07,steps=1,heat=True,eigen=eigen)]
                result,files=mixed.run(ROOT/'plume.e',scans.source(species,axes,option=option))
                self.assertEqual(result.returncode,0,result.stdout+result.stderr)
                rows=next(iter(scans.mode_rows(files).values()))
                front=18+8*ns if eigen else 6
                self.assertEqual(rows.shape[1],(18+22*ns if eigen else 6+14*ns))
                marker=rows[0].copy();marker[front:front+8*ns]=np.arange(8*ns)+1000
                path=self.h.directory/'marker.mode1';path.write_text(' '.join(str(x) for x in marker)+'\n')
                loaded=load_plume_sweep(str(path),nspec=ns,heating=True,eigen=eigen)
                self.assertEqual(loaded['p1'][0],1000)
                expected=1010.;observed=float(loaded['p2ttd_yy'][0])
                self.assertEqual(observed,expected) # Reader correction; historical audit JSON is retained.
                records.append(dict(option=option,eigen=eigen,columns=rows.shape[1],true_p1=1000,read_p1=float(loaded['p1'][0]),
                                    true_species2_yy=expected,read_species2_yy=observed))
        type(self).measurements['reader_columns']=records

    def test_ew_unassigned_and_format_reversion(self):
        exe=self.h.canary_executable();records=[]
        for eigen in [False,True]:
            result,files=mixed.run(exe,scans.source([P,E],[scans.scan(1,.06,.07,heat=True,eigen=eigen)],extra='low_n=.false.'))
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            text=next(t for name,t in files.items() if name.endswith('.mode1'))
            rows=[l.split() for l in text.splitlines() if l.strip()]
            if eigen:
                self.assertTrue(all(np.isnan(float(r[-1])) for r in rows))
                self.assertTrue(all(len(r)==49 for r in rows))
            else:
                self.assertEqual([len(r) for r in rows],[20,1,20,1,20,1])
                self.assertTrue(all(np.isnan(float(r[0])) for r in rows[1::2]))
            records.append(dict(eigen=eigen,row_lengths=[len(r) for r in rows],energy_is_initialized_nan=True))
        type(self).measurements['unassigned_ew']=records

    def test_collisional_absorption_velocity_integral(self):
        # Real laboratory w, Krook shift only in the Vlasov operator. Faraday B
        # uses the real laboratory frequency, not the collisional shifted one.
        e=np.array([1.,.2+.3j,.05-.02j]);w=.73;kz=.1;kp=.3;a=.02;records=[]
        species=[P,(1e-4,1836.,.8,-1.,1.,-.2)]
        for nu in [.01,.1]:
            kn=np.sqrt(P[0]*P[1])/(np.sqrt(2)*nu)
            p=self.h.evaluate(species,[complex(w)],checked=False,kpar=kz,kperp=kp,va=a,extra=f'collision_type=1\nKn={kn:.17g}')[0]
            b=np.cross([kp,0,kz],e)/(a*w)
            for s,chi in zip(species,p['chi']):
                rate=np.sqrt(s[0]*s[1])/(np.sqrt(2)*kn)
                v,n=vlasov_quadrature(s,w+1j*rate,kz,kp,a,e,b,32,32)
                v2,n2=vlasov_quadrature(s,w+1j*rate,kz,kp,a,e,b,48,40)
                current=s[4]/s[3]*(v+np.array([0,0,a*s[5]*n]))
                current2=s[4]/s[3]*(v2+np.array([0,0,a*s[5]*n2]))
                rhs=-1j*w*a*a*(chi@e)
                self.assertLess(error(current,current2),2e-7)
                # The corrected longitudinal row agrees. The oblique transverse
                # response has a separately characterized finite-rate symmetry
                # defect; do not assert full agreement or repair it in an audit.
                self.assertLess(error(current2[2],rhs[2]),2e-7)
                self.assertGreater(error(current2,rhs),1e-7)
                direct=float(np.vdot(e,current2).real/(w*a*a));q=float(absorption(e,chi))
                records.append(dict(collision=rate,mass=s[1],current_scaled=error(current2,rhs),
                                    q_vlasov=direct,q_susceptibility=q,relative_absorption=abs(direct-q)/max(1,abs(q))))
        type(self).measurements['collisional_velocity_integral']=records

    def test_collisional_oblique_symmetry_isolation(self):
        w=.73;kz=.1;kp=.3;a=.02;nu=.1;records=[]
        kn=np.sqrt(P[0])/(np.sqrt(2)*nu)
        for anisotropic,species in [(True,[P,(1e-4,1836.,.8,-1.,1.,-.2)]),
                                    (False,[(.4,1.,1.,1.,1.,0.),(1e-4,1836.,1.,-1.,1.,0.)])]:
            for mode in ['adaptive','fixed128']:
                p=self.h.evaluate(species,[complex(w)],kpar=kz,kperp=kp,va=a,mode=mode,
                                  extra=f'collision_type=1\nKn={kn:.17g}')[0]
                for s,chi in zip(species,p['chi']):
                    rate=np.sqrt(s[0]*s[1])/(np.sqrt(2)*kn);quad=np.zeros((3,3),complex);coarse=quad.copy()
                    for col in range(3):
                        e=np.eye(3,dtype=complex)[col];b=np.cross([kp,0,kz],e)/(a*w)
                        for order,angles,target in [(48,40,coarse),(64,48,quad)]:
                            v,n=vlasov_quadrature(s,w+1j*rate,kz,kp,a,e,b,order,angles)
                            current=s[4]/s[3]*(v+np.array([0,0,a*s[5]*n]));target[:,col]=1j*current/(w*a*a)
                    convergence=error(quad,coarse);self.assertLess(convergence,1e-7)
                    mask=np.ones((3,3),bool);mask[:2,2]=False
                    self.assertLess(error(quad[mask],chi[mask]),1e-7)
                    discrepancy=error(quad[:2,2],chi[:2,2])
                    if anisotropic:self.assertGreater(discrepancy,.01)
                    else:self.assertLess(discrepancy,1e-7)
                    records.append(dict(mode=mode,anisotropic_or_drifting=anisotropic,mass=s[1],nu=rate,
                        convergence_scaled=convergence,matched_components_scaled=error(quad[mask],chi[mask]),
                        xz_yz_scaled=discrepancy,production=np.stack([chi.real,chi.imag],axis=-1).tolist(),
                        vlasov=np.stack([quad.real,quad.imag],axis=-1).tolist()))
        type(self).measurements['collisional_oblique_symmetry']=records

    def test_repeated_frequency_maps_have_no_heating_product(self):
        axes=[scans.scan(1,.06,.08,steps=1,res=1,heat=True)]
        text=scans.source([P,E],axes,option=4)
        text+=('&maps\nloggridw=.false.\nloggridg=.false.\nomi=.035\nomf=.11\n'
               'gami=-.001\ngamf=.001\nnr=12\nni=6\npositive_roots=.true.\n/\n')
        result,files=mixed.run(ROOT/'plume.e',text)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        self.assertFalse(any('.mode' in n for n in files))
        maps=[t for n,t in files.items() if n.endswith('.map')]
        roots=[t for n,t in files.items() if n.endswith('.roots')]
        self.assertEqual(len(maps),2);self.assertEqual(len(roots),2)
        self.assertTrue(all(len(l.split())==5 for t in maps for l in t.splitlines() if l.strip()))
        self.assertTrue(all(len(l.split())==18 for t in roots for l in t.splitlines() if l.strip()))
        rejected=scans.source([CP,E],axes,option=4)+text[text.index('&maps'):]
        failure,_=mixed.run(ROOT/'plume.e',rejected)
        self.assertNotEqual(failure.returncode,0);self.assertIn('heating is unavailable',failure.stdout+failure.stderr)
        type(self).measurements['frequency_maps']=dict(map_files=2,map_columns=5,root_files=2,root_columns=18,
                                                     kinetic_heating_flag_unused=True,mixed_heating_rejected=True)
