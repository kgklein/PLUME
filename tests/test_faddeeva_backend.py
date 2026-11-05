"""Production backend dispatch, physics, failure propagation and compatibility."""
import json
import re
import subprocess
import unittest
from pathlib import Path
import numpy as np
from scipy.optimize import root
from scipy.special import wofz
from faddeeva_build import c_objects
from plasma_dispersion_audit import AuditHarness,MODULES,ROOT,stress_systems,metrics,stix_terms,scipy_z
from mixed_reference import species_tensor,wave_matrix
from collisional_zz_reference import zz_moments
import test_mixed_dispersion as mixed


class FaddeevaBackendTests(unittest.TestCase):
    measurements={'shared_bessel_error':0.,'stress':{},'roots':[],'max_tensor_error':0.,'max_root_difference':0.}
    @classmethod
    def setUpClass(cls):
        subprocess.run(['make','-j2'],cwd=ROOT,check=True,capture_output=True)
        cls.h=AuditHarness(checked=True);cls.addClassCleanup(cls.h.close)
        cls.probes={};cls.direct={}
        for checked in [False,True]:
            directory=cls.h.directory if checked else ROOT
            flags=['-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow'] if checked else ['-O3']
            inc=directory if checked else ROOT/'include'
            objects=[str(directory/f'{m}.o') for m in MODULES]+c_objects(directory)
            count=cls.h.directory/f'count{checked}.o'
            subprocess.run(['cc','-std=c99','-O3','-c',str(ROOT/'tests/faddeeva_count_probe.c'),'-o',str(count)],check=True)
            base=['gfortran','-fdefault-real-8',*flags,'-I',str(inc)]
            for name,source in [('production','faddeeva_production_probe'),('direct','faddeeva_selected_probe')]:
                exe=cls.h.directory/f'{name}{checked}'
                command=base+[str(ROOT/'tests'/f'{source}.f90'),*objects]
                if name=='production':command += [str(count),'-Wl,--wrap=plume_faddeeva_z_value']
                result=subprocess.run(command+['-o',str(exe)],capture_output=True,text=True)
                if result.returncode:raise RuntimeError(result.stderr)
                (cls.probes if name=='production' else cls.direct)[checked]=exe

    @classmethod
    def tearDownClass(cls):
        print('Faddeeva production measurements:',json.dumps(cls.measurements))

    def evaluate(self,species,frequencies,backend='faddeeva',mode='adaptive',checked=False,roots=False,
                 repeats=1,normalization='d_ref',fixed=128,**kw):
        source=mixed.input_text(species,**kw).replace("z_backend='legacy'",f"z_backend='{backend}'")
        if backend is None:source=source.replace("z_backend='None'",'')
        source=source.replace("bessel_mode='fixed'",f"bessel_mode='{mode}'\nbessel_nfixed={fixed}\nbessel_nmax=1024")
        if normalization=='rho_ref':
            beta,alpha=species[0][0],species[0][2];length=np.sqrt(beta*alpha)
            source=source.replace("normalization='d_ref'","normalization='rho_ref'")
            for key,v in [('kpar',kw.get('kpar',.06)*length),('kperp',kw.get('kperp',.05)*length),
                          ('vtp',kw.get('va',.001)*np.sqrt(beta))]:
                source=re.sub(rf'^{key}=.*$',f'{key}={v:.17g}',source,flags=re.M)
            chunks=source.split('&species_')
            for i,s in enumerate(species):chunks[i+1]=chunks[i+1].replace(f'tauS={s[0]:.17g}',f'tauS={beta/s[0]:.17g}')
            source='&species_'.join(chunks)
        stdin=''.join(f'{w.real:.17g} {w.imag:.17g} {int(roots)} {repeats}\n' for w in frequencies)
        result,files=mixed.run(self.probes[checked],source,stdin)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        blocks,found=mixed.parse(result.stdout,len(species))
        index=-1
        for line in result.stdout.splitlines():
            fields=line.split()
            if fields and fields[0]=='DET':
                index+=1;blocks[index]['gamma']={}
            elif fields and fields[0]=='GAMMA':
                blocks[index]['gamma'].setdefault(int(fields[1])-1,[]).append(float(fields[3]))
        return blocks,found,result.stdout,files

    @staticmethod
    def matched_bessel_reference(s,w,args,gamma):
        n=len(gamma)-2
        c,b,extra=stix_terms(s,w,args['kpar'],args['kperp'],args['va'],n,gamma=gamma)
        orders=np.arange(-n,n+1)
        z=scipy_z(args['kpar'],(w-args['kpar']*s[5]-orders*s[1]/s[3])/(args['kpar']*np.sqrt(s[0]*s[1])))
        return np.asarray(np.sum(np.asarray(c+b*z,dtype=np.clongdouble),axis=2)+extra,dtype=complex)

    def test_selected_interface_signed_transition_and_range(self):
        args=[(k,z) for k in [-.31,.31] for z in [0j,3.01+1e-12j,3.01-1e-12j,3.01+.1j,
               2+.01j,2-.01j,1.+1.j,-14j,-26j,-30j,complex(np.nan,0),complex(np.inf,0)]]
        args += [(0.,1+.1j),(float('inf'),1+.1j)]
        for checked,exe in self.direct.items():
            data=''.join(f'faddeeva {k} {z.real} {z.imag}\n' for k,z in args)
            result=subprocess.run([exe],input=data,text=True,capture_output=True,check=True)
            actual=[complex(*map(float,line.split())) for line in result.stdout.splitlines()]
            for (k,z),a in zip(args,actual):
                if k==0 or not np.isfinite(k) or not np.isfinite(z):
                    self.assertTrue(np.isnan(a.real) and np.isnan(a.imag));continue
                with np.errstate(all='ignore'):ref=np.sign(k)*1j*np.sqrt(np.pi)*wofz(np.sign(k)*z)
                if np.isfinite(ref):np.testing.assert_allclose(a,ref,rtol=1e-13,atol=1e-13)
                else:self.assertTrue(np.isnan(a.real) and np.isnan(a.imag))

    def test_stress_point_both_harmonic_modes(self):
        species=stress_systems()[1][1];w=.73+.13j
        args=dict(kpar=-.31,kperp=3.,va=.02)
        expected=np.array([species_tensor(s,w,**args,harmonics=1024) for s in species])
        for checked in [False,True]:
            for mode,fixed in [('fixed',1024),('adaptive',128)]:
                accurate=self.evaluate(species,[w],mode=mode,fixed=fixed,checked=checked,**args)[0][0]
                old=self.evaluate(species,[w],backend='legacy',mode=mode,fixed=fixed,checked=checked,**args)[0][0]
                self.assertEqual(accurate['status'],0)
                error=metrics(accurate['chi'],expected)['scaled'];legacy=metrics(old['chi'],expected)['scaled']
                self.assertLess(error,3e-10);self.assertGreater(legacy,7e-4)
                matched=self.matched_bessel_reference(species[2],w,args,accurate['gamma'][2])
                if mode=='fixed':self.assertLess(metrics(accurate['chi'][2],expected[2])['scaled'],2e-13)
                self.assertLess(metrics(accurate['chi'][2],matched)['scaled'],2e-13)
                type(self).measurements['stress'][f'{checked}/{mode}']={'accurate':error,'legacy':legacy,
                    'per_species':[metrics(accurate['chi'][i],expected[i])['scaled'] for i in range(len(species))],
                    'alpha_matched_gamma':metrics(accurate['chi'][2],matched)['scaled']}

    def test_all_paths_and_independent_bessel_choices(self):
        for name,species in stress_systems():
            for mode in ['adaptive','fixed128']:
                for sign in [-1,1]:
                    for w in [.73+.13j,.73-.08j]:
                        args=dict(kpar=sign*.31,kperp=3.,va=.02)
                        out=self.evaluate(species,[w],mode=mode,**args)[0][0]
                        self.assertEqual(out['status'],0)
                        expected=np.array([species_tensor(s,w,**args,harmonics=1024) for s in species])
                        error=metrics(out['chi'],expected)['scaled']
                        type(self).measurements['max_tensor_error']=max(type(self).measurements['max_tensor_error'],error)
                        # Raw special-function accuracy is recorded separately;
                        # NR fixed128 and the vector seed approximations are unchanged.
                        self.assertTrue(np.isfinite(out['chi']).all())
                        if error==type(self).measurements['max_tensor_error']:
                            type(self).measurements['raw_bessel_floor_case']=[name,mode,sign,[w.real,w.imag]]
                        matched=np.array([self.matched_bessel_reference(s,w,args,out['gamma'][i]) if s[0]>0
                                          else species_tensor(s,w,**args) for i,s in enumerate(species)])
                        matched_error=metrics(out['chi'],matched)['scaled']
                        type(self).measurements['shared_bessel_error']=max(type(self).measurements['shared_bessel_error'],matched_error)
                        np.testing.assert_allclose(out['chi'],matched,rtol=2e-10,atol=2e-9)
                        expected_lam=wave_matrix(species,w,**args,harmonics=1024)+np.sum(matched-expected,axis=0)
                        np.testing.assert_allclose(out['lam'],expected_lam,rtol=2e-10,atol=2e-9)
        sp=stress_systems()[1][1]
        for mode in ['adaptive','fixed128']:
            for backend in ['legacy','faddeeva']:
                d=self.evaluate(sp,[.73+.13j],mode=mode,backend=backend,kperp=3.)[0][0]
                r=self.evaluate(sp,[.73+.13j],mode=mode,backend=backend,kperp=3.,normalization='rho_ref')[0][0]
                np.testing.assert_allclose(d['chi'],r['chi'],rtol=2e-10,atol=1e-8)

    def test_independent_roots_all_paths_and_growth(self):
        p=mixed.PROTON;e=mixed.ELECTRON;cp=(0.,1.,0.,1.,1.,0.);ce=mixed.COLD_E
        beam=[cp,mixed.replace_species(ce,density=1.1,drift=.2/1.1),(.03,.25,1.,.5,.05,2.)]
        cases=[('kinetic',[p,e],.05-.001j,.06,.05),
               ('warm_reference_mixed',[p,ce],.05-.001j,.06,.05),
               ('cold_reference_mixed',[cp,e],.05-.001j,.06,.05),
               ('parallel_cold_reference',[cp,e],.05-.001j,.06,0.),
               ('growing_beam',beam,.3+.08j,.3,.1)]
        for label,sp,guess,kpar,kperp in cases:
            args=dict(kpar=kpar,kperp=kperp,va=.01)
            def residual(z):
                m=wave_matrix(sp,complex(*z),**args,harmonics=512)
                d=np.linalg.det(m/np.max(abs(m)));return [d.real,d.imag]
            solved=root(residual,[guess.real,guess.imag],tol=1e-10)
            self.assertTrue(solved.success,solved.message);expected=complex(*solved.x)
            for checked in [False,True]:
                blocks,roots,stdout,_=self.evaluate(sp,[expected*.999+1e-5j],roots=True,checked=checked,**args)
                self.assertGreater(roots[0][0],0,stdout);actual=roots[0][1]
                delta=abs(actual-expected);self.assertLess(delta,1e-9)
                b=blocks[0];residual_norm=abs(np.linalg.det(b['lam']/np.max(abs(b['lam']))))
                self.assertLess(residual_norm,1e-9)
                type(self).measurements['max_root_difference']=max(type(self).measurements['max_root_difference'],delta)
                if not checked:
                    legacy=self.evaluate(sp,[expected*.999+1e-5j],roots=True,backend='legacy',**args)[1][0][1]
                    type(self).measurements['roots'].append(dict(case=label,accurate=[actual.real,actual.imag],
                        independent=[expected.real,expected.imag],legacy=[legacy.real,legacy.imag],difference=delta,
                        normalized_residual=residual_norm))
        self.assertGreater(type(self).measurements['roots'][-1]['accurate'][1],.01)

    def test_weak_collisions_keep_corrected_zz(self):
        sp=[mixed.PROTON,mixed.ELECTRON]
        for checked in [False,True]:
            for mode in ['adaptive','fixed128']:
                # nu_s = sqrt(mu_s*tau_s)/(sqrt(2)*Kn) in d_ref.
                w=.73+.4j;kp=.31;kt=3.;va=.02;kn=100.
                block=self.evaluate(sp,[w],mode=mode,checked=checked,kpar=kp,kperp=kt,va=va,
                                    extra=f'collision_type=1\nKn={kn}')[0][0]
                self.assertEqual(block['status'],0)
                for i,s in enumerate(sp):
                    nu=np.sqrt(s[0]*s[1])/(np.sqrt(2)*kn)
                    gamma=block['gamma'][i]
                    expected=zz_moments(s,w,kp,kt,va,nu,len(gamma)-2,gamma=gamma)[2]
                    np.testing.assert_allclose(block['chi'][i,2,2],expected,rtol=1e-10,atol=1e-9)

    def test_failed_trial_recovery_and_root_rejection(self):
        sp=[(.4,1.,1.,1.,1.,0.),(.2,1836.,1.,-1.,1.,0.)]
        # Im(zeta_0)=-30 for reference protons: genuinely unrepresentable Z.
        bad=-30j*.06*np.sqrt(.4);good=.07+.005j
        for checked in [False,True]:
            for mode in ['adaptive','fixed128']:
                for species in [sp,[sp[0],mixed.COLD_E],[(0.,1.,0.,1.,1.,0.),sp[1]]]:
                    blocks,_,stdout,_=self.evaluate(species,[0j,bad,good],mode=mode,checked=checked)
                    self.assertNotEqual(blocks[0]['status'],0,stdout)
                    # The electron-only warm case has a smaller |zeta|; use its own range point.
                    if species[0][0]>0:self.assertNotEqual(blocks[1]['status'],0,stdout)
                    self.assertEqual(blocks[2]['status'],0,stdout)
                    self.assertTrue(np.isfinite(blocks[2]['chi']).all())
                cold_ref=[(0.,1.,0.,1.,1.,0.),sp[1]]
                bad_e=-30j*.06*np.sqrt(sp[1][0]*sp[1][1])
                recovered=self.evaluate(cold_ref,[bad_e,good],mode=mode,checked=checked)[0]
                self.assertNotEqual(recovered[0]['status'],0)
                self.assertEqual(recovered[1]['status'],0)
                _,roots,stdout,_=self.evaluate(sp,[bad],mode=mode,checked=checked,roots=True)
                self.assertLess(roots[0][0],0,stdout);self.assertFalse(np.isfinite(roots[0][1]))

    def test_fully_cold_no_z_calls_and_default_selection(self):
        sp=[(0.,1.,0.,1.,1.,0.),mixed.COLD_E]
        _,_,stdout,_=self.evaluate(sp,[.07+.005j],backend=None)
        self.assertIn('Z_CALLS 0',stdout)
        warm=stress_systems()[1][1];args=dict(kpar=-.31,kperp=3.,va=.02)
        default=self.evaluate(warm,[.73+.13j],backend=None,**args)[0][0]
        explicit=self.evaluate(warm,[.73+.13j],**args)[0][0]
        np.testing.assert_array_equal(default['chi'],explicit['chi'])
        for mode in ['adaptive','fixed128']:
            source=mixed.input_text([mixed.PROTON,mixed.ELECTRON]).replace("z_backend='legacy'","z_backend='invalid'")
            result,_=mixed.run(self.probes[False],source,'')
            self.assertNotEqual(result.returncode,0);self.assertIn('z_backend',result.stdout+result.stderr)

    def test_invalid_map_points_and_refinement(self):
        source=mixed.input_text([mixed.PROTON,mixed.ELECTRON],option=0,va=.01).replace("z_backend='legacy'","z_backend='faddeeva'")
        source+=('&maps\nloggridw=.false.\nloggridg=.false.\nomi=0.\nomf=.1\n'
                 'gami=-.005\ngamf=.005\nnr=20\nni=10\npositive_roots=.true.\n/\n')
        result,files=mixed.run(ROOT/'plume.e',source)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        self.assertIn('singular/nonfinite map points',result.stdout)
        self.assertEqual(result.stdout.count('WARNING: excluded'),1)
        records=[list(map(float,line.split())) for line in files['dispersion_test.roots'].splitlines() if line.strip()]
        expected=.06098715236206473-1.270996796077086e-5j
        self.assertTrue(any(abs(complex(*row[4:6])-expected)<5e-6 for row in records),result.stdout)

    def test_examples(self):
        for name in ['faddeeva_accuracy','legacy_numerics']:
            source=(ROOT/'inputs/example'/f'{name}.in').read_text().replace("dataName='example'","dataName='mixed'")
            result,_=mixed.run(ROOT/'plume.e',source)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            self.assertNotIn('refinement failed',result.stdout+result.stderr)

    def test_python_backend_round_trip(self):
        import sys,tempfile
        sys.path.insert(0,str(ROOT))
        from linfpclib.linfpc import plume_input
        for backend in ['legacy','faddeeva']:
            with tempfile.TemporaryDirectory() as tmp:
                obj=plume_input('mixed')
                obj.set_params(.4,.05,.06,.001,2,0,-1,1,False,False,normalization='d_ref',
                               bessel_mode='fixed128',z_backend=backend)
                self.assertEqual(obj.params['z_backend'],backend)
                self.assertEqual(obj.params['bessel_mode'],'fixed128')
                obj.make_species(.4,1.,1.3,1.,1.,0.)
                obj.make_species(.2,1836.,.8,-1.,1.,0.)
                name=str(Path(tmp)/'roundtrip')
                obj.write_input(name,'roundtrip')
                self.assertIn(f"z_backend='{backend}'",Path(name+'.in').read_text())
                loaded=plume_input('mixed');loaded.load_from_file(name+'.in')
                self.assertEqual(loaded.params['z_backend'],backend)
                self.assertEqual(loaded.params['bessel_mode'],'fixed128')
                with self.assertRaises(ValueError):obj.set_params(.4,.05,.06,.001,2,0,-1,1,False,False,z_backend='invalid')
