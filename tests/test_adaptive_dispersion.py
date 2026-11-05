from faddeeva_build import compile_c,c_objects
"""Checkpoint C: production tensors/roots versus converged sums and Stix.

Existing fixtures explicitly select fixed128. These tests exercise the new
(default) adaptive pathway, including diagnostics and both thermal adapters.
"""
from pathlib import Path
import subprocess
import tempfile
import unittest
import numpy as np
from scipy.optimize import root
import test_mixed_dispersion as mixed
from mixed_reference import species_tensor,wave_matrix
ROOT=mixed.ROOT
MODULES=('nrtype','nrutil_trim','vars','plasma_z_faddeeva','plasma_z_legacy','plasma_z_backend','cold_susceptibility','auxiliary_normalization',
         'kinetic_core','harmonic_selection','bessel_cache','adaptive_kinetic','cold_reference','electromagnetic_fields','scan_support',
         'functions','bessel','disprels')

class AdaptiveDispersionTests(unittest.TestCase):
    errors={'truncation':0.,'oracle':0.,'roots':0.}
    @classmethod
    def setUpClass(cls):
        subprocess.run(['make','-j2'],cwd=ROOT,check=True,capture_output=True)
        cls.tmp=tempfile.TemporaryDirectory();cls.addClassCleanup(cls.tmp.cleanup)
        cls.probes={}; cls.core_probes={}
        for checked in [False,True]:
            directory=Path(cls.tmp.name)/str(checked);directory.mkdir()
            flags=['-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow'] if checked else ['-O3']
            sources=[str(ROOT/'src'/f'{m}.f90') for m in MODULES] if checked else [str(ROOT/f'{m}.o') for m in MODULES]
            sources+=compile_c(directory) if checked else c_objects(ROOT)
            exe=directory/'harmonic_probe'
            subprocess.run(['gfortran','-fdefault-real-8',*flags,'-I',str(directory if checked else ROOT/'include'),
                            '-J',str(directory),*sources,str(ROOT/'tests/harmonic_probe.f90'),'-o',str(exe)],
                           cwd=directory,check=True,capture_output=True)
            cls.probes[checked]=exe
            core=directory/'adaptive_core_probe'
            subprocess.run(['gfortran','-fdefault-real-8',*flags,'-I',str(directory if checked else ROOT/'include'),
                            '-J',str(directory),*sources,str(ROOT/'tests/adaptive_core_probe.f90'),'-o',str(core)],
                           cwd=directory,check=True,capture_output=True)
            cls.core_probes[checked]=core

    @classmethod
    def tearDownClass(cls):
        print('Adaptive errors:',cls.errors)

    def evaluate(self,species,frequencies,mode='adaptive',tol=1e-10,fixed=768,cap=1024,
                 checked=False,roots=False,normalization='d_ref',repeats=1,**kw):
        source=mixed.input_text(species,**kw).replace("bessel_mode='fixed'",'')
        source=source.replace('&params',f"&params\nbessel_mode='{mode}'\nbessel_tol={tol}\n"
                              f'bessel_nfixed={fixed}\nbessel_nmax={cap}\nbessel_diagnostics=.true.',1)
        if normalization=='rho_ref':
            beta,alpha=species[0][0],species[0][2]; length=np.sqrt(beta*alpha)
            source=source.replace("normalization='d_ref'","normalization='rho_ref'")
            for key,value in [('kpar',kw.get('kpar',.06)*length),('kperp',kw.get('kperp',.05)*length),
                              ('vtp',kw.get('va',.001)*np.sqrt(beta))]:
                import re
                source=re.sub(rf'^{key}=.*$',f'{key}={value:.17g}',source,flags=re.M)
            chunks=source.split('&species_')
            for i,s in enumerate(species):
                chunks[i+1]=chunks[i+1].replace(f'tauS={s[0]:.17g}',f'tauS={beta/s[0]:.17g}')
            source='&species_'.join(chunks)
        stdin=''.join(f'{w.real} {w.imag} {int(roots)} {repeats}\n' for w in frequencies)
        result,files=mixed.run(self.probes[checked],source,stdin)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        blocks,r=mixed.parse(result.stdout,len(species))
        for b in blocks: self.assertEqual(b['status'],0,result.stdout)
        return blocks,r,result.stdout,files

    def compare(self,a,b,name,tol=2e-7,atol=1e-5):
        err=np.max(abs(a-b)/np.maximum(1.,abs(b)))
        type(self).errors[name]=max(type(self).errors[name],float(err))
        np.testing.assert_allclose(a,b,rtol=tol,atol=atol)

    def test_every_tensor_element_all_paths(self):
        p=(.4,1.,1.7,1.,1.,0.)
        e=(.1,1836.,.8,-1.,1.1,-.2)
        a=(2.,.25,1.3,.5,.05,.4)
        systems=[[p,e],[p,e,a],[p,(0.,1836.,0.,-1.,1.,-.1)],
                 [(0.,1.,0.,1.,1.,0.),e,a],[(0.,1.,0.,1.,1.,0.),(0.,1836.,0.,-1.,1.1,0.),a]]
        for sp in systems:
            for sign in [-1,1]:
                for kp in [.07,3.,30.]:
                    with self.subTest(sp=sp,sign=sign,kp=kp):
                        ws=[.73+.13j,.73-.08j]
                        args=dict(kpar=sign*.31,kperp=kp,va=.02)
                        actual=self.evaluate(sp,ws,**args)[0]
                        reference=self.evaluate(sp,ws,mode='fixed',fixed=1024,**args)[0]
                        for w,out,ref in zip(ws,actual,reference):
                            self.compare(out['chi'],ref['chi'],'truncation',tol=2e-8)
                            independent=np.array([species_tensor(s,w,**args,harmonics=1024) for s in sp])
                            # PLUME's unchanged Z function has a 1e-7 stopping criterion.
                            # At cancellation-sensitive points, measure its error
                            # separately; changing N must not worsen that error.
                            self.compare(out['chi'],independent,'oracle',tol=1e-3)
                            self.assertTrue(np.all(abs(out['chi']-independent) <=
                                abs(ref['chi']-independent)+1e-7*np.maximum(1.,abs(independent))))
                            self.compare(out['lam'],wave_matrix(sp,w,**args,harmonics=1024),'oracle',tol=1e-3)
        # Same physical argument from historical perpendicular gyroradius.
        sp=systems[1]; ws=[.71+.2j]
        left=self.evaluate(sp,ws,kperp=2.)[0][0]
        right=self.evaluate(sp,ws,kperp=2.,normalization='rho_ref')[0][0]
        self.compare(left['chi'],right['chi'],'truncation')

    def test_well_conditioned_independent_tensors(self):
        systems=[[mixed.PROTON,mixed.ELECTRON],
                 [mixed.PROTON,mixed.ELECTRON,mixed.COLD_ALPHA],
                 [(0.,1.,0.,1.,1.,0.),(.2,1836.,.8,-1.,1.,0.)]]
        for sp in systems:
            for sign in [-1,1]:
                args=dict(kpar=sign*.06,kperp=.05,va=.001)
                for w in [.07+.005j,.07-.005j]:
                    actual=self.evaluate(sp,[w],**args)[0][0]['chi']
                    expected=np.array([species_tensor(s,w,**args,harmonics=768) for s in sp])
                    self.compare(actual,expected,'oracle',tol=2e-7)

    def test_checked_and_cap_metadata(self):
        sp=[(1.,1.,1.,1.,1.,0.),(.2,1836.,.8,-1.,1.,.03)]
        out,_,stdout,files=self.evaluate(sp,[.73+.2j]*3,kperp=20.,cap=2,checked=True)
        self.assertEqual(stdout.count('WARNING: Bessel cap:'),2)
        self.assertIn('unresolved=T',files['test.bessel'])
        self.assertIn(' F ',files['test.bessel'])
        self.assertIn('species lambda N capacity',files['test.bessel'])
        self.assertTrue(np.isfinite(out[0]['chi']).all())
        out2=self.evaluate(sp,[.73+.2j],kperp=20.,cap=1024,checked=True)
        self.assertNotIn('WARNING: Bessel cap:',out2[2])
        self.assertIn('unresolved=F',out2[3]['test.bessel'])

    def test_frequency_reuse_and_cold_no_array(self):
        sp=[(0.,1.,0.,1.,1.,0.),(.2,1836.,1.,-1.,1.,.1)]
        _,_,_,files=self.evaluate(sp,[.2+.1j,.3-.02j,.2+.1j])
        rows=[r for r in files['test.bessel'].splitlines()[6:] if not r.startswith('#')]
        self.assertEqual(len(rows),1)
        row=rows[0].split();self.assertEqual(row[0],'2')
        self.assertEqual(row[-3:],['2','1','1'])

    def test_roots_tolerance_and_independent_reference(self):
        cases=[([(1.,1.,1.,1.,1.,0.),(.2,1836.,1.,-1.,1.,0.)],.06,.05,.055-.001j),
               ([(1.,1.,1.,1.,1.,0.),(0.,1836.,0.,-1.,1.,0.)],.06,.05,.055-.001j),
               ([(0.,1.,0.,1.,1.,0.),(.2,1836.,.8,-1.,1.,0.)],.06,.05,.055-.001j),
               ([(0.,1.,0.,1.,1.,0.),(0.,1836.,0.,-1.,1.1,0.),(.01,.25,1.,.5,.05,3.)],.2,.06,.30+.08j)]
        for sp,kpar,kp,guess in cases:
            args=dict(kpar=kpar,kperp=kp,va=.01)
            def residual(z):
                m=wave_matrix(sp,complex(*z),**args,harmonics=768)
                d=np.linalg.det(m/np.max(abs(m)))
                return [d.real,d.imag]
            independent=root(residual,[guess.real,guess.imag],tol=1e-10)
            self.assertTrue(independent.success,independent.message)
            expected=complex(*independent.x)
            ref=self.evaluate(sp,[expected*.999+1e-6j],mode='fixed',roots=True,**args)[1][0]
            for mode,tol,fixed in [('fixed',1e-10,128),('adaptive',1e-6,768),
                                   ('adaptive',1e-10,768),('adaptive',1e-13,768)]:
                blocks,rs,_,_=self.evaluate(sp,[expected*.999+1e-6j],mode=mode,tol=tol,fixed=fixed,roots=True,**args)
                self.assertGreaterEqual(rs[0][0],0)
                self.assertLess(rs[0][0],75)
                w=rs[0][1]
                type(self).errors['roots']=max(type(self).errors['roots'],abs(w-expected))
                np.testing.assert_allclose([w.real,w.imag],[expected.real,expected.imag],rtol=1e-7,atol=1e-9)
                np.testing.assert_allclose([w.real,w.imag],[ref[1].real,ref[1].imag],rtol=1e-7,atol=1e-9)
                m=blocks[-1]['lam']
                self.assertLess(abs(blocks[-1]['det'])/np.max(abs(m))**3,1e-8)

    def test_physical_cache_keys_and_auxiliary_invariance(self):
        # Fields: slot,count,kpar,kperp,a,b*,tau,mu,alpha,Q,D,u,wr,wi,tol.
        base=[2,2,.06,3.,.01,1.,.2,1836.,.8,-1.,1.,.03,.07,.005,1e-10]
        rows=[base.copy() for _ in range(9)]
        for row,b in zip(rows[:4],[.01,.1,1.,10.]): row[5]=b
        rows[4][2]=-.06   # parallel wavevector
        rows[5][10]=.9   # density
        rows[6][11]=-.1  # drift
        rows[7][12]=.08  # complex frequency
        rows[8][1]=3     # species count: preserve still-matching slots
        # Change each perpendicular physical factor: every row invalidates.
        for column,value in [(3,4.),(6,.3),(8,1.4),(7,1000.),(9,-.5)]:
            row=base.copy(); row[column]=value; rows.append(row)
        for checked in [False,True]:
            text=''.join(' '.join(map(str,row))+'\n' for row in rows)
            out=subprocess.run([self.core_probes[checked]],input=text,text=True,capture_output=True)
            self.assertEqual(out.returncode,0,out.stderr)
            blocks=[]; counters=[]
            for line in out.stdout.splitlines():
                f=line.split()
                if not f: continue
                if f[0]=='STATUS':
                    self.assertEqual(f[1],'0'); blocks.append(np.empty((3,3),complex))
                elif f[0]=='CHI': blocks[-1][int(f[1])-1,int(f[2])-1]=complex(float(f[3]),float(f[4]))
                elif f[0]=='COUNTERS': counters.append(tuple(map(int,f[1:])))
            for block in blocks[:4]: self.compare(block,blocks[0],'truncation',tol=5e-9)
            self.assertEqual(counters[8],(8,1,1))
            self.assertEqual(counters[-1],(8,6,6))
            # Cached coefficients agree with a fresh process at each state.
            for row,block in zip(rows[4:],blocks[4:]):
                fresh=subprocess.check_output([self.core_probes[checked]],input=' '.join(map(str,row))+'\n',text=True)
                values=[complex(float(f[3]),float(f[4])) for line in fresh.splitlines()
                        if (f:=line.split()) and f[0]=='CHI']
                np.testing.assert_array_equal(block,np.array(values).reshape(3,3))

    def test_bessel_arithmetic_failure_does_not_poison_cache(self):
        base=[2,2,.06,.05,.01,1.,.2,1836.,.8,-1.,1.,.03,.07,.005,1e-10]
        bad=base.copy();bad[3]=1e100  # finite argument, recurrence order out of range
        for checked in [False,True]:
            records=[base,bad,base]
            out=subprocess.run([self.core_probes[checked]],
                input=''.join(' '.join(map(str,r))+'\n' for r in records),text=True,capture_output=True)
            self.assertEqual(out.returncode,0,out.stderr)
            statuses=[int(line.split()[1]) for line in out.stdout.splitlines() if line.startswith('STATUS')]
            self.assertEqual(statuses,[0,5,0])
            rows=[line.split()[1:] for line in out.stdout.splitlines() if line.startswith('CHI')]
            self.assertEqual(rows[:9],rows[9:])

    def test_order_invariance_and_collisional_shared_body(self):
        sp=[(0.,1.,0.,1.,1.,0.),(.2,1836.,.8,-1.,1.1,.03),(2.,.25,1.3,.5,.05,-.2)]
        a=self.evaluate(sp,[.07+.005j],kperp=3.)[0][0]
        b=self.evaluate([sp[0],sp[2],sp[1]],[.07+.005j],kperp=3.)[0][0]
        self.compare(a['chi'],b['chi'][[0,2,1]],'truncation',tol=5e-9)
        self.compare(a['lam'],b['lam'],'truncation',tol=5e-9)
        warm=[mixed.PROTON,mixed.ELECTRON]
        for checked in [False,True]:
            a=self.evaluate(warm,[.07+.005j],checked=checked,extra='collision_type=1\nKn=100.',kperp=3.)[0][0]
            b=self.evaluate(warm,[.07+.005j],checked=checked,mode='fixed',extra='collision_type=1\nKn=100.',kperp=3.)[0][0]
            self.compare(a['chi'],b['chi'],'truncation',tol=2e-8)

    def test_large_argument_root_beyond_128(self):
        sp=[(1.,1.,1.,1.,1.,0.),(.2,1836.,1.,-1.,1.,0.)]
        args=dict(kpar=.01,kperp=50.,va=.02)
        def residual(z):
            m=wave_matrix(sp,complex(*z),**args,harmonics=1024)
            d=np.linalg.det(m/np.max(abs(m)));return [d.real,d.imag]
        independent=root(residual,[1.01,-.001],tol=1e-10)
        self.assertTrue(independent.success)
        expected=complex(*independent.x)
        for tol in [1e-6,1e-10,1e-13]:
            blocks,rs,_,files=self.evaluate(sp,[expected*.999+1e-6j],roots=True,tol=tol,**args)
            self.assertGreaterEqual(rs[0][0],0);self.assertLess(rs[0][0],75)
            self.assertLess(abs(rs[0][1]-expected),1e-8)
            self.assertLess(abs(blocks[0]['det'])/np.max(abs(blocks[0]['lam']))**3,1e-8)
            first=files['test.bessel'].splitlines()[6].split()
            self.assertGreater(int(first[2]),128)
        reference=self.evaluate(sp,[expected*.999+1e-6j],mode='fixed',fixed=1024,roots=True,**args)[1][0]
        self.assertLess(abs(rs[0][1]-reference[1]),1e-9)
        old=self.evaluate(sp,[expected*.999+1e-6j],mode='fixed',fixed=128,roots=True,**args)[1][0]
        print('Large-lambda roots: independent=',expected,'fixed128=',old,'adaptive=',rs[0])

    def test_python_controls_round_trip_and_default(self):
        import sys
        sys.path.insert(0,str(ROOT))
        from linfpclib.linfpc import plume_input
        wrapped=plume_input('mixed')
        wrapped.set_params(1.,.05,.06,.001,2,0,1,1,False,False,
                           normalization='d_ref',bessel_tol=1e-12,bessel_nmax=1024,
                           bessel_diagnostics=True)
        wrapped.make_species(1.,1.,1.,1.,1.,0.)
        wrapped.make_species(.2,1836.,1.,-1.,1.,0.)
        name=str(Path(self.tmp.name)/'bessel_roundtrip')
        wrapped.write_input(name,'roundtrip')
        loaded=plume_input('mixed');loaded.load_from_file(name+'.in')
        self.assertEqual(loaded.params['bessel_mode'],'adaptive')
        self.assertEqual(loaded.params['bessel_tol'],1e-12)
        self.assertEqual(loaded.params['bessel_nmax'],1024)
        self.assertEqual(loaded.params['bessel_diagnostics'],'.true.')
        source=mixed.input_text([mixed.PROTON,mixed.ELECTRON]).replace("bessel_mode='fixed'",'')
        out,files=mixed.run(self.probes[False],source,'.07 .005 0 1\n')
        self.assertEqual(out.returncode,0,out.stderr)
        self.assertIn('bessel_mode=adaptive',files['test.bessel'])

    def test_adaptive_eigenfunctions_heating_and_scans(self):
        import test_normalization as norm
        probe=Path(self.tmp.name)/'normalization_probe'
        subprocess.run(['gfortran','-O3','-fdefault-real-8','-I',str(ROOT/'include'),
                        str(ROOT/'tests/normalization_probe.f90'),*[str(ROOT/f'{m}.o') for m in MODULES],*c_objects(ROOT),
                        str(ROOT/'fpc.o'),'-o',str(probe)],check=True,capture_output=True)
        species=[(1.,1.,1.3,1.,1.,0.),(.2,1836.,.8,-1.,1.,.02)]
        source=norm.input_text('d_ref',.4,1.3,.05,.06,species)
        adaptive=source.replace("bessel_mode='fixed'","bessel_mode='adaptive'")
        fixed=source.replace("bessel_mode='fixed'","bessel_mode='fixed'\nbessel_nfixed=768")
        a=norm.run_case(probe,adaptive)[0];b=norm.run_case(probe,fixed)[0]
        self.assertEqual(a.returncode,0,a.stderr);self.assertEqual(b.returncode,0,b.stderr)
        norm.NormalizationTests().compare_records(norm.numeric_records(a.stdout),norm.numeric_records(b.stdout))
        for axis in ['kpar','kperp']:
            source=norm.input_text('d_ref',.4,1.3,.05,.06,species,(axis,.09))
            outputs=[]
            for mode in ['adaptive','fixed']:
                text=source.replace("bessel_mode='fixed'",f"bessel_mode='{mode}'\nbessel_nfixed=768")
                out,files=norm.run_case(ROOT/'plume.e',text)
                self.assertEqual(out.returncode,0,out.stdout+out.stderr)
                outputs.append({k:v for k,v in files.items() if '.mode' in k})
            self.assertTrue(outputs[0]);self.assertEqual(outputs[0].keys(),outputs[1].keys())
            for name in outputs[0]:
                a=np.array([[float(x) for x in line.split()] for line in outputs[0][name].splitlines() if line.strip()])
                b=np.array([[float(x) for x in line.split()] for line in outputs[1][name].splitlines() if line.strip()])
                np.testing.assert_allclose(a,b,rtol=2e-6,atol=1e-10)

    def test_invalid_controls(self):
        for extra in ["bessel_mode='bogus'",'bessel_tol=0.','bessel_tol=1.',
                      'bessel_nmax=-1','bessel_nfixed=0','bessel_nfixed=100000']:
            source=mixed.input_text([mixed.PROTON,mixed.ELECTRON],extra=extra)
            out,_=mixed.run(ROOT/'plume.e',source)
            self.assertNotEqual(out.returncode,0)
            self.assertIn('bessel_',out.stderr+out.stdout)
