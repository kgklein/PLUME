"""Stage 4B.2: independent fluid/Vlasov moments and actual diagnostic workflows."""
import json
import subprocess
import unittest
from pathlib import Path
import numpy as np
from plasma_dispersion_audit import ROOT,MODULES
from faddeeva_build import c_objects
from mixed_reference import species_tensor
from cold_moment_reference import fluid_moments,susceptibility_moments,vlasov_quadrature
import test_mixed_scans as scans
import test_mixed_dispersion as mixed
import test_mixed_eigenfunctions as fields

P,E,CP,CE,A,CA=scans.P,scans.E,scans.CP,scans.CE,scans.A,scans.CA


def ztext(z):
    z=complex(z);return f'({z.real:.17g},{z.imag:.17g})'


def values(text):
    a=np.array(list(map(float,text)));return a[::2]+1j*a[1::2]


def scaled(actual,reference):
    return float(np.max(abs(np.asarray(actual)-reference)/np.maximum(1.,abs(np.asarray(reference)))))


class SpeciesMomentTests(unittest.TestCase):
    measures={'cold_moments':0.,'momentum':0.,'continuity':0.,'cold_current':0.,'kinetic_reference':0.,
              'vlasov':0.,'quadrature_convergence':0.,'normalization':0.,'auxiliary':0.,'species_order':0.,
              'legacy_method':0.,'near_resonances':0.,'scan_moments':0.,'field_change':0.}
    @classmethod
    def setUpClass(cls):
        scans.MixedScanTests.setUpClass()
        cls.addClassCleanup(scans.MixedScanTests.tmp.cleanup)
        cls.directory=Path(scans.MixedScanTests.tmp.name);cls.executables=scans.MixedScanTests.executables
        cls.probes={};cls.kernels={};cls.aux={}
        for checked in [False,True]:
            directory=scans.MixedScanTests.checked if checked else ROOT
            flags=['-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow'] if checked else ['-O3']
            base=['gfortran','-fdefault-real-8',*flags,'-I',str(directory if checked else ROOT/'include')]
            for name,target in [('species_moments_probe',cls.kernels),('moments_eigen_probe',cls.probes),('em_auxiliary_probe',cls.aux)]:
                objects=(['nrtype','vars','species_moments','cold_susceptibility'] if target is cls.kernels else MODULES)
                exe=cls.directory/f'{name}{checked}'
                result=subprocess.run(base+[str(ROOT/'tests'/f'{name}.f90'),*[str(directory/f'{m}.o') for m in objects],
                                           *([] if target is cls.kernels else c_objects(directory)),'-o',str(exe)],capture_output=True,text=True)
                if result.returncode:raise RuntimeError(result.stderr)
                target[checked]=exe

    @classmethod
    def tearDownClass(cls):print('Species moment measurements:',json.dumps(cls.measures))

    def record(self,key,actual,reference):
        err=scaled(actual,reference);type(self).measures[key]=max(type(self).measures[key],err);return err

    def kernel(self,s,w,kz,kp,a,e,b,*,checked=False,chi=None):
        text=f'{0 if chi is None else 1} {ztext(w)} {kz} {kp} {a} '+ ' '.join(str(x) for x in s)+'\n'
        text+=' '.join(ztext(z) for z in np.r_[e,b])+'\n'
        if chi is not None:text+='\n'.join(' '.join(ztext(z) for z in row) for row in chi)+'\n'
        result=subprocess.run([self.kernels[checked]],input=text,text=True,capture_output=True)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        f=next(l.split() for l in result.stdout.splitlines() if l.startswith('MOM '))
        chi_out=np.zeros((3,3),complex)
        for line in result.stdout.splitlines():
            q=line.split()
            if q[0]=='CHI':chi_out[int(q[1])-1,int(q[2])-1]=complex(float(q[3]),float(q[4]))
        return int(f[1]),values(f[2:]),chi_out

    def evaluate(self,species,w=.06-.0001j,checked=False,request=1,**kw):
        stdin=f'{w.real:.17g} {w.imag:.17g} {request}\n'
        if request==3:
            e=np.array([1.,.2+.3j,-.07+.02j]);b=np.cross([kw.get('kperp',.05),0.,kw.get('kpar',.06)],e)/(kw.get('va',.01)*w)
            stdin+=' '.join(ztext(z) for z in np.r_[e,b])+'\n'
        result,_=mixed.run(self.probes[checked],scans.source(species,[],**kw),stdin)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        out=fields.parse(result.stdout,len(species))[0]
        out['moments']=np.empty((len(species),4),complex);out['mstatus']=np.ones(len(species),int)
        out['old']=np.empty((len(species),4),complex)
        for line in result.stdout.splitlines():
            f=line.split()
            if not f:continue
            if f[0]=='MOM':out['moments'][int(f[1])-1]=values(f[3:]);out['mstatus'][int(f[1])-1]=int(f[2])
            if f[0]=='OLD':out['old'][int(f[1])-1]=values(f[2:])
        out['stdout']=result.stdout
        return out

    def validate(self,out,sp,kpar=.06,kperp=.05,a=.01):
        self.assertEqual(out['status'],0,out['stdout']);self.assertTrue((out['mstatus']==0).all(),out['stdout'])
        w=out['w'];e=out['e'];b=out['b']
        for s,actual,chi in zip(sp,out['moments'],out['chi']):
            v,eta=actual[:3],actual[3];delta=w-kpar*s[5]
            residual=abs(delta*eta-(kperp*v[0]+kpar*v[2])/a)/(1+abs(delta*eta)+abs((kperp*v[0]+kpar*v[2])/a))
            type(self).measures['continuity']=max(type(self).measures['continuity'],float(residual));self.assertLess(residual,3e-11)
            if s[0]==0:
                ref=np.r_[fluid_moments(s,w,kpar,kperp,a,e,b)]
                self.record('cold_moments',actual,ref);np.testing.assert_allclose(actual,ref,rtol=5e-12,atol=5e-12)
                current=s[4]/s[3]*(v+np.array([0,0,a*s[5]*eta]))
                expected=-1j*w*a*a*(chi@e)
                error=self.record('cold_current',current,expected);self.assertLess(error,2e-11)
            else:
                independent=species_tensor(s,w,kpar,kperp,a,512)
                ref=np.r_[susceptibility_moments(s,w,kpar,kperp,a,e,independent)]
                self.record('kinetic_reference',actual,ref);np.testing.assert_allclose(actual,ref,rtol=3e-7,atol=2e-7)

    def test_cold_offroot_fluid_equations_and_current(self):
        e=np.array([1,.2+.3j,-.07+.02j]);a=.03
        for checked in [False,True]:
            for s in [CP,CE,CA,(0.,.4,-2.,-.7,.3,-.8)]:
                for u in [0.,-.6,.4]:
                    s=mixed.replace_species(s,drift=u)
                    for kz in [-.3,.3]:
                        for kp in [0.,.2]:
                            for w in [.6+.7j,.6-.07j]:
                                b=np.cross([kp,0,kz],e)/(a*w)
                                status,actual,chi=self.kernel(s,w,kz,kp,a,e,b,checked=checked)
                                self.assertEqual(status,0)
                                v,eta=actual[:3],actual[3];ref=np.r_[fluid_moments(s,w,kz,kp,a,e,b)]
                                self.record('cold_moments',actual,ref);np.testing.assert_allclose(actual,ref,rtol=5e-12,atol=5e-12)
                                g=s[1]/s[3];delta=w-kz*u;rhs=g*(e+np.cross(v,[0,0,1])+a*u*np.cross([0,0,1],b))
                                residual=np.linalg.norm(-1j*delta*v-rhs)/(1+np.linalg.norm(delta*v)+np.linalg.norm(rhs))
                                type(self).measures['momentum']=max(type(self).measures['momentum'],float(residual));self.assertLess(residual,5e-12)
                                continuity=delta*eta-np.dot([kp,0,kz],v)/a
                                self.assertLess(abs(continuity)/(1+abs(delta*eta)),5e-12)
                                current=s[4]/s[3]*(v+np.array([0,0,a*u*eta]))
                                expected=-1j*w*a*a*(chi@e)
                                self.assertLess(self.record('cold_current',current,expected),2e-11)

    def test_cold_temperature_and_density_independence(self):
        e=np.array([1,.1j,.01]);w=.6+.7j;kz=.3;kp=.2;a=.03;b=np.cross([kp,0,kz],e)/(a*w)
        for checked in [False,True]:
            baseline=self.kernel(CP,w,kz,kp,a,e,b,checked=checked)[1]
            for alpha,tau,density in [(float('nan'),float('nan'),1.),(-9.,0.,.03),(0.,-99.,2.)]:
                s=(tau,1.,alpha,1.,density,0.)
                status,actual,_=self.kernel(s,w,kz,kp,a,e,b,checked=checked)
                self.assertEqual(status,0);np.testing.assert_array_equal(actual,baseline)

    def test_cold_poles_invalid_range_and_recovery(self):
        e=np.array([1.,.1j,.01]);b=np.zeros(3,complex)
        for checked in [False,True]:
            for w,s,expected in [(0j,CP,2),(.3+0j,mixed.replace_species(CP,drift=1.),3),(1.+0j,CP,4),(-1.+0j,CP,4),
                                 (.6+.7j,mixed.replace_species(CP,density=0.),6),(.6+.7j,mixed.replace_species(CP,density=-1.),1),
                                 (complex(np.nan,0),CP,1),(.6+.7j,mixed.replace_species(CP,mu=0.),1)]:
                status,result,_=self.kernel(s,w,.3,.2,.03,e,b,checked=checked)
                self.assertEqual(status,expected);self.assertTrue(np.isnan(result).all())
            for aa in [0.,-1.,float('inf'),float('nan'),1e-310]:
                status,result,_=self.kernel(CP,.6+.7j,.3,.2,aa,e,b,checked=checked)
                self.assertNotEqual(status,0);self.assertTrue(np.isnan(result).all())
            status,result,_=self.kernel(CP,.6+.7j,.3,.2,.03,e,b,checked=checked)
            self.assertEqual(status,0);self.assertTrue(np.isfinite(result).all())
            # Both evaluations occur in one process, with traps restored between
            # them, rather than relying solely on separate fresh executables.
            batch=''
            for w in [1.+0j,.6+.7j]:
                batch+=f'0 {ztext(w)} .3 .2 .03 0 1 1 1 1 0\n'
                batch+=' '.join(ztext(z) for z in np.r_[e,b])+'\n'
            result=subprocess.run([self.kernels[checked]],input=batch,text=True,capture_output=True)
            self.assertEqual(result.returncode,0,result.stderr)
            blocks=[line.split() for line in result.stdout.splitlines() if line.startswith('MOM ')]
            self.assertEqual([int(f[1]) for f in blocks],[4,0])
            self.assertTrue(np.isfinite(values(blocks[1][2:])).all())

    def test_near_resonances_without_floors(self):
        e=np.array([1,.2j,.1]);a=.03;s=mixed.replace_species(CP,drift=.4)
        for checked in [False,True]:
            for center in [.12,1.12,-.88]:
                for offset in [1e-3,-1e-5,1e-7]:
                    w=center+1j*offset;b=np.cross([.2,0,.3],e)/(a*w)
                    status,actual,_=self.kernel(s,w,.3,.2,a,e,b,checked=checked)
                    self.assertEqual(status,0)
                    ref=np.r_[fluid_moments(s,w,.3,.2,a,e,b)]
                    self.record('near_resonances',actual,ref)
                    np.testing.assert_allclose(actual,ref,rtol=1e-8,atol=1e-7)

    def test_kinetic_velocity_space_quadrature(self):
        e=np.array([1,.2+.3j,-.07+.02j]);a=.03;w=.6+.7j
        for s in [(.4,1.,1.3,1.,.7,.3),(1e-4,1836.,.8,-1.,1.,-.4),(.2,.25,1.8,.5,.05,.8)]:
            for kz in [-.3,.3]:
                kp=.2;b=np.cross([kp,0,kz],e)/(a*w)
                fine=np.r_[vlasov_quadrature(s,w,kz,kp,a,e,b,24,32)]
                coarse=np.r_[vlasov_quadrature(s,w,kz,kp,a,e,b,16,24)]
                self.assertLess(self.record('quadrature_convergence',coarse,fine),1e-9)
                chi=species_tensor(s,w,kz,kp,a,512)
                for checked in [False,True]:
                    status,actual,_=self.kernel(s,w,kz,kp,a,e,b,checked=checked,chi=chi)
                    self.assertEqual(status,0);self.assertLess(self.record('vlasov',actual,fine),1e-9)
                    np.testing.assert_allclose(actual,fine,rtol=1e-9,atol=1e-10)

    def test_kinetic_current_removable_doppler_limit(self):
        s=mixed.replace_species(E,drift=1.);e=np.array([1.,.1j,.01]);w=.3+0j;a=.03
        chi=species_tensor(s,w,.3,.2,a,512)
        for checked in [False,True]:
            status,actual,_=self.kernel(s,w,.3,.2,a,e,np.zeros(3),checked=checked,chi=chi)
            self.assertEqual(status,0)
            ref=np.r_[susceptibility_moments(s,w,.3,.2,a,e,chi)]
            np.testing.assert_allclose(actual,ref,rtol=1e-12,atol=1e-12)

    def test_mixed_single_point_species_dispatch(self):
        for checked in [False,True]:
            for sp in [[CP,CE],[P,CE],[CP,E],[CP,mixed.replace_species(E,density=1.1,drift=.03),A,CA],[P,CA,E]]:
                for sign in [-1,1]:
                    for kp in [0.,.05]:
                        out=self.evaluate(sp,.043-.00002j if len(sp)>2 else .06-.0001j,checked=checked,kpar=sign*.06,kperp=kp)
                        self.validate(out,sp,kpar=sign*.06,kperp=kp)

    def test_complex_offroot_mixed_kinetic_moments(self):
        sp=[CP,mixed.replace_species(E,drift=-.2),A,CA]
        for checked in [False,True]:
            for w in [.6+.7j,.6-.02j]:
                for kz in [-.3,.3]:
                    out=self.evaluate(sp,w,request=3,checked=checked,kpar=kz,kperp=.2,va=.03)
                    self.validate(out,sp,kpar=kz,kperp=.2,a=.03)

    def test_kinetic_legacy_method_and_normalization_equivalence(self):
        sp=[mixed.replace_species(P,drift=.1),mixed.replace_species(E,drift=-.2),A]
        for checked in [False,True]:
            for backend in ['legacy','faddeeva']:
                for mode in ['adaptive','fixed128']:
                    root=self.evaluate(sp,checked=checked,backend=backend,mode=mode)
                    native=self.evaluate(sp,root['w'],request=2,checked=checked,backend=backend,mode=mode)
                    rho=self.evaluate(sp,native['w'],request=2,checked=checked,backend=backend,mode=mode,normalization='rho_ref')
                    for out in [native,rho]:
                        self.record('legacy_method',out['moments'],out['old'])
                        np.testing.assert_allclose(out['moments'],out['old'],rtol=3e-10,atol=2e-10)
                    self.assertLess(self.record('normalization',native['moments'],rho['moments']),3e-10)

    def test_auxiliary_scale_and_species_order_invariance(self):
        sp=[CP,mixed.replace_species(E,density=1.1),A,CA]
        for checked in [False,True]:
            out=self.evaluate(sp,.043-.00002j,checked=checked,mode='fixed128')
            stdin=''.join(f'{out["w"].real:.17g} {out["w"].imag:.17g} {b}\n' for b in [.01,.1,1.,10.])
            result,_=mixed.run(self.aux[checked],scans.source(sp,[],mode='fixed128'),stdin)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            count=0
            for line in result.stdout.splitlines():
                f=line.split()
                if f and f[0]=='AUX_MOM':
                    count+=1;self.assertEqual(int(f[3]),0)
                    self.assertLess(self.record('auxiliary',values(f[4:]),out['moments'][int(f[2])-1]),2e-7)
            self.assertEqual(count,16)
            reordered=self.evaluate([sp[0],sp[2],sp[1],sp[3]],out['w'],request=0,checked=checked)
            self.record('species_order',out['moments'],reordered['moments'][[0,2,1,3]])
            np.testing.assert_allclose(out['moments'],reordered['moments'][[0,2,1,3]],rtol=2e-7,atol=2e-7)

    def test_serial_and_two_dimensional_scan_moments(self):
        for checked in [False,True]:
            for sp in [[CP,CE],[P,CE],[CP,E],[CP,E,A,CA]]:
                for sign in [-1,1]:
                    for option,axes in [(1,[scans.scan(1,sign*.06,sign*.08,eigen=True)]),
                                        (2,[scans.scan(0,.05,.07,steps=1),scans.scan(1,sign*.06,sign*.08,steps=1,eigen=True)])]:
                        result,files=mixed.run(self.executables[checked],scans.source(sp,axes,option=option,kpar=sign*.06,
                            guess=.043-.00002j if len(sp)>2 else .06-.0001j))
                        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
                        for rows in scans.mode_rows(files).values():
                            for row in rows:
                                w=complex(*row[4:6]);kp,kz=row[:2];n=len(sp)
                                v=values(row[18:18+6*n]).reshape(n,3);eta=values(row[18+6*n:18+8*n]);actual=np.column_stack([v,eta])
                                self.assertTrue(np.isfinite(actual).all(),result.stdout)
                                direct=self.evaluate(sp,w,checked=checked,kpar=kz,kperp=kp)
                                self.assertLess(self.record('scan_moments',actual,direct['moments']),5e-6)
                                self.validate(direct,sp,kpar=kz,kperp=kp)

    def test_species_failure_keeps_fields_and_other_species(self):
        sp=[CP,E,mixed.replace_species(A,density=0.)]
        for checked in [False,True]:
            out=self.evaluate(sp,.06-.0001j,checked=checked)
            self.assertEqual(out['status'],0);self.assertTrue(np.isfinite(out['e']).all())
            np.testing.assert_array_equal(out['mstatus'],[0,0,6])
            self.assertTrue(np.isfinite(out['moments'][:2]).all());self.assertTrue(np.isnan(out['moments'][2]).all())
            self.assertIn('species=3; status=6',out['stdout'])
            good=self.evaluate([CP,E],checked=checked)
            self.assertTrue((good['mstatus']==0).all())

    def test_failed_eigenfunction_and_scan_recovery(self):
        for checked in [False,True]:
            bad=self.evaluate([CP,E],0j,request=0,checked=checked)
            self.assertNotEqual(bad['status'],0);self.assertTrue(np.isnan(bad['moments']).all())
            sp=[(1.,1.,1.,1.,1.,0.),CE]
            result,files=mixed.run(self.executables[checked],scans.source(sp,[scans.scan(0,.05,20.,steps=4,eigen=True)],
                                   kpar=.01,guess=.043-.00002j))
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            rows=next(iter(scans.mode_rows(files).values()))
            self.assertTrue(np.isnan(rows[1,4:34]).all());self.assertTrue(np.isfinite(rows[-1,4:34]).all())

    def test_fields_unchanged_from_stage_4b1(self):
        evidence=json.loads((ROOT/'tests/MIXED_EIGENFUNCTION_RESULTS.json').read_text())['measurements']['cases']
        for name,sp in [('cold',[CP,CE]),('warm_reference',[P,CE]),('cold_reference',[CP,E])]:
            old=evidence[f'{name}/False/1/0.05'];out=self.evaluate(sp,complex(*old['w']),request=0)
            self.assertEqual(out['status'],0)
            expected_e=np.array(old['e'])@np.array([1,1j]);expected_b=np.array(old['b'])@np.array([1,1j])
            self.assertLess(self.record('field_change',out['e'],expected_e),2e-14)
            self.record('field_change',out['b'],expected_b)
            np.testing.assert_allclose(out['b'],expected_b,rtol=2e-14,atol=2e-14)

    def test_single_point_moments_output(self):
        for sp in [[CP,CE],[P,CE],[CP,E],[CP,E,A,CA]]:
            out=self.evaluate(sp,.043-.00002j if len(sp)>2 else .06-.0001j)
            result,files=mixed.run(ROOT/'plume.e',scans.source(sp,[],option=-1,guess=out['w']))
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            self.assertIn('ns =',result.stdout);self.assertIn('Us =',result.stdout)
            self.assertIn('Species moment statuses',result.stdout);self.assertNotIn('Ps =',result.stdout)
            self.assertTrue(any('cold-fluid moments' in text for name,text in files.items() if name.startswith('normalization_')))
