"""Independent electromagnetic nullspace/Faraday tests, preserving legacy moments."""
import json
import subprocess
import unittest
from pathlib import Path
import numpy as np
from scipy.optimize import root
from mixed_reference import species_tensor,wave_matrix
from plasma_dispersion_audit import MODULES,ROOT
from faddeeva_build import c_objects
import test_mixed_scans as scans
import test_mixed_dispersion as mixed

P=scans.P;E=scans.E;CP=scans.CP;CE=scans.CE;A=scans.A;CA=scans.CA


def parse(text,nspec):
    blocks=[]
    for line in text.splitlines():
        f=line.split()
        if not f:continue
        if f[0]=='EM':
            blocks.append(dict(status=int(f[1]),w=complex(float(f[2]),float(f[3])),e=np.empty(3,complex),
                               b=np.empty(3,complex),matrix=np.empty((3,3),complex),
                               chi=np.empty((nspec,3,3),complex),legacy_e=np.empty(3,complex),legacy_b=np.empty(3,complex)))
        elif f[0] in ('FIELD','LEGACY_FIELD'):
            prefix='' if f[0]=='FIELD' else 'legacy_'
            blocks[-1][prefix+'e'][int(f[1])-1]=complex(float(f[2]),float(f[3]))
            blocks[-1][prefix+'b'][int(f[1])-1]=complex(float(f[4]),float(f[5]))
        elif f[0]=='LAM':blocks[-1]['matrix'][int(f[1])-1,int(f[2])-1]=complex(float(f[3]),float(f[4]))
        elif f[0]=='CHI':blocks[-1]['chi'][int(f[1])-1,int(f[2])-1,int(f[3])-1]=complex(float(f[4]),float(f[5]))
    return blocks


def distance(actual,reference):
    actual=actual/np.linalg.norm(actual);reference=reference/np.linalg.norm(reference)
    phase=np.vdot(reference,actual)
    return float(np.linalg.norm(actual-reference*phase/abs(phase)))


class MixedEigenfunctionTests(unittest.TestCase):
    measures={'null_residual':0.,'faraday':0.,'reference_angle':0.,'species_error':0.,
              'assembly_error':0.,'polarization':0.,'auxiliary':0.,'scan_fields':0.,'cases':{}}
    @classmethod
    def setUpClass(cls):
        scans.MixedScanTests.setUpClass()
        cls.addClassCleanup(scans.MixedScanTests.tmp.cleanup)
        cls.probes={};cls.matrix_probes={};cls.aux_probes={}
        cls.executables=scans.MixedScanTests.executables
        cls.directory=Path(scans.MixedScanTests.tmp.name)
        for checked in [False,True]:
            directory=scans.MixedScanTests.checked if checked else ROOT
            flags=['-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow'] if checked else ['-O3']
            objects=[str(directory/f'{m}.o') for m in MODULES]+c_objects(directory)
            base=['gfortran','-fdefault-real-8',*flags,'-I',str(directory if checked else ROOT/'include')]
            for name,target in [('em_fields_probe',cls.probes),('em_matrix_probe',cls.matrix_probes),('em_auxiliary_probe',cls.aux_probes)]:
                exe=cls.directory/f'{name}{checked}'
                result=subprocess.run(base+[str(ROOT/'tests'/f'{name}.f90'),*objects,'-o',str(exe)],capture_output=True,text=True)
                if result.returncode:raise RuntimeError(result.stderr)
                target[checked]=exe

    @classmethod
    def tearDownClass(cls):print('Electromagnetic measurements:',json.dumps(cls.measures))

    def evaluate(self,species,w,*,checked=False,request=1,**kw):
        text=scans.source(species,[],**kw)
        result,_=mixed.run(self.probes[checked],text,f'{w.real:.17g} {w.imag:.17g} {request}\n')
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        out=parse(result.stdout,len(species));self.assertEqual(len(out),1,result.stdout)
        return out[0]

    def validate(self,out,species,kpar=.06,kperp=.05,va=.01,reference=True):
        self.assertEqual(out['status'],0,out)
        e=out['e'];b=out['b'];w=out['w'];matrix=out['matrix']
        self.assertEqual(e[0],1.+0j)
        self.assertTrue(np.isfinite(e).all() and np.isfinite(b).all())
        n=np.array([kperp,0.,kpar])/(va*w)
        assembled=np.eye(3)+np.sum(out['chi'],axis=0)+np.outer(n,n)-np.dot(n,n)*np.eye(3)
        err=float(np.max(abs(assembled-matrix)/np.maximum(1.,abs(matrix))))
        type(self).measures['assembly_error']=max(type(self).measures['assembly_error'],err)
        np.testing.assert_allclose(assembled,matrix,rtol=2e-12,atol=1e-7)
        residual=float(np.linalg.norm(matrix@e)/(np.linalg.norm(matrix)*np.linalg.norm(e)))
        type(self).measures['null_residual']=max(type(self).measures['null_residual'],residual)
        self.assertLess(residual,1e-9)
        faraday=np.cross(np.array([kperp,0.,kpar])/va,e)
        fe=float(np.linalg.norm(w*b-faraday)/np.linalg.norm(faraday))
        type(self).measures['faraday']=max(type(self).measures['faraday'],fe)
        self.assertLess(fe,5e-13)
        _,_,vh=np.linalg.svd(matrix)
        angle=distance(e,vh[-1].conj())
        self.assertLess(angle,3e-7)
        if reference:
            independent=np.array([species_tensor(s,w,kpar,kperp,va,512) for s in species])
            species_err=float(np.max(abs(out['chi']-independent)/np.maximum(1.,abs(independent))))
            type(self).measures['species_error']=max(type(self).measures['species_error'],species_err)
            np.testing.assert_allclose(out['chi'],independent,rtol=3e-7,atol=1e-5)
            matrix=wave_matrix(species,w,kpar,kperp,va,512)
            _,_,vh=np.linalg.svd(matrix)
            ref=vh[-1].conj()
            angle=distance(e,ref)
            self.assertLess(angle,3e-7)
            def polar(v):return 2*np.imag(v[0]*v[1].conj())/np.sum(abs(v[:2])**2)
            pe=abs(polar(e)-polar(ref));self.assertLess(pe,3e-7)
            type(self).measures['polarization']=max(type(self).measures['polarization'],float(pe))
        type(self).measures['reference_angle']=max(type(self).measures['reference_angle'],angle)
        return residual,fe,angle

    def test_independent_fields_all_compositions_signed_drifts(self):
        systems={'kinetic':[P,E],'cold':[CP,CE],'warm_reference':[P,CE],'cold_reference':[CP,E],
                 'multi':[CP,mixed.replace_species(E,density=1.1,drift=.03),A,CA]}
        for checked in [False,True]:
            for name,sp in systems.items():
                for sign in [-1,1]:
                    for kp in ([.05] if name=='kinetic' else [0.,.05]):
                        guess=.06-.0001j
                        out=self.evaluate(sp,guess,checked=checked,kpar=sign*.06,kperp=kp)
                        values=self.validate(out,sp,kpar=sign*.06,kperp=kp)
                        type(self).measures['cases'][f'{name}/{checked}/{sign}/{kp}']=dict(
                            w=[out['w'].real,out['w'].imag],e=[[z.real,z.imag] for z in out['e']],
                            b=[[z.real,z.imag] for z in out['b']],null=values[0],faraday=values[1],angle=values[2])

    def test_parallel_circular_polarization(self):
        for checked in [False,True]:
            for kpar in [-.06,.06]:
                out=self.evaluate([CP,CE],.058-.0001j,checked=checked,kpar=kpar,kperp=0.)
                self.validate(out,[CP,CE],kpar=kpar,kperp=0.)
                np.testing.assert_allclose(out['e'],[1.,-1j,0.],rtol=1e-9,atol=1e-9)

    def test_normalization_equivalence_and_legacy_method(self):
        sp=[P,mixed.replace_species(E,drift=.02)]
        for checked in [False,True]:
            for backend in ['legacy','faddeeva']:
                for mode in ['adaptive','fixed128']:
                    native=self.evaluate(sp,.06-.0001j,checked=checked,backend=backend,mode=mode)
                    rho=self.evaluate(sp,native['w'],checked=checked,backend=backend,mode=mode,normalization='rho_ref')
                    self.validate(native,sp,reference=backend=='faddeeva');self.validate(rho,sp,reference=backend=='faddeeva')
                    self.assertLess(distance(native['e'],rho['e']),2e-9)
                    np.testing.assert_allclose(native['b'],rho['b'],rtol=2e-9,atol=1e-9)
                    old=self.evaluate(sp,native['w'],checked=checked,backend=backend,mode=mode,request=2)
                    np.testing.assert_allclose(old['e'],old['legacy_e'],rtol=2e-9,atol=1e-10)
                    np.testing.assert_allclose(old['b'],old['legacy_b'],rtol=2e-9,atol=1e-9)

    def test_growing_and_damped_modes(self):
        sp=[CP,mixed.replace_species(CE,density=1.1,drift=.2/1.1),(.03,.25,1.,.5,.05,2.)]
        for checked in [False,True]:
            for guess in [.3+.08j,.3-.08j]:
                out=self.evaluate(sp,guess,checked=checked,kpar=.3,kperp=.1)
                self.validate(out,sp,kpar=.3,kperp=.1)
                self.assertGreater(abs(out['w'].imag),.01)

    def test_auxiliary_scale_invariance(self):
        sp=[CP,E,A]
        for checked in [False,True]:
            out=self.evaluate(sp,.043-.00002j,checked=checked,mode='fixed128')
            w=out['w']
            text=scans.source(sp,[],mode='fixed128')
            stdin=''.join(f'{w.real:.17g} {w.imag:.17g} {bstar}\n' for bstar in [.01,.1,1.,10.])
            result,_=mixed.run(self.aux_probes[checked],text,stdin)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            for line in result.stdout.splitlines():
                f=line.split()
                if not f or f[0]!='AUX':continue
                self.assertEqual(int(f[1]),0,line)
                values=np.array(list(map(float,f[3:])));values=values[::2]+1j*values[1::2]
                error=distance(values[:3],out['e'])
                type(self).measures['auxiliary']=max(type(self).measures['auxiliary'],error)
                self.assertLess(error,1e-7)
                np.testing.assert_allclose(values[3:],out['b'],rtol=1e-7,atol=1e-8)

    def test_scan_fields_serial_and_two_dimensional(self):
        for checked in [False,True]:
            for sp in [[P,E],[CP,CE],[P,CE],[CP,E],[CP,E,A,CA]]:
                for option,axes in [(1,[scans.scan(1,.06,.08,eigen=True)]),
                                    (2,[scans.scan(0,.05,.07,steps=1),scans.scan(1,.06,.08,steps=1,eigen=True)])]:
                    source=scans.source(sp,axes,option=option,guess=.043-.00002j if len(sp)>2 else .06-.0001j)
                    result,files=mixed.run(self.executables[checked],source)
                    self.assertEqual(result.returncode,0,result.stdout+result.stderr)
                    self.assertNotIn('eigenfunction unavailable',result.stdout,result.stdout)
                    rows=next(iter(scans.mode_rows(files).values()))
                    self.assertEqual(rows.shape[1],18+14*len(sp))
                    for row in rows:
                        kp,kz=row[:2];w=complex(*row[4:6]);values=row[6:18];b=values[:6:2]+1j*values[1:6:2];e=values[6::2]+1j*values[7::2]
                        if any(s[0]==0 for s in sp):self.assertTrue(np.isnan(row[18:18+8*len(sp)]).all())
                        else:self.assertTrue(np.isfinite(row[18:18+8*len(sp)]).all())
                        direct=self.evaluate(sp,w,checked=checked,kpar=kz,kperp=kp)
                        self.assertEqual(direct['status'],0)
                        error=distance(e,direct['e']);self.assertLess(error,3e-6)
                        type(self).measures['scan_fields']=max(type(self).measures['scan_fields'],error)
                        self.assertLess(np.linalg.norm(w*b-np.cross(np.array([kp,0.,kz])/.01,e))/np.linalg.norm(w*b),3e-6)
                        self.validate(direct,sp,kpar=kz,kperp=kp)

    def matrix_request(self,matrix,omega=.1+.01j,checked=False,kpar=.06,kperp=.05,velocity=.01,alpha=1.,mode='d_ref'):
        def complex_text(z):return f'({z.real:.17g},{z.imag:.17g})'
        text=f'{complex_text(omega)} {kpar} {kperp} {velocity} {alpha} {mode}\n'
        text+='\n'.join(' '.join(complex_text(complex(z)) for z in row) for row in matrix)+'\n'
        result=subprocess.run([self.matrix_probes[checked]],input=text,text=True,capture_output=True)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        f=result.stdout.split();v=np.array(list(map(float,f[1:])));v=v[::2]+1j*v[1::2]
        return int(f[0]),v[:3],v[3:]

    def test_cofactor_pivots_phase_and_matrix_scaling(self):
        v=np.array([1.,.3+.8j,.04-.01j]);a=np.array([2.,-1.,0.],complex);a[0]=-np.dot(a[1:],v[1:])
        b=np.array([1j,.2,1.]);b[0]=-np.dot(b[1:],v[1:])
        matrix=np.array([a,b,2*a+3*b])
        for checked in [False,True]:
            for rows in [matrix,matrix[[2,0,1]],matrix*np.array([1e10,1e-8,3.])[:,None],matrix*(1e20+1e20j)]:
                status,e,b=self.matrix_request(rows,checked=checked)
                self.assertEqual(status,0)
                np.testing.assert_allclose(e,v,rtol=1e-12,atol=1e-12)

    def test_invalid_degenerate_gauge_offmode_and_recovery(self):
        for checked in [False,True]:
            for matrix,status in [(np.zeros((3,3),complex),2),(np.ones((3,3),complex),2),
                                  (np.diag([1.,0.,1.]).astype(complex),3),
                                  (np.eye(3,dtype=complex),4),(np.full((3,3),complex(np.nan,0)),1)]:
                code,e,b=self.matrix_request(matrix,checked=checked)
                self.assertEqual(code,status)
                self.assertTrue(np.isnan(e).all() and np.isnan(b).all())
            for omega in [0j,complex(np.nan,0),complex(np.inf,0)]:
                code,e,b=self.matrix_request(np.diag([0.,1.,1.]).astype(complex),omega=omega,checked=checked)
                self.assertEqual(code,1);self.assertTrue(np.isnan(e).all())
            for invalid in [0.,-1.,float('nan'),float('inf')]:
                for keyword in [dict(velocity=invalid),dict(mode='rho_ref',alpha=invalid)]:
                    code,e,b=self.matrix_request(np.diag([0.,1.,1.]).astype(complex),checked=checked,**keyword)
                    self.assertEqual(code,1);self.assertTrue(np.isnan(e).all() and np.isnan(b).all())
            sp=[CP,E]
            result,_=mixed.run(self.probes[checked],scans.source(sp,[]),'0 0 0\n.06 -.0001 1\n')
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            blocks=parse(result.stdout,2);self.assertNotEqual(blocks[0]['status'],0);self.assertEqual(blocks[1]['status'],0)

    def test_fully_cold_perpendicular_x_and_o_gauge(self):
        sp=[CP,CE]
        for checked in [False,True]:
            out=self.evaluate(sp,.2-.001j,checked=checked,kpar=0.,kperp=.2)
            self.validate(out,sp,kpar=0.,kperp=.2)
            # The O mode is exactly Ez-only. Ex=1 is mathematically unavailable.
            plasma=sum(s[4]*s[1]/(s[3]**2*.01**2) for s in sp)
            w=complex(np.sqrt(plasma+(.2/.01)**2))
            ordinary=self.evaluate(sp,w,request=0,checked=checked,kpar=0.,kperp=.2)
            self.assertIn(ordinary['status'],(3,4)) # Finite-precision O-root may leave a tiny nonzero row.
            self.assertTrue(np.isnan(ordinary['e']).all() and np.isnan(ordinary['b']).all())

    def test_cold_reference_anisotropy_and_species_order(self):
        original=[CP,mixed.replace_species(E,density=1.1),A,CA]
        baseline=self.evaluate(original,.043-.00002j)
        for checked in [False,True]:
            for alpha in [0.,-123.,float('nan')]:
                sp=[mixed.replace_species(CP,alpha=alpha),original[1],A,CA]
                out=self.evaluate(sp,baseline['w'],checked=checked)
                self.validate(out,sp)
                self.assertLess(distance(out['e'],baseline['e']),1e-8)
                sp=[sp[0],sp[2],sp[1],sp[3]]
                reordered=self.evaluate(sp,baseline['w'],checked=checked)
                self.validate(reordered,sp)
                self.assertLess(distance(reordered['e'],baseline['e']),1e-8)

    def test_exact_legacy_eigen_scan_files(self):
        revision='8331925c7953dcc68f29f14c1d7b78dd9198bd5a'
        directory=self.directory/'baseline';directory.mkdir()
        names=subprocess.check_output(['git','ls-tree','-r','--name-only',revision,'--','src','third_party','Makefile'],cwd=ROOT,text=True).splitlines()
        for name in names:
            path=directory/name;path.parent.mkdir(parents=True,exist_ok=True)
            path.write_bytes(subprocess.check_output(['git','show',f'{revision}:{name}'],cwd=ROOT))
        (directory/'include').mkdir(exist_ok=True)
        subprocess.run(['make','-j2'],cwd=directory,check=True,capture_output=True)
        for mode in ['rho_ref','d_ref']:
            text=scans.source([P,E],[scans.scan(1,.06,.08,eigen=True),scans.scan(0,.05,.07,eigen=True)],
                              backend='legacy',mode='fixed128',normalization=mode)
            old,expected=mixed.run(directory/'plume.e',text);new,actual=mixed.run(ROOT/'plume.e',text)
            self.assertEqual(old.returncode,0,old.stdout+old.stderr);self.assertEqual(new.returncode,0,new.stdout+new.stderr)
            names=[n for n in expected if n.endswith('.mode1')]
            self.assertEqual(len(names),2)
            for name in names:self.assertEqual(actual[name],expected[name],name)

    def test_failed_scan_fields_and_python_roundtrip(self):
        import sys
        sys.path.insert(0,str(ROOT))
        from linfpclib.linfpc import plume_input,load_plume_sweep
        for checked in [False,True]:
            sp=[(1.,1.,1.,1.,1.,0.),CE]
            text=scans.source(sp,[scans.scan(0,.05,20.,steps=4,eigen=True)],kpar=.01,guess=.043-.00002j)
            result,files=mixed.run(self.executables[checked],text)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            rows=next(iter(scans.mode_rows(files).values()))
            self.assertTrue(np.isnan(rows[1,4:18]).all())
            self.assertTrue(np.isfinite(rows[-1,4:18]).all(),result.stdout)
        obj=plume_input('mixed')
        obj.set_params(0.,.05,.06,.01,2,1,1,1,False,False,normalization='d_ref')
        obj.make_species(*CP);obj.make_species(*E);obj.make_guess(.06,-.0001)
        obj.make_scan(1,0,.06,.08,False,2,2,False,True)
        name=str(self.directory/'written')
        obj.write_input(name,'written')
        result,files=mixed.run(ROOT/'plume.e',Path(name+'.in').read_text())
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        mode=next(t for n,t in files.items() if n.endswith('.mode1'))
        path=self.directory/'field.mode1';path.write_text(mode)
        loaded=load_plume_sweep(str(path),nspec=2,heating=False,eigen=True)
        self.assertTrue(np.isfinite(loaded['exr']).all())
        self.assertTrue(np.isnan(loaded['ux1r']).all() and np.isnan(loaded['n2r']).all())
        obj.scan_inputs[0]['heating']='.true.'
        with self.assertRaisesRegex(ValueError,'Stage 4B.3'):obj.check_input_for_errors()

    def test_single_point_output_and_unavailable_moments(self):
        for sp in [[CP,CE],[P,CE],[CP,E]]:
            out=self.evaluate(sp,.06-.0001j)
            text=scans.source(sp,[],option=-1,guess=out['w'])
            result,_=mixed.run(ROOT/'plume.e',text)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            self.assertIn('ef =',result.stdout);self.assertIn('bf =',result.stdout)
            self.assertIn('species moments, wave energy, heating and FPC unavailable',result.stdout)
            self.assertNotIn('ns =',result.stdout);self.assertNotIn('Ps =',result.stdout)
