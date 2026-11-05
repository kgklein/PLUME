"""Stage 4A: actual scan loops versus standalone and independent dispersion roots."""
import json
import re
import subprocess
import tempfile
import unittest
from pathlib import Path
import numpy as np
from scipy.optimize import root
from mixed_reference import wave_matrix
from plasma_dispersion_audit import MODULES,ROOT
from faddeeva_build import compile_c,c_objects
import test_mixed_dispersion as mixed

P=mixed.PROTON;E=mixed.ELECTRON
CP=(0.,1.,1.,1.,1.,0.);CE=mixed.COLD_E
A=(.03,.25,1.2,.5,.05,.3)
CA=(0.,.25,1.,.5,.05,.3)


def scan(axis,first,last,*,style=0,steps=2,res=2,log=False,tensor=False,heat=False,eigen=False):
    return dict(scan_style=style,scan_type=axis,swi=first,swf=last,ns=steps,nres=res,
                swlog=log,tensor=tensor,heating=heat,eigen=eigen)


def source(species,scans,*,option=1,kpar=.06,kperp=.05,va=.01,backend='faddeeva',mode='adaptive',
           guess=.06-.0001j,normalization='d_ref',extra=''):
    text=mixed.input_text(species,option=option,kpar=kpar,kperp=kperp,va=va,guess=guess,extra=extra)
    text=text.replace('nscan=0',f'nscan={len(scans)}').replace("z_backend='legacy'",f"z_backend='{backend}'")
    text=text.replace("bessel_mode='fixed'",f"bessel_mode='{mode}'\nbessel_diagnostics=.true.\nbessel_nmax=1024")
    scans=[s.copy() for s in scans]
    if normalization=='rho_ref':
        beta=species[0][0];length=np.sqrt(beta*species[0][2])
        text=text.replace("normalization='d_ref'","normalization='rho_ref'")
        for key,val in [('kpar',kpar*length),('kperp',kperp*length),('vtp',va*np.sqrt(beta))]:
            text=re.sub(rf'^{key}=.*$',f'{key}={val:.17g}',text,flags=re.M)
        chunks=text.split('&species_')
        for i,s in enumerate(species):chunks[i+1]=chunks[i+1].replace(f'tauS={s[0]:.17g}',f'tauS={beta/s[0]:.17g}')
        text='&species_'.join(chunks)
        for s in scans:
            if s['scan_style']==0 or s['scan_type']==0:
                s['swi']*=length;s['swf']*=length
            elif s['scan_type']==2:s['swf']*=length
    for i,s in enumerate(scans,1):
        text+=f'&scan_input_{i}\n'
        text+='\n'.join(f'{k}='+('.true.' if v else '.false.') if isinstance(v,bool) else f'{k}={v:.17g}' for k,v in s.items())+'\n/\n'
    return text


def mode_rows(files):
    return {name:np.atleast_2d(np.loadtxt(text.splitlines())) for name,text in files.items()
            if name.endswith('.mode1') and '_tensor_' not in name}


class MixedScanTests(unittest.TestCase):
    measurements={'rows':0,'single_real':0.,'single_imag':0.,'independent_real':0.,'independent_imag':0.,
                  'residual':0.,'rounded_residual':0.,'final_real':0.,'final_imag':0.,'cases':{},'cache':{}}
    @classmethod
    def setUpClass(cls):
        subprocess.run(['make','-j2'],cwd=ROOT,check=True,capture_output=True)
        cls.tmp=tempfile.TemporaryDirectory();cls.addClassCleanup(cls.tmp.cleanup)
        cls.checked=Path(cls.tmp.name)/'checked';cls.checked.mkdir();compile_c(cls.checked)
        flags=['-O0','-g','-fdefault-real-8','-fcheck=all','-ffpe-trap=invalid,zero,overflow']
        base=['gfortran',*flags,'-I',str(cls.checked),'-J',str(cls.checked)]
        result=subprocess.run(base+['-c',*[str(ROOT/'src'/f'{m}.f90') for m in (*MODULES,'fpc')]],
                              cwd=cls.checked,capture_output=True,text=True)
        if result.returncode:raise RuntimeError(result.stderr)
        cls.executables={False:ROOT/'plume.e'};cls.probes={};cls.direct={};cls.failure={}
        for checked in [False,True]:
            directory=cls.checked if checked else ROOT
            objects=[str(directory/f'{m}.o') for m in MODULES]+c_objects(directory)
            compiler=['gfortran',*(flags if checked else ['-O3','-fdefault-real-8']),'-I',str(cls.checked if checked else ROOT/'include')]
            programs=[('scan','mixed_scan_probe'),('direct','mixed_probe'),('failure','scan_failure_probe')]
            if checked:programs.append(('plume','plume'))
            for name,program in programs:
                exe=Path(cls.tmp.name)/f'{name}{checked}'
                file=ROOT/('src' if name=='plume' else 'tests')/f'{program}.f90'
                result=subprocess.run(compiler+[str(file),*objects,*([str(directory/'fpc.o')] if name=='plume' else []),'-o',str(exe)],
                                      capture_output=True,text=True)
                if result.returncode:raise RuntimeError(result.stderr)
                {'scan':cls.probes,'direct':cls.direct,'failure':cls.failure,'plume':cls.executables}[name][checked]=exe

    @classmethod
    def tearDownClass(cls):print('Mixed scan measurements:',json.dumps(cls.measurements))

    def run_scan(self,species,scans,checked=False,probe=False,**kw):
        result,files=mixed.run(self.probes[checked] if probe else self.executables[checked],source(species,scans,**kw))
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        self.assertNotIn('scan root failed',result.stdout,result.stdout)
        if probe:self.assertIn('PHYSICAL_STATE_UNCHANGED',result.stdout)
        rows=mode_rows(files);self.assertTrue(rows,result.stdout)
        return rows,result.stdout,files

    def validate_rows(self,species,rows,*,checked=False,normalization='d_ref',backend='faddeeva',independent=True):
        errors={k:0. for k in ['single_real','single_imag','independent_real','independent_imag','residual']}
        for row in rows:
            self.assertEqual(len(row),6+6*len(species))
            kp,kz,beta,velocity=row[:4];w=complex(*row[4:6])
            self.assertTrue(np.isfinite(w));self.assertNotEqual(kz,0.)
            expected_params=np.array(species).copy()
            if normalization=='rho_ref':
                length=np.sqrt(species[0][0]*species[0][2]);kp/=length;kz/=length
                va=velocity/np.sqrt(beta)
                expected_params[:,0]=beta/expected_params[:,0]
            else:va=velocity
            np.testing.assert_allclose(row[6:].reshape(-1,6),expected_params,rtol=6e-7,atol=1e-13)
            self.assertEqual(beta,species[0][0]);self.assertLess(abs(va-.01),6e-9)
            va=.01 # Actual input, rather than its rounded output representation.
            # A new process solves one point, independently of any scan state.
            text=source(species,[],kpar=kz,kperp=kp,va=va,backend=backend,normalization=normalization)
            guess=w*.9999+1e-8j
            out,_=mixed.run(self.direct[checked],text,f'{guess.real:.17g} {guess.imag:.17g} 1\n')
            self.assertEqual(out.returncode,0,out.stdout+out.stderr)
            blocks,rs=mixed.parse(out.stdout,len(species));self.assertGreater(rs[0][0],0,out.stdout)
            direct=rs[0][1]
            for key,actual,reference in [('single_real',w.real,direct.real),('single_imag',w.imag,direct.imag)]:
                error=abs(actual-reference);errors[key]=max(errors[key],error)
                self.assertLess(error,1e-11+6e-7*abs(reference),(key,row,direct))
            rounded=wave_matrix(species,w,kz,kp,va,512)
            rounded_res=abs(np.linalg.det(rounded/np.max(abs(rounded))))
            type(self).measurements['rounded_residual']=max(type(self).measurements['rounded_residual'],float(rounded_res))
            self.assertLess(rounded_res,1e-8)
            m=blocks[0]['lam'];res=abs(np.linalg.det(m/np.max(abs(m))))
            errors['residual']=max(errors['residual'],float(res));self.assertLess(res,1e-8)
            if independent:
                def residual(z):
                    matrix=wave_matrix(species,complex(*z),kz,kp,va,512)
                    det=np.linalg.det(matrix/np.max(abs(matrix)));return [det.real,det.imag]
                solution=root(residual,[direct.real,direct.imag],tol=1e-10)
                expected=complex(*solution.x)
                self.assertTrue(solution.success or np.linalg.norm(residual(solution.x))<1e-10,solution.message)
                for key,a,b in [('independent_real',w.real,expected.real),('independent_imag',w.imag,expected.imag)]:
                    error=abs(a-b);errors[key]=max(errors[key],error)
                    self.assertLess(error,3e-10+1e-6*abs(b),(key,row,expected))
            type(self).measurements['rows']+=1
        for k,v in errors.items():type(self).measurements[k]=max(type(self).measurements[k],v)
        return errors

    def test_compositions_parallel_scan_both_directions(self):
        systems={'kinetic':[P,E],'cold':[CP,CE],'warm_reference':[P,CE],'cold_reference':[CP,E],
                 'multi_drift':[CP,mixed.replace_species(E,density=1.1,drift=.03),A,CA]}
        for checked in [False,True]:
            for name,sp in systems.items():
                for sign in [-1,1]:
                    scans=[scan(1,sign*.06,sign*.09)]
                    rows,stdout,files=self.run_scan(sp,scans,kpar=sign*.06,checked=checked,probe=True)
                    data=next(iter(rows.values()));self.assertEqual(len(data),3)
                    e=self.validate_rows(sp,data,checked=checked)
                    fields=next(line.split() for line in stdout.splitlines() if line.startswith('FINAL_ROOT'))
                    raw=complex(float(fields[3]),float(fields[4]));kp,kz=data[-1,:2]
                    def reference(z):
                        m=wave_matrix(sp,complex(*z),kz,kp,.01,512)
                        d=np.linalg.det(m/np.max(abs(m)));return [d.real,d.imag]
                    independent=root(reference,[raw.real,raw.imag],tol=1e-10)
                    delta=[abs(raw.real-independent.x[0]),abs(raw.imag-independent.x[1])]
                    self.assertLess(max(delta),2e-9)
                    type(self).measurements['final_real']=max(type(self).measurements['final_real'],delta[0])
                    type(self).measurements['final_imag']=max(type(self).measurements['final_imag'],delta[1])
                    e['full_precision_final_real']=delta[0];e['full_precision_final_imag']=delta[1]
                    e['final_root']=[raw.real,raw.imag];e['independent_root']=independent.x.tolist()
                    type(self).measurements['cases'][f'{name}/{checked}/{sign}']=e
                    if sp[0][0]==0 and any(s[0]>0 for s in sp):
                        self.assertIn('b_star=',files['normalization_test.txt'])
                    self.assertLessEqual(stdout.count('WARNING: Nonzero equilibrium current'),1)

    def test_perpendicular_component_forward_reverse_and_normalizations(self):
        for sp in [[P,E],[P,CE],[CP,E],[CP,CE]]:
            for checked in [False,True]:
                left=self.run_scan(sp,[scan(0,.05,.09)],checked=checked)[0]
                right=self.run_scan(sp,[scan(0,.09,.05)],kperp=.09,checked=checked)[0]
                a=next(iter(left.values()));b=next(iter(right.values()))[::-1]
                np.testing.assert_allclose(a[:,4:6],b[:,4:6],rtol=2e-5,atol=1e-10)
                self.validate_rows(sp,a,checked=checked)
            if all(s[0]>0 for s in sp):
                rho=self.run_scan(sp,[scan(0,.05,.09)],normalization='rho_ref',backend='legacy',mode='fixed128')[0]
                self.validate_rows(sp,next(iter(rho.values())),normalization='rho_ref',backend='legacy')

    def test_two_dimensional_cartesian(self):
        for checked in [False,True]:
            for sp in [[CP,CE],[P,CE],[CP,E],[CP,E,A,CA]]:
                rows,_,_=self.run_scan(sp,[scan(0,.05,.07,steps=1),scan(1,.06,.08,steps=1)],option=2,
                                       checked=checked,probe=True)
                data=next(iter(rows.values()));self.assertEqual(len(data),4)
                np.testing.assert_allclose(data[:,:2],[[.05,.06],[.05,.08],[.07,.06],[.07,.08]])
                self.validate_rows(sp,data,checked=checked)

    def test_compound_vector_angle_magnitude_and_polar_grid(self):
        for checked in [False,True]:
            for s in [scan(0,.08,.09,style=-1),scan(1,55.,0.,style=-1),scan(2,0.,.12,style=-1,log=True)]:
                for sp in [[P,CE],[CP,E],[CP,CE]]:
                    rows,_,_=self.run_scan(sp,[s],checked=checked)
                    self.validate_rows(sp,next(iter(rows.values())),checked=checked)
            rows,_,_=self.run_scan([CP,E],[scan(1,55.,0.,style=-1,steps=1),scan(2,0.,.12,style=-1,steps=1)],
                                   option=2,checked=checked)
            self.validate_rows([CP,E],next(iter(rows.values())),checked=checked)

    def test_parallel_limit_and_initial_point(self):
        for checked in [False,True]:
            for sp in [[P,CE],[CP,E],[CP,CE]]:
                # Input kpar differs: initialize the supplied guess at the first actual scan point.
                rows,_,_=self.run_scan(sp,[scan(1,.06,.09)],kperp=0.,kpar=.4,checked=checked)
                self.validate_rows(sp,next(iter(rows.values())),checked=checked)
                rows,_,_=self.run_scan(sp,[scan(2,0.,.1,style=-1,log=True)],kperp=0.,checked=checked)
                self.validate_rows(sp,next(iter(rows.values())),checked=checked)

    def test_map_initialized_continuation(self):
        for checked in [False,True]:
            for sp in [[P,CE],[CP,E],[CP,CE]]:
                text=source(sp,[scan(1,.06,.09)],kpar=.4).replace('use_map=.false.','use_map=.true.')
                text+=('&maps\nloggridw=.false.\nloggridg=.false.\nomi=.035\nomf=.11\n'
                       'gami=-.001\ngamf=.001\nnr=36\nni=10\npositive_roots=.true.\n/\n')
                result,files=mixed.run(self.executables[checked],text)
                self.assertEqual(result.returncode,0,result.stdout+result.stderr)
                self.assertNotIn('scan root failed',result.stdout)
                data=next(iter(mode_rows(files).values()))
                self.validate_rows(sp,data,checked=checked)

    def test_logarithmic_and_multiple_serial_scans(self):
        for checked in [False,True]:
            rows,_,_=self.run_scan([CP,E],[scan(1,.06,.12,log=True),scan(0,.05,.08)],checked=checked,probe=True)
            self.assertEqual(len(rows),2)
            for values in rows.values():self.validate_rows([CP,E],values,checked=checked)

    def test_backend_harmonic_independence(self):
        for backend in ['legacy','faddeeva']:
            for mode in ['adaptive','fixed','fixed128']:
                rows,_,_=self.run_scan([CP,E,A],[scan(1,.06,.08)],mode=mode,backend=backend,guess=.043-.00002j)
                self.validate_rows([CP,E,A],next(iter(rows.values())),backend=backend,independent=backend=='faddeeva')

    def test_cache_growth_reuse_and_large_argument_roots(self):
        sp=[(1.,1.,1.,1.,1.,0.),CE]
        for checked in [False,True]:
            rows,stdout,files=self.run_scan(sp,[scan(0,10.,50.,steps=4),scan(0,50.,10.,steps=4)],
                                          kperp=10.,kpar=.01,guess=1.05-.005j,checked=checked,probe=True)
            for values in rows.values():self.validate_rows(sp,values,checked=checked)
            stats=next(t for n,t in files.items() if n.endswith('.bessel'))
            fields=next(l.split() for l in stats.splitlines() if l.strip() and l.split()[0]=='1')
            self.assertGreater(int(fields[3]),263-1);self.assertLess(int(fields[2]),int(fields[3]))
            self.assertGreater(int(fields[8]),1);self.assertGreater(int(fields[7]),int(fields[8]))
            self.assertEqual(len([l for l in stats.splitlines() if l.strip() and l.split()[0].isdigit()]),1)
            type(self).measurements['cache'][f'growth/{checked}']=dict(active=int(fields[2]),capacity=int(fields[3]),
                                                                   hits=int(fields[7]),misses=int(fields[8]))
            _,_,files=self.run_scan([CP,E],[scan(1,.06,.09)],checked=checked,probe=True)
            stats=next(t for n,t in files.items() if n.endswith('.bessel'))
            fields=next(l.split() for l in stats.splitlines() if l.strip() and l.split()[0]=='2')
            self.assertEqual(int(fields[8]),1)
            type(self).measurements['cache'][f'parallel/{checked}']=dict(hits=int(fields[7]),misses=int(fields[8]))

    def test_actual_failed_point_then_recovery(self):
        sp=[(1.,1.,1.,1.,1.,0.),CE]
        for checked in [False,True]:
            result,files=mixed.run(self.probes[checked],source(sp,[scan(0,.05,20.,steps=4)],kpar=.01,
                                                                  guess=.043-.00002j))
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            data=next(iter(mode_rows(files).values()))
            self.assertTrue(np.isnan(data[1,4:6]).all(),result.stdout)
            self.assertTrue(np.isfinite(data[2:,4:6]).all(),result.stdout)
            self.assertIn('retaining last valid continuation seed',result.stdout)
            self.assertIn('PHYSICAL_STATE_UNCHANGED',result.stdout)
            self.validate_rows(sp,data[2:],checked=checked)
            type(self).measurements['cases'][f'coarse_scan_recovery/{checked}']=data[:,[0,4,5]].tolist()

    def test_frequency_map_scans(self):
        for checked in [False,True]:
            for sp in [[P,CE],[CP,E],[mixed.replace_species(CP,alpha=-123.),CE]]:
                text=source(sp,[scan(1,.06,.08,steps=1,res=1)],option=4)
                text+=('&maps\nloggridw=.false.\nloggridg=.false.\nomi=.035\nomf=.11\n'
                       'gami=-.001\ngamf=.001\nnr=36\nni=10\npositive_roots=.true.\n/\n')
                result,files=mixed.run(self.executables[checked],text)
                self.assertEqual(result.returncode,0,result.stdout+result.stderr)
                maps=[np.loadtxt(t.splitlines()) for n,t in files.items() if n.endswith('.map')]
                self.assertEqual(len(maps),2)
                for m in maps:
                    self.assertAlmostEqual(min(m[:,0]),.035);self.assertAlmostEqual(max(m[:,0]),.11)
                roots=[np.atleast_2d(np.loadtxt(t.splitlines())) for n,t in files.items() if n.endswith('.roots')]
                self.assertEqual(len(roots),2)
                for values in roots:
                    for row in values:
                        kp,kz=row[:2];w=complex(*row[4:6])
                        def residual(z):
                            m=wave_matrix(sp,complex(*z),kz,kp,.01,512)
                            d=np.linalg.det(m/np.max(abs(m)));return [d.real,d.imag]
                        independent=root(residual,[w.real,w.imag],tol=1e-10)
                        self.assertLess(abs(w-complex(*independent.x)),6e-5*abs(w))

    def test_failure_guard_recovers(self):
        for exe in self.failure.values():
            result=subprocess.run([exe],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            self.assertIn('ATTEMPT 2 FLAG -1',result.stdout)
            self.assertIn('ATTEMPT 3 FLAG -3',result.stdout)
            self.assertIn('ATTEMPT 4 FLAG -2',result.stdout)
            self.assertIn('ATTEMPT 5 FLAG 3',result.stdout)
            self.assertEqual(result.stdout.count('WARNING: scan root failed'),3)

    def test_reject_zero_parallel_and_unsupported_requests(self):
        for checked in [False,True]:
            for sp in [[P,E],[CP,CE],[P,CE],[CP,E]]:
                for s in [scan(1,.06,0.),scan(1,-.06,.06),scan(1,-.1,.2,steps=3,res=1),scan(0,.07,0.,style=-1),scan(1,90.,0.,style=-1)]:
                    result,_=mixed.run(self.executables[checked],source(sp,[s]))
                    self.assertNotEqual(result.returncode,0,result.stdout+result.stderr)
                    self.assertIn('Stage 4A',result.stdout+result.stderr)
                if any(s[0]==0 for s in sp):
                    for s,word in [(scan(2,.4,.8),'wavevector'),(scan(0,.05,.1,heat=True),'Stage 4B'),
                                   (scan(0,.05,.1,eigen=True,heat=True),'Stage 4B')]:
                        result,_=mixed.run(self.executables[checked],source(sp,[s]))
                        self.assertNotEqual(result.returncode,0);self.assertIn(word,result.stdout+result.stderr)

    def test_examples_and_python_wrapper(self):
        import sys
        sys.path.insert(0,str(ROOT))
        from linfpclib.linfpc import plume_input
        for name in ['mixed_scan_warm_reference','mixed_scan_cold_reference','fully_cold_scan']:
            text=(ROOT/'inputs/example'/f'{name}.in').read_text().replace("dataName='example'","dataName='mixed'")
            result,files=mixed.run(ROOT/'plume.e',text)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            self.assertEqual(len(next(iter(mode_rows(files).values()))),9)
        for sp in [[P,CE],[CP,E],[CP,CE]]:
            obj=plume_input('mixed')
            obj.set_params(sp[0][0],.05,.06,.01,len(sp),1,1,1,False,False,normalization='d_ref')
            for species in sp:obj.make_species(*species)
            obj.make_scan(1,0,.06,.09,False,2,2,False,False)
            obj.make_guess(.06,-.0001)
            with tempfile.TemporaryDirectory() as tmp:
                name=str(Path(tmp)/'case')
                obj.write_input(name,'wrapped_scan')
                result,files=mixed.run(ROOT/'plume.e',Path(name+'.in').read_text())
                self.assertEqual(result.returncode,0,result.stdout+result.stderr)
                self.assertEqual(len(next(iter(mode_rows(files).values()))),3)
                obj.params['heating_convention']='legacy'
                obj.scan_inputs[0]['heating']='.true.'
                with self.assertRaisesRegex(ValueError,'Stage 4B'):obj.check_input_for_errors()

    def test_exact_kinetic_scan_compatibility(self):
        # Pinned immediately before Stage 4A, rather than HEAD/self-comparison.
        revision='8331925c7953dcc68f29f14c1d7b78dd9198bd5a'
        directory=Path(self.tmp.name)/'baseline';directory.mkdir()
        tracked=subprocess.check_output(['git','ls-tree','-r','--name-only',revision,'--','src','third_party','Makefile'],cwd=ROOT,text=True).splitlines()
        for name in tracked:
            file=directory/name;file.parent.mkdir(parents=True,exist_ok=True)
            file.write_bytes(subprocess.check_output(['git','show',f'{revision}:{name}'],cwd=ROOT))
        (directory/'include').mkdir(exist_ok=True)
        subprocess.run(['make','-j2'],cwd=directory,check=True,capture_output=True)
        for mode in ['rho_ref','d_ref']:
            text=source([P,E],[scan(1,.06,.09),scan(0,.05,.09)],backend='legacy',mode='fixed128',normalization=mode)
            old,old_files=mixed.run(directory/'plume.e',text.replace("heating_convention='legacy'\n",''))
            new,new_files=mixed.run(ROOT/'plume.e',text)
            self.assertEqual(old.returncode,0,old.stdout+old.stderr)
            self.assertEqual(new.returncode,0,new.stdout+new.stderr)
            names=[n for n in old_files if n.endswith('.mode1')]
            self.assertEqual(len(names),2)
            for name in names:self.assertEqual(old_files[name],new_files[name],name)

    def test_native_state_cold_anisotropy_and_tensor_output(self):
        sp=[CP,E,A,CA]
        for checked in [False,True]:
            rows,_,files=self.run_scan(sp,[scan(1,.06,.09,tensor=True)],checked=checked,probe=True)
            altered=[mixed.replace_species(s,alpha=-123.) if s[0]==0 else s for s in sp]
            new,_,_=self.run_scan(altered,[scan(1,.06,.09)],checked=checked,probe=True)
            np.testing.assert_allclose(next(iter(rows.values()))[:,4:6],next(iter(new.values()))[:,4:6],rtol=1e-8,atol=1e-12)
            tensors=[np.loadtxt(text.splitlines()) for name,text in files.items() if '_tensor_' in name]
            self.assertEqual(len(tensors),1);self.assertEqual(tensors[0].shape,(3,4+18*len(sp)))
            self.assertTrue(np.isfinite(tensors[0]).all())
