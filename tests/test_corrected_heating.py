"""Corrected diagnostics: independent Stix/fluid/65-digit energy references.

No historical expected numerical values are replaced. Shared helpers choose
legacy explicitly for the prior audit; this module requests corrected explicitly.
"""
import json
import subprocess
import sys
import unittest
from pathlib import Path
import numpy as np
from plasma_dispersion_audit import ROOT,MODULES
from faddeeva_build import c_objects
from heating_energy_reference import energy_reference,absorption,channels,harmonic_reference,mp
from cold_moment_reference import fluid_moments
import test_mixed_scans as scans
import test_mixed_dispersion as mixed
from test_heating_energy_audit import BENCHMARKS,error
sys.path.insert(0,str(ROOT))
from linfpclib.linfpc import load_plume_energy,load_plume_sweep

P,E,CP,CE,A,CA=scans.P,scans.E,scans.CP,scans.CE,scans.A,scans.CA
HISTORICAL=json.loads((ROOT/'tests/MIXED_HEATING_ENERGY_AUDIT_RESULTS.json').read_text())

def source(species,axes=(),**kw):
    return scans.source(species,list(axes),**kw).replace("heating_convention='legacy'","heating_convention='corrected'")

def parse(text,ns):
    points=[]
    for line in text.splitlines():
        f=line.split()
        if not f:continue
        if f[0]=='MODE':points.append(dict(w=complex(*map(float,f[1:])),species=np.zeros((ns,8)),raw=np.zeros((ns,7)),channels=np.zeros((ns,7)),e=np.empty(3,complex),b=np.empty(3,complex),moments=np.empty((ns,4),complex)))
        elif f[0]=='STATUS':points[-1]['status']=list(map(int,f[1:]))
        elif f[0]=='ENERGY':points[-1]['energy']=np.array(list(map(float,f[1:])))
        elif f[0]=='SPEC':points[-1]['species'][int(f[1])-1]=list(map(float,f[2:]))
        elif f[0] in ['RAW_CHANNEL','CHANNEL']:points[-1]['raw' if f[0]=='RAW_CHANNEL' else 'channels'][int(f[1])-1]=list(map(float,f[2:]))
        elif f[0]=='MOM':
            a=np.array(list(map(float,f[2:])));points[-1]['moments'][int(f[1])-1]=a[::2]+1j*a[1::2]
        elif f[0]=='FIELD':
            a=list(map(float,f[2:]));points[-1]['e'][int(f[1])-1]=complex(*a[:2]);points[-1]['b'][int(f[1])-1]=complex(*a[2:])
    return points

class CorrectedHeatingTests(unittest.TestCase):
    measurements={};cases={}
    @classmethod
    def setUpClass(cls):
        scans.MixedScanTests.setUpClass();cls.tmp=scans.MixedScanTests.tmp;cls.directory=Path(cls.tmp.name)
        cls.executables=scans.MixedScanTests.executables;cls.checked=scans.MixedScanTests.checked;cls.probes={};cls.fractions={}
        cls.addClassCleanup(cls.tmp.cleanup)
        for checked in [False,True]:
            directory=cls.checked if checked else ROOT
            flags=['-O0','-g','-fcheck=all','-ffpe-trap=invalid,zero,overflow'] if checked else ['-O3']
            for name,program,target in [('heating','corrected_heating_probe',cls.probes),('fraction','heating_fraction_probe',cls.fractions)]:
                exe=cls.directory/f'{name}{checked}'
                result=subprocess.run(['gfortran','-fdefault-real-8',*flags,'-fno-fast-math','-I',str(directory if checked else ROOT/'include'),
                    str(ROOT/'tests'/f'{program}.f90'),*[str(directory/f'{m}.o') for m in MODULES],*c_objects(directory),'-o',str(exe)],capture_output=True,text=True)
                if result.returncode:raise RuntimeError(result.stderr)
                target[checked]=exe
    @classmethod
    def tearDownClass(cls):print('Corrected heating measurements:',json.dumps(cls.measurements))

    def evaluate(self,sp,w,checked=False,mode=0,e=None,aux=None,harmonic='adaptive',**kw):
        frequencies=[w] if np.isscalar(w) else w
        stdin=''
        for z in frequencies:
            stdin+=f'{z.real:.17g} {z.imag:.17g} {mode}\n'
            if mode==2:stdin+=' '.join(f'({v.real:.17g},{v.imag:.17g})' for v in e)+'\n'
            if mode==3:stdin+=f'{aux:.17g}\n'
        out,_=mixed.run(self.probes[checked],source(sp,mode=harmonic,**kw),stdin)
        self.assertEqual(out.returncode,0,out.stdout+out.stderr)
        result=parse(out.stdout,len(sp));self.assertEqual(len(result),len(frequencies))
        return result[0] if np.isscalar(w) else result

    def case(self,name,checked=False):
        key=f'{name}/{checked}'
        if key not in type(self).cases:
            sp,kw=BENCHMARKS[name];kw=kw.copy();kw.pop('guess',None);a=kw.pop('a',.01)
            z=complex(*HISTORICAL['cases'][f'{name}/False']['w'])
            type(self).cases[key]=self.evaluate(sp,z,checked=checked,va=a,**kw)
        return type(self).cases[key]

    def test_independent_benchmark_energies_and_species(self):
        records={}
        for name in ['kinetic','warm_reference_mixed','cold_reference_mixed','fully_cold','interleaved','weak_growth']:
            for checked in [False,True]:
                p=self.case(name,checked);self.assertEqual(p['status'][:3],[0,0,0],(name,p['status']))
                ref=HISTORICAL['cases'][f'{name}/False'];H=p['energy'][0]
                totalerr=abs(H-ref['energy_reference'])/abs(ref['energy_reference'])
                specieserr=error(p['species'][:,1],ref['species_energy'])
                self.assertLess(totalerr,1e-6);self.assertLess(specieserr,2e-6)
                warm=np.array([s[0]>0 for s in BENCHMARKS[name][0]])
                np.testing.assert_allclose(p['species'][warm,3],np.array(ref['reference_q'])[warm],rtol=3e-8,atol=2e-7)
                self.assertAlmostEqual(H,sum(p['energy'][1:4]),delta=1e-11*abs(H))
                self.assertAlmostEqual(sum(p['energy'][4:7]),1.,delta=1e-12)
                self.assertAlmostEqual(p['energy'][7],p['energy'][1]/(sum(abs(p['energy'][1:3]))+sum(abs(p['species'][:,1]))),delta=1e-14)
                sp=BENCHMARKS[name][0]
                for i,s in enumerate(sp):
                    if s[0]==0:
                        self.assertEqual(p['species'][i,0],1);self.assertTrue(np.isnan(p['species'][i,4]));self.assertTrue(np.isnan(p['channels'][i]).all())
                    else:
                        self.assertEqual(p['species'][i,0],0)
                        self.assertAlmostEqual(p['species'][i,4],p['species'][i,3]/H,delta=1e-12)
                        self.assertLess(p['species'][i,5],1e-6);self.assertLess(p['species'][i,6],1e-6)
                records[f'{name}/{checked}']=dict(H=H,H_species=p['species'][:,1].tolist(),total_relative_error=totalerr,species_scaled_error=specieserr,
                    f_E=p['energy'][4],f_E_abs=p['energy'][7],uncertainty=p['energy'][8],gamma_est=p['energy'][10],derivative_steps=p['species'][:,7].tolist())
        type(self).measurements['benchmarks']=records

    def test_independent_seven_channels_and_omitted_scope(self):
        sp=[P,E,A];e=np.array([1,.3+.7j,.25-.1j]);maximum=0.;scope=[]
        for checked in [False,True]:
            for kz in [-.7,.7]:
                for mode in ['fixed128','adaptive']:
                    p=self.evaluate(sp,1.3+0j,checked=checked,mode=2,e=e,kpar=kz,kperp=3.,harmonic=mode)
                    terms=[harmonic_reference(s,1.3,kz,3.,.01,128)[0] for s in sp]
                    expected=np.array([channels(e,t[127:130]) for t in terms])
                    raw=np.array([sum(absorption(e,t)) for t in terms])
                    maximum=max(maximum,error(p['raw'],expected));np.testing.assert_allclose(p['raw'],expected,rtol=3e-8,atol=1e-8)
                    np.testing.assert_allclose(p['species'][:,3],raw,rtol=3e-8,atol=1e-8)
                    np.testing.assert_allclose(p['raw'][:,:4].sum(axis=1),p['raw'][:,4],rtol=1e-12,atol=1e-8)
                    scope.append(float(abs(raw[0]-sum(p['raw'][0,4:]))/abs(raw[0])))
        self.assertGreater(max(scope),.01);type(self).measurements.update(seven_channels_scaled_error=maximum,omitted_fraction=max(scope))

    def test_fresh_high_precision_drifting_reference(self):
        e=np.array([1.,.2+.3j,.01-.02j]);maximum=0.;raw_error=0.
        systems=[[P,mixed.replace_species(CE,drift=-.07)],
                 [CP,mixed.replace_species(E,drift=.07),A,CA]]
        for sp in systems:
            for kz in [-.31,.31]:
                p=self.evaluate(sp,.73+0j,mode=2,e=e,kpar=kz,kperp=.7,va=.02)
                ref=energy_reference(sp,.73,kz,.7,.02,e)
                self.assertEqual(p['status'][0],0)
                maximum=max(maximum,error(p['species'][:,1],ref['contributions']))
                np.testing.assert_allclose(p['species'][:,1],ref['contributions'],rtol=2e-6,atol=1e-6)
                warm=np.array([s[0]>0 for s in sp]);raw_error=max(raw_error,error(p['species'][warm,3],ref['q'][warm]))
                np.testing.assert_allclose(p['species'][warm,3],ref['q'][warm],rtol=3e-8,atol=1e-8)
        type(self).measurements['fresh_reference_species_scaled_error']=maximum
        type(self).measurements['fresh_reference_Q_scaled_error']=raw_error

    def test_cold_kinetic_energy_identity(self):
        e=np.array([1,.2+.3j,.1-.07j]);maximum=0.
        for checked in [False,True]:
            for a in [.01,.3]:
                p=self.evaluate([CP,CE],.73+0j,mode=2,e=e,checked=checked,va=a,kpar=.31,kperp=.7)
                b=np.cross([.7,0,.31],e)/(a*.73)
                particles=np.array([s[4]/(s[1]*a*a)*np.vdot(fluid_moments(s,.73,.31,.7,a,e,b)[0],fluid_moments(s,.73,.31,.7,a,e,b)[0]).real for s in [CP,CE]])
                np.testing.assert_allclose(p['species'][:,1],particles,rtol=4e-13,atol=1e-8)
                maximum=max(maximum,error(p['species'][:,1],particles));self.assertTrue(np.isnan(p['energy'][9]))
        type(self).measurements['cold_identity_scaled_error']=maximum

    def test_electric_vacuum_substantial(self):
        p=self.case('high_a_cold');self.assertEqual(p['status'][:3],[0,0,0]);self.assertGreater(p['energy'][4],.49);self.assertLess(p['energy'][4],.51)
        type(self).measurements['moderate_a']=dict(H=p['energy'][0],H_E=p['energy'][1],f_E=p['energy'][4],f_E_abs=p['energy'][7])

    def test_weak_damping_and_growth_consistency(self):
        maximum=0.
        for name in ['kinetic','warm_reference_mixed','cold_reference_mixed','interleaved','weak_growth']:
            p=self.case(name);self.assertEqual(p['status'][1],0);difference=abs(p['energy'][10]-p['w'].imag)/abs(p['w'].imag)
            maximum=max(maximum,difference);self.assertLess(difference,2e-4)
        type(self).measurements['weak_balance_relative_error']=maximum

    def test_reactive_applicability_negative_energy(self):
        for name in ['reactive_mixed','reactive_cold']:
            p=self.case(name);self.assertEqual(p['status'][0],0);self.assertEqual(p['status'][1],3)
            self.assertLess(p['energy'][0],0.);self.assertLess(p['energy'][4],0.);self.assertGreater(p['energy'][7],0.)
            self.assertTrue(np.isnan(p['species'][:,4]).all());self.assertTrue(np.isnan(p['channels']).all());self.assertTrue(np.isfinite(p['e']).all())
            self.assertAlmostEqual(sum(p['energy'][4:7]),1.,delta=1e-12)

    def test_signed_fractions_and_cancellation(self):
        for checked in [False,True]:
            result=subprocess.run([str(self.fractions[checked])],input='1 2 -10 3 0\n1 2 -1 -2 0\n1 2 -1 -2.00000001 0.00001\n1 2 3 4 0\n',text=True,capture_output=True)
            self.assertEqual(result.returncode,0,result.stderr);rows=np.loadtxt(result.stdout.splitlines())
            self.assertEqual(rows[0,0],0);self.assertLess(rows[0,2],0.);self.assertGreater(rows[0,4],1.)
            self.assertEqual(rows[1,0],6);self.assertTrue(np.isnan(rows[1,2:5]).all());self.assertTrue(np.isfinite(rows[1,5]))
            self.assertEqual(rows[2,0],6);self.assertEqual(rows[3,0],0);self.assertAlmostEqual(sum(rows[3,2:5]),1.)

    def test_normalization_and_backend_independence(self):
        maximum=0.
        for backend in ['faddeeva','legacy']:
            for harmonic in ['adaptive','fixed128']:
                # bessel mode and probe mode are deliberately separate arguments.
                p=self.evaluate([P,E],self.case('kinetic')['w'],harmonic=harmonic,backend=backend)
                q=self.evaluate([P,E],p['w'],normalization='rho_ref',harmonic=harmonic,backend=backend)
                self.assertEqual(p['status'][0],0);self.assertEqual(q['status'][0],0)
                maximum=max(maximum,error(p['energy'][:8],q['energy'][:8]));np.testing.assert_allclose(p['energy'][:8],q['energy'][:8],rtol=2e-6,atol=1e-10)
        type(self).measurements['normalization_scaled_error']=maximum

    def test_parallel_signed_wavevectors(self):
        for checked in [False,True]:
            for sp in [[P,CE],[CP,E],[CP,E,A,CA]]:
                for kz in [-.06,.06]:
                    p=self.evaluate(sp,.06-.00001j,checked=checked,mode=1,kpar=kz,kperp=0.)
                    self.assertEqual(p['status'][0],0,p['status']);self.assertTrue(np.isfinite(p['energy'][:8]).all())
                    for i,s in enumerate(sp):
                        if s[0]>0:np.testing.assert_allclose(p['raw'][i,:5],0.,atol=1e-10)

    def test_auxiliary_scale_invariance(self):
        sp=[CP,E,A,CA];w=self.case('interleaved')['w'];values=[]
        for b in [.01,.1,1.,10.]:
            p=self.evaluate(sp,w,mode=3,aux=b);self.assertEqual(p['status'][0],0);values.append(p)
        maximum=max(error(p['energy'][:8],values[0]['energy'][:8]) for p in values)
        self.assertLess(maximum,2e-6);type(self).measurements['auxiliary_scaled_error']=maximum

    def test_failures_and_recovery_without_solver_state_changes(self):
        e=np.array([1,.2+.3j,.1-.07j])
        for checked in [False,True]:
            # At complex mode w=i or 1+i the accepted cold dispersion is finite,
            # but the real-frequency energy snapshot has wr=0 or a gyro pole.
            results=self.evaluate([CP,CE],[1j,1+1j,.73+0j],mode=2,e=e,checked=checked,kpar=.31,kperp=.7)
            self.assertEqual(results[0]['status'][0],3);self.assertEqual(results[1]['status'][0],4);self.assertEqual(results[2]['status'][0],0)
            results=self.evaluate([CP,E],[.73+0j,.73+0j],mode=2,e=np.array([0,1,0],complex),checked=checked,kpar=.31,kperp=.7)
            self.assertTrue(all(p['status'][0]==1 for p in results))

    def test_unresolved_derivative_preserves_other_diagnostics(self):
        for checked in [False,True]:
            w=self.case('cold_reference_mixed')['w']
            stdin=f'{w.real:.17g} {w.imag:.17g} 4\n{w.real:.17g} {w.imag:.17g} 0\n'
            out,_=mixed.run(self.probes[checked],source([CP,E]),stdin)
            self.assertEqual(out.returncode,0,out.stdout+out.stderr)
            p,recovered=parse(out.stdout,2)
            self.assertEqual(p['status'][0],5,p['status'])
            self.assertTrue(np.isfinite(p['e']).all());self.assertTrue(np.isfinite(p['b']).all())
            self.assertTrue(np.isnan(p['species'][:,4]).all());self.assertTrue(np.isnan(p['channels']).all())
            self.assertTrue(np.isfinite(p['moments']).all())
            self.assertEqual(recovered['status'][0],0)

    def test_serial_and_two_dimensional_companion_outputs(self):
        records=0;maximum=0.;failed=0;recovered=0.
        for checked in [False,True]:
            for sp in [[P,E],[P,CE],[CP,E],[CP,CE],[CP,E,A,CA]]:
                for option in [1,2]:
                    axes=[scans.scan(1,.06,.062,steps=1,heat=True,eigen=True)] if option==1 else [scans.scan(0,.05,.052,steps=1),scans.scan(1,.06,.062,steps=1,heat=True,eigen=True)]
                    name={tuple(map(tuple,v[0])):k for k,v in BENCHMARKS.items()}.get(tuple(map(tuple,sp)))
                    guess=complex(*HISTORICAL['cases'][f'{name}/False']['w'])
                    out,files=mixed.run(self.executables[checked],source(sp,axes,option=option,guess=guess))
                    self.assertEqual(out.returncode,0,out.stdout+out.stderr)
                    path=self.directory/'scan.energy';path.write_text(next(t for n,t in files.items() if n.endswith('.energy')));energy=load_plume_energy(path)
                    self.assertTrue((energy['root_valid']==1).any(),out.stdout)
                    self.assertEqual(len(energy['record']),2 if option==1 else 4)
                    for i in range(len(energy['record'])):
                        if energy['root_valid'][i]==0:
                            failed+=1
                            self.assertEqual(energy['energy_status'][i],8)
                            self.assertTrue(np.isnan(energy['H_total'][i]))
                            self.assertTrue(np.isnan(energy['Ps_s1'][i]))
                            recovered+=int(np.any(energy['root_valid'][i+1:]==1))
                            continue
                        self.assertEqual(energy['energy_status'][i],0,out.stdout)
                        w=complex(energy['omega_r'][i],energy['gamma'][i]);p=self.evaluate(sp,w,kperp=energy['kperp'][i],kpar=energy['kpar'][i],checked=checked)
                        maximum=max(maximum,error(energy['H_total'][i],p['energy'][0]));self.assertLess(error(energy['H_total'][i],p['energy'][0]),2e-6)
                        name=Path(energy['source_file'][i]).name;self.assertIn(name,files)
                        rows=np.atleast_2d(np.loadtxt(files[name].splitlines()));row=rows[np.argmin(abs(rows[:,1]-energy['kpar'][i])+abs(rows[:,0]-energy['kperp'][i]))]
                        self.assertEqual(len(row),18+22*len(sp));start=18+8*len(sp)
                        np.testing.assert_allclose(row[start:start+len(sp)],p['species'][:,4],rtol=8e-6,atol=1e-12,equal_nan=True,err_msg=str((sp,option,checked,i,energy['kperp'][i],energy['kpar'][i],row[:6])))
                        np.testing.assert_allclose(row[start+len(sp):start+8*len(sp)].reshape(-1,7),p['channels'],rtol=8e-6,atol=1e-12,equal_nan=True,err_msg=str((sp,option,checked,i,energy['kperp'][i],energy['kpar'][i],row[:6])))
                        records+=1
        type(self).measurements['scan_records']=records;type(self).measurements['scan_single_scaled_error']=maximum;type(self).measurements['scan_failed_records']=failed;type(self).measurements['scan_failure_recoveries']=recovered

    def test_single_point_companion_and_default(self):
        text=source([CP,E],option=-1).replace("heating_convention='corrected'\n",'')
        out,files=mixed.run(ROOT/'plume.e',text);self.assertEqual(out.returncode,0,out.stderr)
        self.assertIn('Heating convention: corrected',out.stdout);self.assertTrue(any(n.endswith('.energy') for n in files))

    def test_input_validation_and_legacy_opt_in(self):
        tests=[(source([P,E],option=-1).replace("heating_convention='corrected'","heating_convention='wrong'"),'heating_convention'),
            (source([P,E],option=-1,extra='collision_type=1\nKn=100'),'collision'),
            (scans.source([P,CE],[scans.scan(1,.06,.07,heat=True)]),'legacy'),
            (source([CP,E],[scans.scan(2,.4,.5,heat=True)]),'wavevector')]
        for text,word in tests:
            out,_=mixed.run(ROOT/'plume.e',text);self.assertNotEqual(out.returncode,0);self.assertIn(word,(out.stdout+out.stderr).lower())
        out,files=mixed.run(ROOT/'plume.e',scans.source([P,E],[scans.scan(1,.06,.062,heat=True,eigen=True)],backend='legacy',mode='fixed128'))
        self.assertEqual(out.returncode,0,out.stdout+out.stderr);self.assertFalse(any(n.endswith('.energy') for n in files))

    def test_low_n_false_format_and_unassigned_ew(self):
        for eigen in [False,True]:
            out,files=mixed.run(ROOT/'plume.e',source([CP,E],[scans.scan(1,.06,.062,heat=True,eigen=eigen)],extra='low_n=.false.'))
            self.assertEqual(out.returncode,0,out.stdout+out.stderr)
            rows=next(iter(scans.mode_rows(files).values()));self.assertEqual(rows.shape[1],49 if eigen else 20)
            if eigen:self.assertTrue(np.isnan(rows[:,-1]).all())

    def test_companion_reader_rejects_wrong_schema(self):
        path=self.directory/'bad.energy';path.write_text('# schema=wrong\n')
        with self.assertRaises(ValueError):load_plume_energy(path)
