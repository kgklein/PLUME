"""Corrected rho_ref parameter scans: state, diagnostics, validation and provenance."""
from pathlib import Path
import re
import shlex
import subprocess
import tempfile
import unittest

import numpy as np
import test_mixed_dispersion as mixed
import test_mixed_scans as scans
from heating_energy_reference import energy_reference
from plasma_dispersion_audit import MODULES
from faddeeva_build import c_objects

ROOT=mixed.ROOT
SPECIES=[(1.,1.,1.3,1.,1.,0.),(2.,1836.,.8,-1.,1.,0.)]
BETA=.4
VTP=.01*np.sqrt(BETA)
KP=.05*np.sqrt(BETA*1.3)
KZ=.06*np.sqrt(BETA*1.3)


def assign(text,key,value):
    result,count=re.subn(rf'(?m)^{key}=.*$',f'{key}={value}',text)
    if not count:raise AssertionError(key)
    return result


def source(axes,*,species=SPECIES,option=1,backend='faddeeva',harmonic='adaptive',low=True,roots=1):
    text=mixed.input_text(species,kperp=KP,kpar=KZ,va=VTP,option=option,guess=.061-.000013j,
                          extra=f'low_n=.{str(low).lower()}.\ncollision_type=0')
    text=assign(text,'normalization',"'rho_ref'")
    text=assign(text,'heating_convention',"'corrected'")
    text=assign(text,'z_backend',repr(backend));text=assign(text,'bessel_mode',repr(harmonic))
    text=assign(text,'betap',BETA);text=assign(text,'nscan',len(axes))
    text=assign(text,'nroot_max',roots)
    if roots==2:text+='&guess_2\ng_om=.061\ng_gam=-.000013\n/\n'
    for i,axis in enumerate(axes,1):
        text+=f'&scan_input_{i}\n'
        text+='\n'.join(f'{k}={".true." if v is True else ".false." if v is False else v}' for k,v in axis.items())+'\n/\n'
    return text


def axis(kind,first,last,*,style=0,steps=1,res=2,log=False,heat=True,eigen=True):
    return scans.scan(kind,first,last,style=style,steps=steps,res=res,log=log,heat=heat,eigen=eigen)


def records(text):
    columns=next(line[11:].split() for line in text.splitlines() if line.startswith('# columns: '))
    data=[];state={};populations={}
    for line in text.splitlines():
        values=dict(re.findall(r'(\w+)=\s*([^\s]+)',line))
        if line.startswith('# point_state '):state[int(values.pop('record'))]={k:float(v) for k,v in values.items()}
        elif line.startswith('# species_state '):
            number=int(values.pop('record'));s=int(values.pop('species'))
            populations.setdefault(number,{})[s]=tuple(float(values[k]) for k in ('tauS','muS','alphS','Qs','Ds','vvS'))
        elif line.strip() and not line.startswith('#'):
            fields=shlex.split(line);assert len(fields)==len(columns)
            data.append(dict(zip(columns,fields)))
    return data,state,populations


class ParameterHeatingScanTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        scans.MixedScanTests.setUpClass()
        cls.temp=scans.MixedScanTests.tmp;cls.addClassCleanup(cls.temp.cleanup)
        cls.executables=scans.MixedScanTests.executables.copy()
        cls.directory=Path(cls.temp.name)
        # Deliberate numerical failure at an otherwise valid parameter point.
        # No production test hooks or reference outputs are changed.
        text=(ROOT/'src/disprels.f90').read_text()
        text=text.replace('complex function disp(om)\n',
             'complex function disp(om)\n  use vars, only: betap\n  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan\n',1)
        start=text.index('complex function disp(om)')
        insertion=text.index('  if (has_cold_reference()) then',start)
        text=text[:insertion]+'''  if (betap>.45 .and. betap<.55) then
     disp=cmplx(ieee_value(0.,ieee_quiet_nan),0.)
     return
  endif
'''+text[insertion:]
        probe=cls.directory/'failure_disprels.f90';probe.write_text(text)
        checked=scans.MixedScanTests.checked
        subprocess.run(['gfortran','-O0','-g','-fdefault-real-8','-fcheck=all',
                        '-ffpe-trap=invalid,zero,overflow','-I'+str(checked),'-I'+str(ROOT/'src'),
                        '-J'+str(cls.directory),'-c',str(probe),'-o',str(cls.directory/'failure_disprels.o')],check=True,capture_output=True)
        cls.failure=cls.directory/'failure'
        objects=[str(cls.directory/'failure_disprels.o') if m=='disprels' else str(checked/f'{m}.o') for m in MODULES]
        subprocess.run(['gfortran','-fdefault-real-8',str(ROOT/'src/plume.f90'),'-I'+str(checked),
                        *objects,str(checked/'fpc.o'),*c_objects(checked),'-o',str(cls.failure)],check=True,capture_output=True)

    def run_scan(self,text,checked=False,executable=None):
        result,files=mixed.run(executable or self.executables[checked],text)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        energy=next(v for k,v in files.items() if k.endswith('.energy'))
        data,states,populations=records(energy)
        self.assertIn('# scan_axis ',energy)
        self.assertEqual(len(states),len(data));self.assertEqual(len(populations),len(data))
        for record in data:
            n=int(record['record']);name=Path(record['source_file']).name
            self.assertIn(name,files)
            self.assertEqual(sorted(populations[n]),list(range(1,len(SPECIES)+1)))
        return data,states,populations,files

    def compare_points(self,result,text,checked=False):
        data,states,populations,files=result
        offsets={name:0 for name in files if name.endswith('.mode1') or name.endswith('.mode2')}
        for record in data:
            n=int(record['record']);name=Path(record['source_file']).name
            rows=np.atleast_2d(np.loadtxt(files[name].splitlines()))
            row=rows[offsets[name]];offsets[name]+=1
            self.assertEqual(len(row),18+22*len(SPECIES))
            np.testing.assert_allclose(row[:4],[float(record['kperp']),float(record['kpar']),states[n]['betap'],states[n]['vtp']],rtol=6e-7)
            native=list(populations[n].values())
            np.testing.assert_allclose(row[-6*len(native):],np.array(native).flatten(),rtol=6e-7,atol=1e-13)
            if int(record['root_valid'])==0:
                self.assertEqual(record['energy_status'],'8');self.assertTrue(np.isnan(row[4:6]).all());continue
            single=source([],species=native,option=-1,
                          backend='legacy' if "z_backend='legacy'" in text else 'faddeeva',
                          harmonic='fixed128' if "bessel_mode='fixed128'" in text else 'fixed' if "bessel_mode='fixed'" in text else 'adaptive')
            for key,value in [('betap',states[n]['betap']),('vtp',states[n]['vtp']),('kperp',record['kperp']),
                              ('kpar',record['kpar']),('g_om',record['omega_r']),('g_gam',record['gamma'])]:
                single=assign(single,key,value)
            out,direct=mixed.run(self.executables[checked],single)
            self.assertEqual(out.returncode,0,out.stdout+out.stderr)
            reference=records(next(v for k,v in direct.items() if k.endswith('.energy')))[0][0]
            for key in ('energy_status','applicability','field_status'):
                self.assertEqual(record[key],reference[key],(key,record))
            keys=['H_total','H_E','H_B','H_species_total','f_E','f_B','f_species','f_E_abs']
            keys += [f'{prefix}_s{s}' for s in (1,2) for prefix in ('H','Ps','f','Q','C1','C2','C3','C4','C5','C6','C7')]
            for key in keys:
                np.testing.assert_allclose(float(record[key]),float(reference[key]),rtol=2e-6,atol=1e-10,equal_nan=True,err_msg=key)
            # Normalized powers/channels in the old mode block match the companion.
            expected=[float(record[f'Ps_s{s}']) for s in (1,2)]
            expected += [float(record[f'C{c}_s{s}']) for s in (1,2) for c in range(1,8)]
            np.testing.assert_allclose(row[34:50],expected,rtol=8e-6,atol=1e-10,equal_nan=True)

    def test_all_scalar_controls_and_reference_species(self):
        cases=[axis(2,.4,.44),axis(3,VTP,VTP*1.05,log=True)]
        for s,values in enumerate(SPECIES,1):
            for kind,start in enumerate(values):
                end=start*1.03 if start else .01
                cases.append(axis(kind,start,end,style=s,log=kind<3))
        for checked in (False,True):
            for a in cases:
                with self.subTest(checked=checked,style=a['scan_style'],kind=a['scan_type']):
                    text=source([a]);result=self.run_scan(text,checked)
                    self.assertEqual(len(result[0]),2)
                    self.compare_points(result,text,checked)

    def test_grids_both_axis_orders_flags_and_serial_segments(self):
        pairs=[(axis(0,KP,KP*1.1,heat=False),axis(2,.4,.44)),
               (axis(2,.4,.44,heat=False),axis(2,.8,.85,style=2)),
               (axis(0,1.,1.05,style=1,heat=False),axis(5,0.,.01,style=2)),
               (axis(1,1836.,1900.,style=2,heat=False),axis(4,1.,1.02,style=1)),
               (axis(1,50.,0.,style=-1,heat=False),axis(2,.4,.44)),
               (axis(2,0.,.07,style=-1,heat=False),axis(2,.4,.44)),
               (axis(0,KP*1.1,KZ*1.1,style=-1,heat=False),axis(2,.4,.44))]
        for checked in (False,True):
            for pair in pairs:
                for reverse in (False,True):
                    axes=[a.copy() for a in (pair[::-1] if reverse else pair)]
                    axes[0]['heating']=False;axes[1]['heating']=True
                    text=source(axes,option=2,roots=2)
                    result=self.run_scan(text,checked)
                    self.assertEqual(len(result[0]),8);self.compare_points(result,text,checked)
            axes=[axis(2,.4,.44),axis(2,.8,.85,style=2),axis(0,KP,KP*1.1)]
            text=source(axes);result=self.run_scan(text,checked)
            self.assertEqual(len(result[0]),6);self.compare_points(result,text,checked)
            # The last segment retains changes made by the preceding segments.
            self.assertAlmostEqual(result[1][6]['betap'],.44)
            self.assertAlmostEqual(result[2][6][2][2],.85)

    def test_reverse_traversal_harmonics_and_z_backends(self):
        for backend in ('faddeeva','legacy'):
            for mode in ('adaptive','fixed','fixed128'):
                for first,last in ((.4,.44),(.44,.4)):
                    text=source([axis(2,first,last)],backend=backend,harmonic=mode)
                    self.compare_points(self.run_scan(text),text)

    def test_invalid_trajectories_and_unsupported_scope(self):
        invalid=[(axis(2,'NaN',.44),'finite'),(axis(2,0.,1.),'positive'),(axis(3,VTP,0.),'positive'),
                 (axis(0,1.,0.,style=2),'tauS=0'),(axis(1,1.,-1.,style=2),'positive'),
                 (axis(2,.8,-.1,style=2),'positive'),(axis(3,-1.,1.,style=2),'cross zero'),
                 (axis(3,-1.,0.,style=2),'reach zero'),(axis(4,1.,-.1,style=2),'nonnegative'),
                 (axis(5,-.01,.01,style=2,log=True),'logarithmic'),
                 (axis(2,.4,.44,steps=0),'positive ns'),(axis(2,.4,.44,res=0),'positive ns'),
                 (axis(2,.4,.44,steps=2147483647,res=2),'integer range'),
                 (axis(2,.4,.44,style=3),'species/style'),(axis(0,.4,.44,style=3),'species/style'),
                 (axis(6,.4,.44,style=2),'scan type')]
        for checked in (False,True):
            for a,message in invalid:
                out,_=mixed.run(self.executables[checked],source([a]))
                self.assertNotEqual(out.returncode,0);self.assertIn(message,out.stdout+out.stderr)
            for axes in ([axis(2,.4,.44),axis(2,.4,.45)],
                         [axis(1,-.06,.06,steps=3),axis(2,.4,.44)]):
                out,_=mixed.run(self.executables[checked],source(axes,option=2))
                self.assertNotEqual(out.returncode,0)
            # An ignored first-axis flag cannot evade the workflow restriction.
            axes=[axis(2,.4,.44,heat=False),axis(0,KP,KP*1.1)]
            for change in ("normalization='d_ref'",'collision_type=1'):
                text=source(axes,option=2)
                text=assign(text,*change.split('=',1))
                if change.startswith('normalization'):text=assign(text,'betap',1.)
                out,_=mixed.run(self.executables[checked],text)
                self.assertNotEqual(out.returncode,0)

    def test_absent_population_and_warning_only_equilibria(self):
        for a in (axis(4,1.,0.,style=2),axis(3,-1.,-1.05,style=2),axis(5,0.,.01,style=1)):
            text=source([a]);result=self.run_scan(text)
            self.compare_points(result,text)
            warnings=result[3]['test.warn']
            self.assertIn('WARNING: parameter scan',warnings)
            if a['scan_type']==4:
                self.assertEqual(result[0][-1]['status_s2'],'2.0000000000000000E+000')

    def test_failure_recovery_and_comments_on_failed_records(self):
        text=source([axis(2,.4,.6,steps=2,res=1)])
        result=self.run_scan(text,executable=self.failure)
        self.assertEqual([r['root_valid'] for r in result[0]],['1','0','1'])
        self.compare_points(result,text)

    def test_output_flags_low_n_and_existing_energy_reader(self):
        from linfpclib.linfpc import load_plume_energy
        for low in (False,True):
            for eigen in (False,True):
                axes=[axis(2,.4,.44,heat=False),axis(0,KP,KP*1.1,eigen=eigen)]
                text=source(axes,option=2,low=low)
                result=self.run_scan(text)
                rows=np.atleast_2d(np.loadtxt(next(v for k,v in result[3].items() if k.endswith('.mode1')).splitlines()))
                expected=(62 if eigen else 34) if low else (49 if eigen else 20)
                self.assertEqual(rows.shape,(4,expected))
                path=self.directory/'parameter.energy'
                path.write_text(result[3]['test.energy'])
                loaded=load_plume_energy(path)
                self.assertEqual(len(loaded['record']),4)
                self.assertIn('H_total',loaded)
                if not low and eigen:self.assertTrue(np.isnan(rows[:,-1]).all())
        # First-axis flags never request grid diagnostics by themselves.
        text=source([axis(2,.4,.44),axis(0,KP,KP*1.1,heat=False)],option=2)
        out,files=mixed.run(self.executables[False],text)
        self.assertEqual(out.returncode,0,out.stderr)
        self.assertNotIn('test.energy',files)

    def test_tensor_snapshots_after_parameter_heating(self):
        from mixed_reference import species_tensor
        axes=[axis(2,.4,.44,heat=False),axis(2,.8,.85,style=2)]
        for a in axes:a['tensor']=True
        text=source(axes,option=2)
        data,states,populations,files=self.run_scan(text)
        rows=np.atleast_2d(np.loadtxt(next(v for k,v in files.items() if '_tensor_' in k).splitlines()))
        self.assertEqual(rows.shape,(4,40))
        components=[(0,0),(0,1),(0,2),(1,1),(1,2),(2,2)]
        for r,row in zip(data,rows):
            n=int(r['record']);beta=states[n]['betap'];sp=np.array(list(populations[n].values()))
            length=np.sqrt(beta*sp[0,2]);sp[:,0]=beta/sp[:,0]
            w=complex(float(r['omega_r']),float(r['gamma']))
            chi=[species_tensor(s,w,float(r['kpar'])/length,float(r['kperp'])/length,
                                states[n]['vtp']/np.sqrt(beta)) for s in sp]
            expected=np.array([[c[i,j] for c in chi] for i,j in components])
            actual=row[4:28].reshape(6,2,2)
            np.testing.assert_allclose(actual[:,:,0]+1j*actual[:,:,1],expected,rtol=2e-6,atol=1e-7)

    def test_independent_energy_for_a_parameter_grid(self):
        text=source([axis(2,.4,.44,heat=False),axis(2,.8,.85,style=2)],option=2)
        result=self.run_scan(text)
        data,states,populations,files=result
        rows=np.atleast_2d(np.loadtxt(next(v for k,v in files.items() if k.endswith('.mode1')).splitlines()))
        for r,row in zip(data,rows):
            n=int(r['record']);beta=states[n]['betap'];sp=np.array(list(populations[n].values()))
            length=np.sqrt(beta*sp[0,2]);sp[:,0]=beta/sp[:,0]
            a=states[n]['vtp']/np.sqrt(beta)
            e=row[12:18:2]+1j*row[13:18:2]
            ref=energy_reference(sp,float(r['omega_r']),float(r['kpar'])/length,float(r['kperp'])/length,a,e)
            self.assertEqual(r['energy_status'],'0')
            np.testing.assert_allclose(float(r['H_total']),ref['energy'],rtol=2e-6)
            for s in (1,2):
                np.testing.assert_allclose(float(r[f'Ps_s{s}']),ref['q'][s-1]/ref['energy'],rtol=8e-6,atol=1e-10)


if __name__=='__main__':unittest.main()
