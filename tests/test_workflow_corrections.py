"""Stage 5B.0: schema markers, real files, scalar roundtrips and map characterization.

All solver outputs live in TemporaryDirectory; no golden numerical data change.
"""
import contextlib
import io
import subprocess
import sys
import tempfile
import unittest
import warnings
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from contours.plot_plume_resonances import read_plume_input,read_dispersion_output
from linfpclib.linfpc import (plume_input,loadlinfpccepar,loadlinfpcceperp,
    loadlinfpccart,loadlinfpc3d,load_plume_energy)
from linfpclib._namelist import read_scalar_namelists
import test_mixed_dispersion as mixed
from test_mixed_scans import source,scan,P,E,CP,CE


def quiet_write(obj,path):
    with contextlib.redirect_stdout(io.StringIO()),warnings.catch_warnings():
        warnings.simplefilter('ignore')
        obj.write_input(str(path),obj.params.get('outputName','roundtrip'))


class ContourSchemaTests(unittest.TestCase):
    def test_unique_marker_layouts_and_unavailable_values(self):
        labels=('TTD_yy','TTD_yz','LD_zy','LD_zz','n0','n+1','n-1')
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);inp=p/'case.in';file=p/'marker.mode1'
            for n in (1,2,4):
                for eigen in (False,True):
                    for heat in (False,True):
                        for low in (False,True):
                            sp=[P]*n
                            text=source(sp,[scan(0,.05,.06,eigen=eigen,heat=heat)])
                            text=text.replace('nscan=1','low_n='+('.true.' if low else '.false.')+'\nnscan=1')
                            inp.write_text(text);meta=read_plume_input(inp)
                            prefix=18+8*n if eigen else 6
                            count=(8*n if low else n) if heat else 0
                            params=prefix+count
                            row=np.arange(params+6*n+int(eigen and heat and not low),dtype=float)+1000
                            if heat:row[prefix]=np.nan
                            np.savetxt(file,[row,row])
                            result=read_dispersion_output(file,meta,1)
                            self.assertEqual(result.eigen_available,eigen)
                            self.assertEqual(result.heating_available,heat)
                            self.assertEqual(result.params_start,params)
                            np.testing.assert_array_equal(result.data[:,params:params+6*n],np.tile(row[params:params+6*n],(2,1)))
                            if heat:
                                self.assertTrue(np.isnan(result.data[0,result.columns['P_species_1']]))
                                if low:
                                    for s in range(n):
                                        for c,label in enumerate(labels):
                                            idx=result.columns[f'P_species_{s+1}_{label}']
                                            self.assertEqual(idx,prefix+n+7*s+c)
                                            self.assertEqual(result.data[0,idx],row[idx])
                            else:
                                self.assertIsNone(result.heat_start)
                                self.assertFalse(any(k.startswith('P_species') for k in result.columns))

    def test_ambiguous_or_mismatched_layout_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);inp=p/'case.in';out=p/'case.mode1'
            inp.write_text(source([P,E],[scan(0,.05,.06,eigen=True),scan(1,.06,.07,heat=True)]))
            meta=read_plume_input(inp);np.savetxt(out,[np.arange(46)]) # eigen-only, including parameters
            with self.assertRaisesRegex(ValueError,'ambiguous.*scan_index'):
                read_dispersion_output(out,meta,1)
            self.assertFalse(read_dispersion_output(out,meta,1,scan_index=1).heating_available)
            with self.assertRaisesRegex(ValueError,'expected .*columns'):
                read_dispersion_output(out,meta,1,scan_index=2)
            np.savetxt(out,[np.arange(47)])
            with self.assertRaisesRegex(ValueError,'extra columns are not inferred'):
                read_dispersion_output(out,meta,1,scan_index=1)

    def test_two_dimensional_uses_second_axis_flags(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);inp=p/'case.in';out=p/'case.mode1'
            inp.write_text(source([P,E],[scan(0,.05,.06,heat=True),scan(1,.06,.07,eigen=True)],option=2))
            np.savetxt(out,[np.arange(46)])
            data=read_dispersion_output(out,read_plume_input(inp),1)
            self.assertTrue(data.eigen_available);self.assertFalse(data.heating_available)
            self.assertEqual(data.params_start,34)

    def test_nonrectangular_legacy_orphan_records_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);inp=p/'case.in';out=p/'case.mode1'
            inp.write_text(source([P,E],[scan(0,.05,.06,eigen=True,heat=True)]))
            out.write_text(' '.join(map(str,range(62)))+'\n0.0\n')
            with self.assertRaisesRegex(ValueError,'nonrectangular.*orphan'):
                read_dispersion_output(out,read_plume_input(inp),1)


class WrapperRoundtripTests(unittest.TestCase):
    def test_writeout_and_map_are_independent(self):
        for use_map in (False,True):
            for writeout in (False,True):
                obj=plume_input('test');obj.set_params(1,.1,.2,.01,2,0,1,1,use_map,writeout)
                self.assertEqual(obj.params['use_map'],'.true.' if use_map else '.false.')
                self.assertEqual(obj.params['writeOut'],'.true.' if writeout else '.false.')

    def test_map_switches_and_tensor_setters_are_independent(self):
        for positive in (False,True):
            for logarithmic in (False,True):
                obj=plume_input('test');obj.set_maps(logarithmic,.1,1,-.1,-1,positive,
                                                     loggridg=True,nr=7,ni=9)
                self.assertEqual(obj.maps['positive_roots'],'.true.' if positive else '.false.')
                self.assertEqual(obj.maps['loggridw'],'.true.' if logarithmic else '.false.')
                self.assertEqual(obj.maps['loggridg'],'.true.')
                self.assertEqual((obj.maps['nr'],obj.maps['ni']),(7,9))
        obj.make_scan(0,0,.1,.2,False,2,1,False,True,tensor=True)
        self.assertEqual(obj.scan_inputs[0]['tensor'],'.true.')

    def roundtrip(self,text):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);original=p/'original.in';copy=p/'copy'
            original.write_text(text);obj=plume_input('constructor_name');obj.load_from_file(original)
            quiet_write(obj,copy)
            old=read_scalar_namelists(text);new=dict(read_scalar_namelists(copy.with_suffix('.in').read_text()))
            for group,values in old:
                for name,value in values.items():self.assertEqual(new[group][name],value,(group,name))
            second=plume_input('new_name');second.load_from_file(copy.with_suffix('.in'))
            self.assertEqual(obj.params,second.params)
            self.assertEqual(obj.species,second.species);self.assertEqual(obj.maps,second.maps)
            self.assertEqual(obj.scan_inputs,second.scan_inputs);self.assertEqual(obj.fpc,second.fpc)
            self.assertEqual(obj.guesses,second.guesses);self.assertEqual(obj.dataname,second.dataname)
            # A repeated load replaces indexed groups instead of appending them.
            second.load_from_file(copy.with_suffix('.in'));self.assertEqual(obj.species,second.species)

    def test_scalar_optional_controls_map_scan_fpc_roundtrips(self):
        scans=[scan(0,.05,.07,eigen=True,heat=True,tensor=True)]
        for species in ([P,E],[P,CE],[CP,E],[CP,E,P,CE]):
            text=source(species,scans).replace("heating_convention='legacy'","heating_convention='corrected'")
            self.roundtrip(text.replace('nscan=1','low_n=.false.\nnscan=1'))
        text=source([P,E],[scan(0,.05,.07,eigen=True,tensor=True)],option=4)
        text+='&maps loggridw=.true., loggridg=.true., positive_roots=.false., omi=.1, omf=1, gami=-1d-4, gamf=-1d-2, nr=7, ni=9 /\n'
        self.roundtrip(text)
        for option in (6,7):
            text=source([P,E],[],option=option)
            text+='&fpc vperpmin=0, vperpmax=3, vparmin=-4, vparmax=5, delv=.25, vxmin=-1, vxmax=2, vymin=-3, vymax=4, vzmin=-5, vzmax=6, elecdircontribution=.2 /\n'
            self.roundtrip(text)

    def test_quoted_paths_case_and_d_exponents(self):
        text=source([P,E],[]).replace("dataName='mixed'","DATANAME='a path/with ! punctuation'")
        text=text.replace("outputName='test'","outputName='it''s a mode'").replace('vtp=0.001','VTP=1d-3')
        self.roundtrip(text)

    def test_logical_scalar_spellings(self):
        for true,false in [('T','F'),('.TRUE.','.FALSE.'),('.t.','.f.')]:
            values=dict(read_scalar_namelists(f'&maps loggridw={true} loggridg={false} /'))['maps']
            self.assertEqual(values,dict(loggridw='.true.',loggridg='.false.'))

    def test_unknown_or_unsupported_controls_fail_before_overwriting(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'input.in'
            for bad in ('&maps nr_typo=3 /','&radial_input n_radius=3 /','&fpc vx_bounds(1)=0 /'):
                path.write_text(source([P,E],[])+bad)
                with self.assertRaises(ValueError):plume_input('test').load_from_file(path)
            path.write_text(source([P,E],[]));obj=plume_input('test');obj.load_from_file(path)
            obj.maps['critical_unknown']=True;original=path.read_text()
            with self.assertRaisesRegex(ValueError,'refusing to discard'):
                quiet_write(obj,path.with_suffix(''))
            self.assertEqual(path.read_text(),original)

    def test_repository_fixtures_roundtrip(self):
        names=['normalization_rho_ref','normalization_d_ref','mixed_scan_cold_reference',
               'mixed_scan_warm_reference','heating_corrected_cold','heating_corrected_kinetic']
        for name in names:
            with self.subTest(example=name):self.roundtrip((ROOT/'inputs/example'/f'{name}.in').read_text())


class FPCMissingValueTests(unittest.TestCase):
    header='tau bi kpar kperp vti mu omega.r omega.i\n1 1 .1 .1 .01 1 .1 -.001\n'

    def test_gyro_missing_zero_and_nonfinite_are_distinct(self):
        with tempfile.TemporaryDirectory() as d:
            file=Path(d)/'test.cpar.specie01.mode01'
            file.write_text(self.header+'vperpmin vperpmax vparmin vparmax delv\n0 1 0 2 1\n---\nNaN 0 -2\n3 Inf -Inf\n')
            expected=np.array([[np.nan,0,-2],[3,np.inf,-np.inf]])
            for reader in (loadlinfpccepar,loadlinfpcceperp):
                key='CEpar' if reader==loadlinfpccepar else 'CEperp'
                data=reader(str(file))[key];np.testing.assert_array_equal(data,expected)
                viz=reader(str(file),visualization_nan_to_zero=True)[key]
                self.assertEqual(viz.shape,data.shape);self.assertEqual(viz[0,0],0)
                np.testing.assert_array_equal(viz[1],expected[1])

    def test_cartesian_projection_missing_zero_and_signs(self):
        with tempfile.TemporaryDirectory() as d:
            file=Path(d)/'test.cparcart.specie01.mode01'
            blocks=['NaN 0\n-2 3\n','4 NaN\n0 -5\n','Inf -Inf\n6 0\n']
            file.write_text(self.header+'vxmin vxmax vymin vymax vzmin vzmax delv\n0 1 0 1 0 1 1\n---\n'+'\n'.join(blocks))
            data=loadlinfpccart(str(file));viz=loadlinfpccart(str(file),visualization_nan_to_zero=True)
            for key,block in zip(('CEparvxvy','CEparvxvz','CEparvyvz'),blocks):
                expected=np.loadtxt(io.StringIO(block));np.testing.assert_array_equal(data[key],expected)
                self.assertEqual(data[key].shape,(2,2));self.assertEqual(viz[key].shape,data[key].shape)
                finite=np.isfinite(expected);np.testing.assert_array_equal(viz[key][finite],expected[finite])
                self.assertTrue(np.all(viz[key][np.isnan(expected)]==0))

    def test_cartesian_3d_no_silent_zero_or_clipping(self):
        with tempfile.TemporaryDirectory() as d:
            file=Path(d)/'test.df1.real.specie01.mode01'
            file.write_text(self.header+'bounds\n0 1 0 1 0 1 1\n---\nNaN 0\n-2 3\n1e-120 1e120\nInf -Inf\n')
            data=loadlinfpc3d(str(file));viz=loadlinfpc3d(str(file),visualization_nan_to_zero=True)
            # The storage order is historical (z,y,x); only sample handling changes.
            array=data['fs1'];self.assertEqual(array.shape,(2,2,2))
            self.assertEqual(np.isnan(array).sum(),1);self.assertEqual(np.isinf(array).sum(),2)
            for value in (0,-2,3,1e-120,1e120):self.assertIn(value,array)
            self.assertEqual(np.isnan(viz['fs1']).sum(),0)


class ProductionWorkflowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        subprocess.run(['make','-j2'],cwd=ROOT,check=True,capture_output=True)

    def test_actual_eigen_heating_and_energy_outputs(self):
        for species in ([P,E],[P,CE],[CP,E],[CP,CE]):
            for eigen,heat in ((True,False),(False,True),(True,True)):
                text=source(species,[scan(0,.05,.06,steps=1,res=1,eigen=eigen,heat=heat)])
                text=text.replace("heating_convention='legacy'","heating_convention='corrected'")
                result,files=mixed.run(ROOT/'plume.e',text)
                self.assertEqual(result.returncode,0,result.stdout+result.stderr)
                with tempfile.TemporaryDirectory() as d:
                    p=Path(d);inp=p/'case.in';inp.write_text(text)
                    file=p/'test.mode1';file.write_text(next(v for k,v in files.items() if k.endswith('.mode1') and '_tensor_' not in k))
                    data=read_dispersion_output(file,read_plume_input(inp),1)
                    self.assertEqual(data.heating_available,heat)
                    self.assertEqual(data.params_start,(18+8*len(species) if eigen else 6)+(8*len(species) if heat else 0))
                    for s,sp in enumerate(species):
                        np.testing.assert_allclose(data.data[:,data.params_start+6*s:data.params_start+6*s+6],np.tile(sp,(2,1)),rtol=6e-7)
                        if heat and sp[0]==0:self.assertTrue(np.isnan(data.data[:,data.columns[f'P_species_{s+1}']]).all())
                    if heat:
                        energy=p/'test.energy';energy.write_text(next(v for k,v in files.items() if k.endswith('.energy')))
                        self.assertEqual(len(load_plume_energy(energy)['record']),2)

    def test_actual_fpc_files_preserve_samples_and_injected_nan(self):
        from test_normalization import input_text,run_case
        sp=[(1,1,1.3,1,1,0),(.5,1836,1,-1,1,.02)]
        for option in (6,7):
            text=input_text('rho_ref',.4,1.3,.05,.06,sp).replace('option=-1',f'option={option}')
            text+='&fpc\nvperpmin=0\nvperpmax=.5\nvparmin=-.5\nvparmax=.5\nvxmin=-.5\nvxmax=.5\nvymin=-.5\nvymax=.5\nvzmin=-.5\nvzmax=.5\ndelv=.5\n/\n&guess_1\ng_om=.08\ng_gam=-.001\n/\n'
            result,files=run_case(ROOT/'plume.e',text)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            suffix='.cpar.specie01.mode01' if option==6 else '.cparcart.specie01.mode01'
            name=next(k for k in files if k.endswith(suffix));contents=files[name]
            with tempfile.TemporaryDirectory() as d:
                file=Path(d)/name;file.write_text(contents)
                reader=loadlinfpccepar if option==6 else loadlinfpccart
                key='CEpar' if option==6 else 'CEparvxvy'
                actual=reader(str(file))[key]
                start=5;lines=contents.splitlines()
                self.assertGreater(len(lines[start].split()),1)
                width=actual.shape[1]
                expected=np.array([[float(x) for x in line.split()[:width]] for line in lines[start:start+actual.shape[0]]])
                np.testing.assert_array_equal(actual,expected)
                values=lines[start].split();values[0]='NaN';values[1]='0.0'
                lines[start]=' '.join(values);file.write_text('\n'.join(lines)+'\n')
                loaded=reader(str(file))[key]
                self.assertEqual(loaded.shape,actual.shape)
                self.assertTrue(np.isnan(loaded[0,0]));self.assertEqual(loaded[0,1],0)
                np.testing.assert_array_equal(loaded[1:],actual[1:])

    def test_log_map_actual_historical_coordinates(self):
        for real_log,imag_log in ((True,False),(False,True),(True,True)):
            for sign in (-1,1):
                g0=sign*1e-4;g1=sign*1e-2
                text=source([CP,CE],[],option=0)
                text+=f'&maps loggridw={str(real_log).replace("True",".true.").replace("False",".false.")}, loggridg={str(imag_log).replace("True",".true.").replace("False",".false.")}, omi=.1, omf=1, gami={g0}, gamf={g1}, nr=2, ni=2 /\n'
                result,files=mixed.run(ROOT/'plume.e',text)
                self.assertEqual(result.returncode,0,result.stdout+result.stderr)
                self.assertEqual(result.stdout.count('WARNING: historical logarithmic frequency maps'),1)
                actual=np.loadtxt(io.StringIO(next(v for k,v in files.items() if k.endswith('.map'))))
                wr=10**(.1+np.arange(3)*.5) if real_log else np.linspace(.1,1,3)
                wi=-10**(abs(g0)+np.arange(3)) if imag_log else np.linspace(g0,g1,3)
                np.testing.assert_allclose(np.unique(actual[:,0]),np.sort(wr),rtol=6e-7)
                np.testing.assert_allclose(np.unique(actual[:,1]),np.sort(wi),rtol=6e-7)

    def test_ignored_map_requests_warn_without_schema_change(self):
        text=source([P,E],[scan(0,.05,.06,steps=1,res=1,eigen=True,heat=True,tensor=True)],option=4)
        text+='&maps loggridw=.false.,loggridg=.false.,omi=.04,omf=.08,gami=-.001,gamf=.001,nr=2,ni=2 /\n'
        result,files=mixed.run(ROOT/'plume.e',text)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        self.assertEqual(result.stdout.count('eigen/tensor requests are ignored'),1)
        self.assertEqual(result.stdout.count('legacy heating request is ignored'),1)
        for name,contents in files.items():
            if name.endswith('.map'):
                self.assertEqual(np.loadtxt(io.StringIO(contents)).shape[1],5)
        self.assertFalse(any(k.endswith('.energy') for k in files))

    def test_all_kinetic_parallel_warns_but_mixed_parallel_does_not(self):
        for species in ([P,E],[P,CE],[CP,E],[CP,CE]):
            text=source(species,[],option=0,kperp=0)
            text+='&maps loggridw=.false.,loggridg=.false.,omi=.04,omf=.08,gami=.001,gamf=.002,nr=1,ni=1 /\n'
            result,_=mixed.run(ROOT/'plume.e',text)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            self.assertEqual(result.stdout.count('WARNING: exact-parallel all-kinetic'),int(all(s[0]>0 for s in species)))


class CheckedWorkflowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from faddeeva_build import compile_c,c_objects
        from plasma_dispersion_audit import MODULES
        cls.temp=tempfile.TemporaryDirectory();cls.addClassCleanup(cls.temp.cleanup)
        directory=Path(cls.temp.name);compile_c(directory)
        flags=['-O0','-g','-fdefault-real-8','-fcheck=all','-ffpe-trap=invalid,zero,overflow']
        base=['gfortran',*flags,'-I',str(directory),'-J',str(directory)]
        out=subprocess.run(base+['-c',*[str(ROOT/'src'/f'{m}.f90') for m in (*MODULES,'fpc')]],
                           cwd=directory,text=True,capture_output=True)
        if out.returncode:raise RuntimeError(out.stderr)
        objects=[str(directory/f'{m}.o') for m in MODULES]+c_objects(directory)
        cls.input_probe=directory/'input_probe';cls.plume=directory/'plume'
        for executable,path,extra in [(cls.input_probe,ROOT/'tests/workflow_input_probe.f90',[]),
                                       (cls.plume,ROOT/'src/plume.f90',[str(directory/'fpc.o')])]:
            out=subprocess.run(base+[str(path),*objects,*extra,'-o',str(executable)],cwd=directory,text=True,capture_output=True)
            if out.returncode:raise RuntimeError(out.stderr)

    def test_checked_log_map_matches_production(self):
        text=source([CP,CE],[],option=0)
        text+='&maps loggridw=.true.,loggridg=.true.,omi=.1,omf=1,gami=-.0001,gamf=-.01,nr=2,ni=2 /\n'
        normal,files=mixed.run(ROOT/'plume.e',text);checked,checked_files=mixed.run(self.plume,text)
        self.assertEqual(normal.returncode,0,normal.stderr);self.assertEqual(checked.returncode,0,checked.stderr)
        name=next(k for k in files if k.endswith('.map'))
        a=np.loadtxt(io.StringIO(files[name]));b=np.loadtxt(io.StringIO(checked_files[name]))
        np.testing.assert_array_equal(a[:,:2],b[:,:2])
        np.testing.assert_allclose(a[:,2:],b[:,2:],rtol=1e-12,atol=0)

    def test_checked_input_warnings_do_not_evaluate_unsupported_kernel(self):
        text=source([P,E],[scan(0,0,.01,eigen=True,heat=True,tensor=True)],option=4,kperp=0)
        text+='&maps loggridw=.true.,loggridg=.true.,omi=.1,omf=1,gami=-.0001,gamf=-.01,nr=2,ni=2 /\n'
        result,_=mixed.run(self.input_probe,text)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        for warning in ('WARNING: exact-parallel all-kinetic','eigen/tensor requests are ignored',
                        'legacy heating request is ignored','WARNING: historical logarithmic frequency maps'):
            self.assertEqual(result.stdout.count(warning),1,result.stdout)
        self.assertIn('INPUT_OK',result.stdout)
