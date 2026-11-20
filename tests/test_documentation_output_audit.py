"""Frozen Stage 5A read-only characterization, not desired future bug behavior.

When fixing the named issues, update/retire the defect-characterization tests;
these are separate from the solver's historical numerical golden values.
"""
import contextlib
import io
from pathlib import Path
import tempfile
import unittest
import numpy as np
from documentation_output_audit import ROOT,namelists,parse_input
from documentation_workflow_probes import markers
import sys
sys.path.insert(0,str(ROOT))
from linfpclib.linfpc import load_plume_energy,loadlinfpccepar,plume_input

class DocumentationAuditTests(unittest.TestCase):
    def test_all_namelist_groups_are_inventoried(self):
        groups=namelists();self.assertEqual(len(groups),14)
        self.assertEqual(len(set(p.lower() for g in groups for p in g['parameters'])),89)
        self.assertIn('heating_convention',groups[0]['parameters'])
    def test_example_inventory_parser_preserves_signed_parameters(self):
        groups=parse_input("&params normalization='d_ref' betap=0 /\n&species_1 tauS=0 Qs=-1 vvS=-0.1 /\n")
        self.assertEqual(groups[1]['values'],dict(taus=0,qs=-1,vvs=-.1))
    def test_current_and_debug_reader_unique_markers(self):
        with tempfile.TemporaryDirectory() as d:
            with contextlib.redirect_stdout(io.StringIO()):data=markers(ROOT,Path(d))
        records=[r for r in data['markers'] if 'passed' in r]
        self.assertEqual(len(records),33);self.assertTrue(all(r['passed'] for r in records))
    def test_low_n_false_reader_is_currently_unsupported(self):
        with tempfile.TemporaryDirectory() as d:
            with contextlib.redirect_stdout(io.StringIO()):data=markers(ROOT,Path(d))
        records=[r for r in data['markers'] if r.get('case')=='low_n=false']
        self.assertEqual(len(records),4);self.assertTrue(all('IndexError' in r['outcome'] for r in records))
    def test_writeout_is_independent_of_use_map(self):
        obj=plume_input('audit');obj.set_params(.4,.05,.06,.01,2,0,1,1,False,True,normalization='d_ref')
        self.assertEqual(obj.params['writeOut'],'.true.')
    def test_fpc_reader_preserves_nan(self):
        with tempfile.TemporaryDirectory() as d:
            file=Path(d)/'test.cpar.specie01.mode01'
            file.write_text('tau bi kpar kperp vti mu omega.r omega.i\n1 1 .1 .1 .01 1 .1 -.001\nvperpmin vperpmax vparmin vparmax delv\n0 1 -1 1 1\n---\nNaN 2 3\n4 5 6\n')
            actual=loadlinfpccepar(str(file))
            self.assertTrue(np.isnan(actual['CEpar'][0,0]))
    def test_energy_reader_preserves_failed_rows_and_quoted_paths(self):
        with tempfile.TemporaryDirectory() as d:
            file=Path(d)/'test.energy';file.write_text('# schema=plume-energy-v1\n# columns: record workflow scan outer point mode root_valid field_status energy_status applicability response_code H_total source_file\n1 2 0 0 0 1 0 -1 8 0 0 NaN "data/a path/mode.mode1"\n')
            actual=load_plume_energy(file);self.assertEqual(actual['root_valid'][0],0);self.assertTrue(np.isnan(actual['H_total'][0]))
            self.assertEqual(actual['source_file'][0],'data/a path/mode.mode1')
    def test_legacy_debug_reader_has_no_general_species_schema(self):
        import inspect
        from linfpclib.linfpc import loadeigen
        self.assertIn("eigendict['Pe'] = complex(0,0)",inspect.getsource(loadeigen))
