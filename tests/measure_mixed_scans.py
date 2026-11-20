"""Reproduce Stage 4A scan comparisons; writes new results, never golden data."""
import argparse
import contextlib
import io
import json
import math
import platform
import subprocess
import unittest
from pathlib import Path
from test_mixed_scans import MixedScanTests


def finite_json(value):
    if isinstance(value,dict):return {k:finite_json(v) for k,v in value.items()}
    if isinstance(value,list):return [finite_json(v) for v in value]
    if isinstance(value,float) and not math.isfinite(value):return None
    return value


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    stream=io.StringIO()
    with contextlib.redirect_stdout(stream):
        suite=unittest.defaultTestLoader.loadTestsFromTestCase(MixedScanTests)
        result=unittest.TextTestRunner(stream=stream,verbosity=2).run(suite)
    if not result.wasSuccessful():raise RuntimeError(stream.getvalue())
    data={'compiler':subprocess.check_output(['gfortran','--version'],text=True).splitlines()[0],
          'platform':platform.platform(),'tests':result.testsRun,'measurements':MixedScanTests.measurements}
    args.output.write_text(json.dumps(finite_json(data),indent=2,allow_nan=False)+'\n')
    print('Wrote',args.output)


if __name__=='__main__':main()
