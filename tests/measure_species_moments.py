"""Reproduce Stage 4B.2 measurements; no regression expectations are rewritten."""
import argparse
import contextlib
import io
import json
import platform
import subprocess
import unittest
from pathlib import Path
import numpy as np
import scipy
from test_species_moments import SpeciesMomentTests


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    log=io.StringIO()
    with contextlib.redirect_stdout(log):
        suite=unittest.defaultTestLoader.loadTestsFromTestCase(SpeciesMomentTests)
        result=unittest.TextTestRunner(stream=log,verbosity=2).run(suite)
    if not result.wasSuccessful():raise RuntimeError(log.getvalue())
    data=dict(compiler=subprocess.check_output(['gfortran','--version'],text=True).splitlines()[0],
              platform=platform.platform(),numpy=np.__version__,scipy=scipy.__version__,tests=result.testsRun,
              measurements=SpeciesMomentTests.measures)
    args.output.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n');print('Wrote',args.output)


if __name__=='__main__':main()
