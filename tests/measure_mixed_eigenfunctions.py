"""Reproduce Stage 4B.1 numerical comparisons without rewriting golden data."""
import argparse
import contextlib
import io
import json
import platform
import subprocess
import unittest
from pathlib import Path
from test_mixed_eigenfunctions import MixedEigenfunctionTests


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    stream=io.StringIO()
    with contextlib.redirect_stdout(stream):
        suite=unittest.defaultTestLoader.loadTestsFromTestCase(MixedEigenfunctionTests)
        result=unittest.TextTestRunner(stream=stream,verbosity=2).run(suite)
    if not result.wasSuccessful():raise RuntimeError(stream.getvalue())
    data=dict(compiler=subprocess.check_output(['gfortran','--version'],text=True).splitlines()[0],
              platform=platform.platform(),tests=result.testsRun,measurements=MixedEigenfunctionTests.measures)
    args.output.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    print('Wrote',args.output)


if __name__=='__main__':main()
