"""Shared C99 build support for standalone and checked Fortran probes."""
from pathlib import Path
import subprocess
ROOT=Path(__file__).resolve().parents[1]
C_MODULES=('faddeeva','faddeeva_bridge')
def c_objects(directory):
    return [str(Path(directory)/f'{name}.o') for name in C_MODULES]
def compile_c(directory,pic=False):
    directory=Path(directory)
    for name,source in [('faddeeva',ROOT/'third_party/faddeeva/Faddeeva.cc'),
                        ('faddeeva_bridge',ROOT/'src/faddeeva_bridge.c')]:
        subprocess.run(['cc','-O3','-std=c99','-fno-fast-math',*(['-fPIC'] if pic else []),
                        '-I',str(ROOT/'third_party/faddeeva'),'-x','c','-c',str(source),
                        '-o',str(directory/f'{name}.o')],check=True,capture_output=True)
    return c_objects(directory)
