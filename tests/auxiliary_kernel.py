"""Build-time access to the unchanged collisionless kinetic expressions.

This is an isolated experiment, not a second maintained kinetic implementation.
Extract current disp0, bind normalization anisotropy explicitly, and disable its
legacy cache (which implicitly couples normalization to species 1). Production
source, objects and callers are untouched. Fail loudly if source structure moves.
"""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]


def experimental_kernel():
    source = (ROOT/'src/disprels.f90').read_text()
    start = source.index('   complex function disp0(om)')
    end = source.index('   end function disp0', start)+len('   end function disp0')
    body = source[start:end]
    assert body.count('alphp = spec(1)%alph_s') == 1
    cache_start = body.index('     if ((kperp .eq. kperp_last )')
    cache_end = body.index('     !Arguments of Bessel functions', cache_start)
    # reuse_bessel is already false; remove only the legacy cache lookup.
    body = body[:cache_start]+body[cache_end:]
    body = body.replace('disp0(om)', 'probe_disp0(om, normalization_alpha)', 1)
    body = body.replace('     complex :: om\n',
                        '     complex :: om\n     real, intent(in) :: normalization_alpha\n', 1)
    body = body.replace('alphp = spec(1)%alph_s', 'alphp = normalization_alpha')
    body = body.replace('     disp0=', '     probe_disp0=', 1)
    body = body.replace('end function disp0', 'end function probe_disp0', 1)
    count = re.search(r'integer, parameter :: nbrack=(\d+)', source).group(1)
    return ('! Generated for tests only; never modify or commit the generated body.\n'
            'module auxiliary_kinetic_kernel\n'
            '  use disprels, only: bessel, zet_in\n'
            '  implicit none\n'
            f'  integer, parameter :: nbrack={count}\n'
            'contains\n'+body+'\nend module auxiliary_kinetic_kernel\n')
