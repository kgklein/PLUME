"""Small scalar namelist reader for PLUME's existing Python wrapper.

Not a general Fortran parser: unsupported groups/arrays/expressions fail instead
of silently dropping controls. Field names are case insensitive; output uses
canonical PLUME spelling. No new runtime dependency is required.
"""
import re

_FIELDS = {
    'params': {
        **dict.fromkeys(['heating_convention','z_backend','bessel_mode','normalization','dataName','outputName'], str),
        **dict.fromkeys(['bessel_nmax','bessel_nfixed','nspec','nscan','option','nroot_max','collision_type'], int),
        **dict.fromkeys(['bessel_diagnostics','use_map','low_n','writeOut'], bool),
        **dict.fromkeys(['bessel_tol','betap','kperp','kpar','vtp','Kn'], float)},
    'species':dict.fromkeys(['tauS','muS','alphS','Qs','Ds','vvS'], float),
    'maps': {**dict.fromkeys(['loggridw','loggridg','positive_roots'], bool),
             **dict.fromkeys(['omi','omf','gami','gamf'], float), 'nr':int,'ni':int},
    'scan_input': {**dict.fromkeys(['scan_type','scan_style','ns','nres'], int),
                   **dict.fromkeys(['swi','swf'], float),
                   **dict.fromkeys(['swlog','heating','eigen','tensor'], bool)},
    'guess':dict.fromkeys(['g_om','g_gam'], float),
    'fpc':dict.fromkeys(['vperpmin','vperpmax','vparmin','vparmax','delv','vxmin','vxmax','vymin','vymax','vzmin','vzmax','elecdircontribution'], float),
}
_TOKEN = re.compile(r'''\s+|![^\n]*|'(?:[^']|'')*'|"(?:[^"]|"")*"|&[A-Za-z_]\w*|\.(?:true|false|t|f)\.|[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[ed][+-]?\d+)?|[A-Za-z_]\w*|[=,/]''', re.I)


def _schema(group):
    base=re.sub(r'_\d+$','',group.lower())
    if base not in _FIELDS:
        raise ValueError(f'Unsupported namelist &{group}; the Python wrapper supports params, species_N, maps, scan_input_N, guess_N and fpc')
    return _FIELDS[base]


def _canonical(group,key):
    schema=_schema(group);names={name.lower():name for name in schema}
    if key.lower() not in names:
        raise ValueError(f'Unrecognized input control {key!r} in &{group}; refusing to discard it')
    name=names[key.lower()]
    return name,schema[name]


def _value(raw,kind,group,key):
    try:
        if kind is str:
            if raw.startswith(("'",'"')):return raw[1:-1].replace(raw[0]*2,raw[0])
            return raw
        if kind is bool:
            logical=raw.lower().strip('.')
            if logical not in ('true','false','t','f'):raise ValueError('expected a Fortran logical')
            return '.true.' if logical.startswith('t') else '.false.'
        if kind is int:return int(raw)
        return float(raw.lower().replace('d','e'))
    except ValueError as exc:
        raise ValueError(f'Invalid scalar {raw!r} for {key} in &{group}') from exc


def read_scalar_namelists(text):
    tokens=[];position=0
    while position<len(text):
        match=_TOKEN.match(text,position)
        if not match:raise ValueError(f'Unsupported namelist syntax near {text[position:position+35]!r}')
        token=match.group();position=match.end()
        if not token.isspace() and not token.startswith('!'):tokens.append(token)
    groups=[];i=0
    while i<len(tokens):
        if not tokens[i].startswith('&'):raise ValueError(f'Expected namelist header, found {tokens[i]!r}')
        group=tokens[i][1:].lower();_schema(group);i+=1;values={}
        while i<len(tokens) and tokens[i].lower() not in ('/','&end'):
            if tokens[i]==',':i+=1;continue
            key,kind=_canonical(group,tokens[i]);i+=1
            if i>=len(tokens) or tokens[i]!='=':raise ValueError(f'Expected = after {key} in &{group}')
            i+=1
            if i>=len(tokens) or tokens[i] in ('/',',') or tokens[i].startswith('&'):raise ValueError(f'Missing scalar for {key} in &{group}')
            values[key]=_value(tokens[i],kind,group,key);i+=1
        if i>=len(tokens):raise ValueError(f'Unterminated namelist &{group}')
        i+=1;groups.append((group,values))
    return groups


def validate_controls(group,values,derived=()):
    """Check names before a writer can delete/replace an existing input file."""
    for key in values:
        if key not in derived:_canonical(group,key)


def scalar_text(group,key,value):
    """Serialize scalar controls without losing logicals or quoted path spaces."""
    _,kind=_canonical(group,key)
    if kind is str:return "'"+str(value).replace("'","''")+"'"
    if kind is bool:
        if isinstance(value,bool):return '.true.' if value else '.false.'
        raw=str(value).lower()
        if raw not in ('.true.','.false.'):raise ValueError(f'{key} requires a logical value')
        return raw
    return str(value)
