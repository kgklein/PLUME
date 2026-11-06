"""Read-only manual coverage/link/math checks; never stages or deploys a site.

python3 tests/check_user_manual.py --site /tmp/<preview>/docs_out --output /tmp/manual-checks.json
External literature URLs are not certified by this local check. Repository URLs
are checked against this checkout, including downloadable input sources.
"""
import argparse
from html.parser import HTMLParser
import json
from pathlib import Path
import re
from urllib.parse import unquote, urlsplit

ROOT=Path(__file__).resolve().parents[1]
PAGES=['README.md','INSTALL.md','tutorial.md','input.md','normalization.md',
       'models.md','modes.md','eigenfunctions.md','heating.md','output.md',
       'examples.md','numerics.md','limitations.md','reproducibility.md',
       'development.md','README-JETPLUME.md']


def repository_target(url):
    path=unquote(url.path)
    if url.netloc=='raw.githubusercontent.com' and path.startswith('/kgklein/PLUME/main/'):
        return path.removeprefix('/kgklein/PLUME/main/')
    if url.netloc=='github.com' and path.startswith('/kgklein/PLUME/blob/main/'):
        return path.removeprefix('/kgklein/PLUME/blob/main/')
    return None


class Document(HTMLParser):
    def __init__(self,text):
        super().__init__(convert_charrefs=True)
        self.ids=set();self.links=[];self.math=[];self._math=None;self.tables=0
        self._table_depth=0;self.table_emphasis=0
        self.feed(text)

    def handle_starttag(self,tag,attrs):
        attributes=dict(attrs)
        if tag=='table':self.tables+=1;self._table_depth+=1
        if tag=='em' and self._table_depth:self.table_emphasis+=1
        if 'id' in attributes:self.ids.add(attributes['id'])
        if 'href' in attributes:self.links.append(attributes['href'])
        if tag in ('img','script') and 'src' in attributes:self.links.append(attributes['src'])
        if tag=='script' and attributes.get('type','').startswith('math/tex'):
            self._math=[]

    def handle_data(self,data):
        if self._math is not None:self._math.append(data)

    def handle_endtag(self,tag):
        if tag=='table':self._table_depth-=1
        if tag=='script' and self._math is not None:
            self.math.append(''.join(self._math));self._math=None


def source_checks():
    from documentation_output_audit import namelists
    actual={name.lower() for group in namelists() for name in group['parameters']}
    text=(ROOT/'input.md').read_text()
    documented=re.findall(r'<a id="param-([\w_]+)"></a>',text)
    errors=[];downloads=set();links=0
    if set(documented)!=actual or len(documented)!=len(actual):
        errors.append(dict(input_missing=sorted(actual-set(documented)),
                           input_extra=sorted(set(documented)-actual),records=len(documented)))
    for name in PAGES:
        source=(ROOT/name).read_text()
        if source.count(r'\[')!=source.count(r'\]') or source.count(r'\(')!=source.count(r'\)'):
            errors.append(dict(page=name,error='Unbalanced math delimiters'))
        for target in re.findall(r'!?\[[^\]]*\]\(([^\s)]+)\)',source):
            url=urlsplit(target);path=unquote(url.path);links+=1
            if url.netloc=='raw.githubusercontent.com' and path.startswith('/kgklein/PLUME/main/'):
                relative=path.removeprefix('/kgklein/PLUME/main/')
            elif url.netloc=='github.com' and path.startswith('/kgklein/PLUME/blob/main/'):
                relative=path.removeprefix('/kgklein/PLUME/blob/main/')
            elif not url.scheme and not url.netloc:
                relative=path
            else:continue
            if not relative:continue
            resolved=ROOT/relative
            if not (resolved.is_file() or (resolved.is_dir() and (resolved/'README.md').is_file())):
                errors.append(dict(page=name,missing_source=relative))
            if relative.startswith('inputs/'):
                downloads.add(relative)
    return dict(parameter_names=len(actual),parameter_records=len(documented),manual_pages=len(PAGES),
                source_links=links,downloadable_inputs=sorted(downloads),errors=errors)


def site_checks(site):
    site=site.resolve();cache={};broken=[];missing=[];math={};tables={};api_links=set()
    def read(path):
        if path not in cache:cache[path]=Document(path.read_text())
        return cache[path]
    for name in ['ford_project.md',*PAGES]:
        page=site/'index.html' if name=='ford_project.md' else site/'page'/('index.html' if name=='README.md' else Path(name).with_suffix('.html'))
        if not page.is_file():missing.append(name);continue
        parsed=read(page);math[name]=len(parsed.math);tables[name]=parsed.tables
        source=(ROOT/name).read_text()
        if name=='input.md':
            expected_ids=set(re.findall(r'param-[\w_]+',source))
            if not expected_ids<=parsed.ids:
                broken.append(dict(page=name,error='Input anchors lost in rendering',missing=sorted(expected_ids-parsed.ids)))
            if parsed.table_emphasis:
                broken.append(dict(page=name,error='Unintended emphasis in input formula cells; escape multiplication asterisks'))
        intended=len(re.findall(r'^\|.*\|\n\|[-: |]+\|',source,re.M))
        if intended!=parsed.tables:
            broken.append(dict(page=name,error='Markdown table not rendered',expected=intended,actual=parsed.tables))
        for target in parsed.links:
            url=urlsplit(target)
            if url.scheme or url.netloc:
                repository=repository_target(url)
                if repository is not None and not (ROOT/repository).is_file():
                    broken.append(dict(page=name,target=target,error='missing repository target after staging'))
                continue
            destination=(page.parent/unquote(url.path)).resolve() if url.path else page
            if destination.is_dir():destination=destination/'index.html'
            if not destination.is_file():
                broken.append(dict(page=name,target=target,error='missing file'));continue
            if destination.suffix=='.html':
                if url.fragment and unquote(url.fragment) not in read(destination).ids:
                    broken.append(dict(page=name,target=target,error='missing anchor'))
                if destination.parent.name in ('proc','module','lists'):
                    api_links.add(str(destination.relative_to(site)))
        if not parsed.math and name in ('normalization.md','eigenfunctions.md','heating.md','numerics.md'):
            broken.append(dict(page=name,error='Scientific page has no TeX math scripts'))
    return dict(site=str(site),missing_pages=missing,broken_links=broken,math_scripts=math,
                tables=tables,valid_api_destinations=sorted(api_links))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--site',type=Path);parser.add_argument('--output',type=Path)
    args=parser.parse_args();result={'source':source_checks()}
    if args.site:result['site']=site_checks(args.site)
    text=json.dumps(result,indent=2)
    if args.output:args.output.write_text(text+'\n')
    print(text)
    # Missing staging pages are a real reported failure, not a false pass.
    raise SystemExit(bool(result['source']['errors'] or
                          (args.site and (result['site']['missing_pages'] or result['site']['broken_links']))))
