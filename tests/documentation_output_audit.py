"""Stage 5A read-only inventory/crawler and isolated characterization helpers.

Writes only the explicit --output/--cache destinations. No existing source,
input, output, docs, reader or test is edited. Use a temporary checkout for
builds/example execution. Optional audit dependencies: requests, bs4, NumPy.
"""
import argparse
from collections import deque
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import re
import subprocess
from urllib.parse import urljoin,urlsplit,urldefrag,unquote

ROOT=Path(__file__).resolve().parents[1]
BASE='https://kgklein.github.io/PLUME/'


def documentation_sources():
    """Current project overview, FORD configuration and manual source pages."""
    return [ROOT/'README.md', ROOT/'ford_project.md', *sorted((ROOT/'docs').glob('*.md'))]


def crawl(output,cache):
    import requests
    from bs4 import BeautifulSoup
    cache.mkdir(parents=True,exist_ok=True)
    expected={urljoin(BASE,'page/'+p.stem+'.html') for p in (ROOT/'docs').glob('*.md')}
    expected|={urljoin(BASE,'page/'+p.stem+'.html') for p in (ROOT/'tests').glob('*REPORT.md')}
    expected|={urljoin(BASE,'module/'+name.lower()+'.html') for p in (ROOT/'src').glob('*.f90') for name in re.findall(r'(?im)^\s*module\s+(?!procedure\b)(\w+)',p.read_text())}
    pending={BASE}|expected;seen=set();pages={};resources={};total=0
    def fetch(url):
        try:
            r=requests.get(url,timeout=25)
            return url,r.status_code,r.url,r.headers,r.text,None
        except requests.RequestException as e:return url,None,url,{},'',str(e)
    with ThreadPoolExecutor(max_workers=6) as pool:
        while pending:
            urls=sorted(pending-seen);pending=set()
            if not urls:break
            for url,code,final,headers,html,error in pool.map(fetch,urls):
                seen.add(url);total+=1
                item=dict(expected_current_page=url in expected,url=url,status=code,final_url=final,error=error,last_modified=headers.get('Last-Modified'),etag=headers.get('ETag'))
                pages[url]=item
                if code!=200:continue
                digest=hashlib.sha256(html.encode()).hexdigest();name=hashlib.sha256(url.encode()).hexdigest()[:20]
                (cache/(name+'.html')).write_text(html)
                soup=BeautifulSoup(html,'html.parser');body=soup.find('body') or soup
                for tag in body.find_all(['script','style']):tag.decompose()
                text=body.get_text('\n',strip=True);(cache/(name+'.txt')).write_text(text)
                ids=[tag['id'] for tag in soup.find_all(id=True)]
                links=[]
                for tag in soup.find_all(['a','img','script','link']):
                    href=tag.get('href') or tag.get('src')
                    if not href or href.startswith(('mailto:','javascript:','data:')):continue
                    target=urljoin(url,href);clean,fragment=urldefrag(target)
                    links.append(dict(text=tag.get_text(' ',strip=True)[:120],url=target,tag=tag.name,fragment=fragment))
                    if clean.startswith(BASE):
                        if clean.endswith(('.html','/')):
                            if clean not in seen:pending.add(clean)
                        else:resources.setdefault(clean,dict(referenced_by=[]))['referenced_by'].append(url)
                item.update(title=soup.title.get_text(' ',strip=True) if soup.title else None,
                    headings=[h.get_text(' ',strip=True) for h in soup.find_all(re.compile('^h[1-6]$'))],
                    cache=name,sha256=digest,words=len(text.split()),ids=ids,links=links,
                    math_delimiters=dict(display=html.count('\\['),inline=html.count('\\('),dollars=html.count('$$')),
                    mathjax='mathjax' in html.lower(),code_blocks=len(soup.find_all('pre')),
                    capabilities={key:key.lower() in text.lower() for key in ['d_ref','rho_ref','tauS','heating_convention','Faddeeva','bessel_mode','wave energy','energy-v1','cold-fluid','relativistic','parallelised']})
                if total>1600:raise RuntimeError('Unexpected crawl size: refuse unbounded crawl')
        # Same-site non-HTML resources, including erroneously deployed .md links.
        for url,code,final,headers,content,error in pool.map(fetch,sorted(resources)):
            resources[url].update(status=code,final_url=final,error=error)
    broken=[]
    for url,item in pages.items():
        for link in item.get('links',[]):
            target,fragment=urldefrag(link['url'])
            if not target.startswith(BASE):continue
            record=pages.get(target) or resources.get(target)
            reason=None
            if record and record.get('status')!=200:reason='HTTP '+str(record.get('status'))
            elif fragment and target in pages and unquote(fragment) not in pages[target].get('ids',[]):reason='missing anchor '+fragment
            if reason:broken.append(dict(source=url,target=link['url'],text=link['text'],reason=reason))
    result=dict(retrieved_utc=datetime.now(timezone.utc).isoformat(),base=BASE,pages=list(pages.values()),resources=resources,broken_internal_links=broken,
        limitations=['Static HTML and math configuration audited; no browser/JavaScript rendering engine used.','External links are inventoried but not exhaustively fetched.'])
    output.write_text(json.dumps(result,indent=2)+'\n');print('Pages',len(pages),'resources',len(resources),'broken links',len(broken))


def namelists():
    source=(ROOT/'src/functions.f90').read_text();groups=[]
    logical=[]
    for number,line in enumerate(source.splitlines(),1):
        code=line.split('!')[0].strip()
        if logical and logical[-1]['continuation']:
            logical[-1]['text']+=' '+code.lstrip('&').rstrip('&');logical[-1]['continuation']=code.endswith('&')
        elif code:logical.append(dict(line=number,text=code.rstrip('&'),continuation=code.endswith('&')))
    for line in logical:
        m=re.match(r'(?i)namelist\s*/(\w+)/\s*(.*)',line['text'])
        if m:groups.append(dict(group=m[1],parameters=[p.strip() for p in m[2].split(',') if p.strip()],source='src/functions.f90',line=line['line']))
    return groups


def parse_input(text):
    # Inventory parser, not a replacement for a Fortran namelist reader.
    text=re.sub(r'!.*','',text)
    groups=[]
    for m in re.finditer(r'(?is)&(\w+)\b(.*?)(?<!\w)/(?:\s|$)',text):
        values={}
        for a in re.finditer(r"(?im)\b([a-z]\w*)\s*=\s*('[^']*'|\"[^\"]*\"|[^,\s/]+)",m[2]):
            raw=a[2];v=raw.strip("'\"")
            if v.lower() in ['.true.','.false.']:v=v.lower()=='.true.'
            else:
                try:v=float(v.lower().replace('d','e'));v=int(v) if v.is_integer() else v
                except ValueError:pass
            values[a[1].lower()]=v
        groups.append(dict(group=m[1].lower(),values=values))
    return groups


def inventory(output):
    tracked=set(subprocess.check_output(['git','ls-files'],cwd=ROOT,text=True).splitlines())
    defaults=dict(normalization='rho_ref',z_backend='faddeeva',heating_convention='corrected',bessel_mode='adaptive',bessel_tol=1e-10,bessel_nmax=512,bessel_nfixed=128,bessel_diagnostics=False,nscan=0,collision_type=0,kn=1e30,low_n=True,positive_roots=False,nr=128,ni=128)
    docs={str(p.relative_to(ROOT)):p.read_text() for p in documentation_sources()}
    examples=[]
    for p in sorted((ROOT/'inputs').rglob('*.in')):
        if str(p.relative_to(ROOT)) not in tracked and p.parent != ROOT/'inputs/example':continue
        groups=parse_input(p.read_text(errors='replace'));params=next((g['values'] for g in groups if g['group']=='params'),{})
        species=[g['values'] for g in groups if g['group'].startswith('species_')]
        scans=[g['values'] for g in groups if g['group'].startswith('scan_input_')]
        effective={**defaults,**params};warm=sum(s.get('taus',-1)>0 for s in species);cold=sum(s.get('taus',-1)==0 for s in species)
        examples.append(dict(path=str(p.relative_to(ROOT)),tracked=str(p.relative_to(ROOT)) in tracked,params=params,effective=effective,species=species,scans=scans,
            composition='mixed' if warm and cold else 'cold' if cold else 'kinetic' if warm else 'unknown',reference='cold' if species and species[0].get('taus')==0 else 'warm/unknown',
            documented_in=[str(k) for k,v in docs.items() if p.name in v],execution='not run',scientific_interpretation='not inferred from return code'))
    namelist_groups=namelists();allnames=sorted(set(p.lower() for g in namelist_groups for p in g['parameters']))
    params=[]
    for name in allnames:
        locations=[dict(group=g['group'],source=g['source'],line=g['line']) for g in namelist_groups if name in [p.lower() for p in g['parameters']]]
        params.append(dict(name=name,groups=locations,default=defaults.get(name,'REQUIRED/uninitialized: no documented safe default established'),
            documentation=[k for k,v in docs.items() if re.search(r'(?i)\b'+re.escape(name)+r'\b',v)],
            examples=[e['path'] for e in examples if name in e['params'] or any(name in g for g in e['species']+e['scans'])][:5]))
    plotting=[]
    for folder in ['linfpclib','plotter','SWIFT','contours']:
        for p in sorted((ROOT/folder).rglob('*')):
            if p.is_file() and p.suffix in ['.py','.plt','.gp','.sh'] and (str(p.relative_to(ROOT)) in tracked or p.parent in [ROOT/'plotter/example',ROOT/'linfpclib']):
                text=p.read_text(errors='replace');plotting.append(dict(path=str(p.relative_to(ROOT)),tracked=str(p.relative_to(ROOT)) in tracked,
                    functions=re.findall(r'(?m)^def (\w+)\(',text),lines=len(text.splitlines()),
                    schemas=sorted(set(re.findall(r'\.(?:mode[1-9]?|energy|roots|map|fpc|eig|sweep)\b',text))),
                    hardcoded_evidence=[dict(line=i,text=l.strip()[:220]) for i,l in enumerate(text.splitlines(),1) if re.search(r'\b(?:nspec\s*=\s*[234]|use_ps_split|load_plume_|using\s+[0-9]|\$[1-9]|p[123](?:ld|ttd|cd)|rho_ref|betap|heating_convention)',l)][:25]))
    output.write_text(json.dumps(dict(namelist_groups=namelist_groups,input_parameters=params,examples=examples,python_and_plotting=plotting),indent=2)+'\n')
    print(len(params),'distinct input names;',len(namelist_groups),'namelist definitions;',len(examples),'input files;',len(plotting),'interfaces/plot scripts')


def main():
    p=argparse.ArgumentParser();p.add_argument('action',choices=['site','inventory']);p.add_argument('--output',type=Path,required=True);p.add_argument('--cache',type=Path)
    a=p.parse_args()
    if a.action=='site':
        if not a.cache:p.error('--cache is required for the website crawl')
        crawl(a.output,a.cache)
    else:inventory(a.output)

if __name__=='__main__':main()
