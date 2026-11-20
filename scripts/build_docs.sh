#!/usr/bin/env bash
# Build and check the current manual without modifying its Markdown sources.
# Requires Python 3, FORD (tested with 7.0.9), and Graphviz; installs nothing.
set -euo pipefail
for dependency in python3 ford dot; do
    if ! command -v "$dependency" >/dev/null 2>&1; then
        echo "Required documentation dependency not found: $dependency" >&2
        exit 1
    fi
done
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
python3 - "$script_dir/.." <<'PY'
from pathlib import Path
from urllib.parse import urlsplit, urlunsplit
import json
import re
import shutil
import subprocess
import sys
import tempfile

root = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(root / 'tests'))
import check_user_manual as checker


def local_target(target):
    """Change only local Markdown destinations; keep query strings/fragments."""
    url = urlsplit(target)
    if url.scheme or url.netloc or not url.path.endswith('.md'):
        return target
    path = 'index.html' if url.path in ('README.md', './README.md') else url.path[:-3] + '.html'
    return urlunsplit(('', '', path, url.query, url.fragment))


def convert_links(text):
    # Skip fenced and inline code. Support inline, reference and raw HTML links.
    fence = None
    lines = []
    for line in text.splitlines(keepends=True):
        marker = re.match(r'^\s{0,3}(`{3,}|~{3,})', line)
        if marker:
            token = marker.group(1)
            if fence is None:
                fence = token
            elif token[0] == fence[0] and len(token) >= len(fence):
                fence = None
            lines.append(line)
            continue
        if fence:
            lines.append(line)
            continue
        chunks = re.split(r'(`+[^`]*`+)', line)
        for i in range(0, len(chunks), 2):
            chunk = chunks[i]
            chunk = re.sub(r'(\]\(<?)([^\s)>]+)',
                           lambda m: m[1] + local_target(m[2]), chunk)
            chunk = re.sub(r'(^\s{0,3}\[[^\]]+\]:\s*<?)([^\s>]+)',
                           lambda m: m[1] + local_target(m[2]), chunk)
            chunk = re.sub(r'((?:href|src)=["\'])([^"\']+)',
                           lambda m: m[1] + local_target(m[2]), chunk)
            chunks[i] = chunk
        lines.append(''.join(chunks))
    return ''.join(lines)


def project_config(staging, output):
    _, metadata, body = (root / 'ford_project.md').read_text().split('---', 2)
    # These historical settings are unsupported or overridden for local preview.
    replaced = {'srd_dir', 'base_url', 'ordered_subpage', 'page_dir',
                'output_dir', 'src_dir', 'favicon', 'project_url'}
    kept = []
    skip = False
    for line in metadata.splitlines():
        match = re.match(r'^(\w+):', line)
        if match:
            skip = match[1] in replaced
        if not skip and line.strip():
            kept.append(line)
    kept.extend([f'src_dir: {root / "src"}', f'page_dir: {staging}',
                 f'output_dir: {output}', f'favicon: {root / "favicon.ico"}',
                 'project_url:'])
    return '---\n' + '\n'.join(kept) + '\n---\n' + body


def api_issues(site):
    """Report generated API defects separately from the required manual check."""
    cache = {}
    issues = []
    for directory in ('proc', 'module', 'interface', 'type', 'program', 'sourcefile', 'lists'):
        for page in (site / directory).glob('*.html'):
            parsed = checker.Document(page.read_text())
            for target in parsed.links:
                url = urlsplit(target)
                if url.scheme or url.netloc:
                    continue
                from urllib.parse import unquote
                destination = (page.parent / unquote(url.path)).resolve() if url.path else page
                if destination.is_dir():
                    destination /= 'index.html'
                error = None
                if not destination.is_file():
                    error = 'missing file'
                elif destination.suffix == '.html' and url.fragment:
                    if destination not in cache:
                        cache[destination] = checker.Document(destination.read_text()).ids
                    if unquote(url.fragment) not in cache[destination]:
                        error = 'missing anchor'
                if error:
                    issues.append(dict(page=str(page.relative_to(site)), target=target, error=error))
    return issues


try:
    with tempfile.TemporaryDirectory(prefix='.docs-build-', dir=root) as temporary:
        work = Path(temporary)
        pages = work / 'pages'
        pages.mkdir()
        site = work / 'site'
        for name in checker.MANUAL_PAGES:
            source = root / name
            text = source.read_text()
            if name != 'docs/index.md':
                title = re.search(r'^# (.+)$', text, re.M)
                if not title:
                    raise RuntimeError(f'Missing page title: {name}')
                label = 'Plasma Models' if source.name == 'models.md' else title[1]
                text = f'title: {label}\n\n' + text
            (pages / source.name).write_text(convert_links(text))
        config = work / 'ford_project.md'
        config.write_text(project_config(pages, site))
        log = work / 'build.log'
        print('Building the complete PLUME manual and Fortran API with FORD...', flush=True)
        with log.open('w') as stream:
            build = subprocess.run(['ford', str(config)], cwd=root, stdout=stream,
                                   stderr=subprocess.STDOUT)
        if build.returncode:
            print(log.read_text(), file=sys.stderr)
            raise RuntimeError(f'FORD failed with exit status {build.returncode}')
        shutil.copytree(root / 'docs/assets', site / 'page/assets')
        result = {'source': checker.source_checks(), 'site': checker.site_checks(site)}
        if result['source']['errors'] or result['site']['missing_pages'] or result['site']['broken_links']:
            print(json.dumps(result, indent=2), file=sys.stderr)
            raise RuntimeError('Manual validation failed; previous docs_out was preserved')
        issues = api_issues(site)
        # Store stable final paths in the report rather than the temporary build path.
        result['site']['site'] = str(root / 'docs_out')
        result['api_issues'] = issues
        (site / 'validation.json').write_text(json.dumps(result, indent=2) + '\n')
        shutil.copy2(log, site / 'build.log')
        destination = root / 'docs_out'
        backup = None
        if destination.exists() or destination.is_symlink():
            backup = Path(tempfile.mkdtemp(prefix='docs_out.backup-', dir=root))
            backup.rmdir()
            destination.rename(backup)
        try:
            site.rename(destination)
        except BaseException:
            if backup is not None:
                backup.rename(destination)
            raise
        if backup is not None:
            print(f'Previous output preserved: {backup}')
        print(f'Validated {len(checker.MANUAL_PAGES)} manual pages: no broken manual links or missing assets.')
        print(f'Generated API link issues: {len(issues)} (see docs_out/validation.json).')
        print('FORD warnings and build details: docs_out/build.log')
        print('MathJax uses an external CDN; external website availability is not checked.')
        print(f'Homepage: {destination / "page/index.html"}')
        print('Open from the repository root: xdg-open docs_out/page/index.html')
except (OSError, RuntimeError, ValueError) as error:
    sys.exit(f'Documentation build failed: {error}')
PY
