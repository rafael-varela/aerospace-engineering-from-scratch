#!/usr/bin/env python3
"""Validate references/sources.yaml and the per-phase reference maps.

Checks (always):
  * every source has the required fields and a valid access_tier;
  * every `access_tier: required-free` or `supplementary-free` source has a license recorded;
  * every source id referenced in references/phase-*.md exists in sources.yaml;
  * every source has a verification date no older than 18 months.

With --links, additionally sends a HEAD/GET request to every url and reports non-2xx/3xx
responses. A working link is reachability only, not content verification.
"""
import argparse
import datetime as dt
import re
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
SOURCES = ROOT / 'references' / 'sources.yaml'
REQUIRED = {'id', 'title', 'creator', 'url', 'license', 'access_tier', 'verified'}
TIERS = {'required-free', 'supplementary-free', 'optional-paid'}
MAX_AGE_DAYS = 548


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--links', action='store_true', help='also check that every url responds')
    args = ap.parse_args()

    data = yaml.safe_load(SOURCES.read_text()) or {}
    sources = data.get('sources', [])
    errors = []
    ids = set()
    today = dt.date.today()

    for src in sources:
        sid = src.get('id', '<missing id>')
        missing = REQUIRED - set(src)
        if missing:
            errors.append(f'{sid}: missing fields {sorted(missing)}')
        if sid in ids:
            errors.append(f'{sid}: duplicate id')
        ids.add(sid)
        if src.get('access_tier') not in TIERS:
            errors.append(f'{sid}: access_tier must be one of {sorted(TIERS)}')
        if src.get('access_tier') == 'optional-paid' and not src.get('value_added'):
            errors.append(f'{sid}: optional-paid sources must say what value they add (value_added)')
        v = src.get('verified')
        if isinstance(v, dt.date) and (today - v).days > MAX_AGE_DAYS:
            errors.append(f'{sid}: verification date {v} is older than {MAX_AGE_DAYS} days; re-verify')

    for md in sorted((ROOT / 'references').glob('phase-*.md')) + sorted((ROOT / 'references').glob('bridge-*.md')):
        for ref in set(re.findall(r'`([a-z0-9][a-z0-9-]+)`', md.read_text())):
            if ref.startswith('p') and re.fullmatch(r'p\d\d-\d\d', ref):
                continue  # lesson ids
            if ref not in ids and re.fullmatch(r'[a-z0-9-]+-(?:v\d|\d{4}|2e|3e|ocw|nasa|nist|bipm|docs)[a-z0-9-]*|[a-z]+-[a-z0-9-]+', ref):
                # only flag things that look like source ids (contain a hyphen) and are not known
                if ref in KNOWN_NON_SOURCES:
                    continue
                errors.append(f'{md.name}: references unknown source id `{ref}`')

    if args.links:
        import requests  # optional dependency
        for src in sources:
            url = src.get('url')
            try:
                r = requests.head(url, allow_redirects=True, timeout=20, headers={'User-Agent': 'aefs-linkcheck'})
                if r.status_code >= 400:
                    r = requests.get(url, allow_redirects=True, timeout=30, headers={'User-Agent': 'aefs-linkcheck'})
                if r.status_code >= 400:
                    errors.append(f'{src["id"]}: {url} returned {r.status_code}')
            except Exception as e:  # noqa: BLE001
                errors.append(f'{src["id"]}: {url} failed: {e}')

    print(f'{len(sources)} sources checked')
    for e in errors:
        print('ERROR', e)
    if errors:
        sys.exit(1)
    print('reference manifest OK')


# identifiers that appear in backticks in the maps but are not source ids
KNOWN_NON_SOURCES = {'required-free', 'supplementary-free', 'optional-paid', 'sources-yaml', 'check-references', 'lesson-md'}

if __name__ == '__main__':
    main()
