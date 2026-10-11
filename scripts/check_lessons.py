#!/usr/bin/env python3
"""Check every lesson under curriculum/ against the structure LESSON_TEMPLATE.md requires.

For each curriculum/phase-XX/lesson-YY-<slug>/ (templates excluded) this verifies:
  * the eight files exist (lesson.md, problems.md, starter/main.py, tests/test_main.py and, in
    solutions/phase-XX/lesson-YY-<slug>/, hints.md, method.md, worked-solutions.md, solution.py);
  * lesson.md has every metadata row, the six part headings in order, an adversarial self-checks
    section and a "Check yourself" section, and says "Time | unmeasured";
  * every `@id` cited in lesson.md exists in references/sources.yaml and is not optional-paid in a
    required position (a line containing "required" and an optional-paid id fails);
  * lesson.md never claims a credential ("degree-equivalent", "certified", "accredited" outside a
    negation are flagged for a human look);
  * the practical-activity row names one of the four modes.
Exit status 1 on any error. Warnings do not fail the check.
"""
import re
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
META_ROWS = ['Learning outcomes', 'Prerequisites', 'Engineering question', 'Assumptions and limits', 'Analytical work',
             'Computational deliverable', 'Practical activity', 'Equipment', 'Safety', 'Evidence', 'Assessment',
             'References', 'Extensions', 'Time']
PARTS = ['## 1. The engineering problem', '## 2. Understand and derive', '## 3. Build the computational model',
         '## 4. Experiment at home, or investigate real data', '## 5. Validate and improve', '## 6. Defend your work']
MODES = ['physical experiment', 'computational experiment', 'data-based investigation', 'design-build']
CLAIM_WORDS = re.compile(r'degree-equivalent|accredited|certif(?:ied|icate)', re.I)


def main():
    sources = yaml.safe_load((ROOT / 'references' / 'sources.yaml').read_text())['sources']
    tiers = {s['id']: s['access_tier'] for s in sources}
    errors, warnings = [], []
    lesson_dirs = sorted(p for p in (ROOT / 'curriculum').glob('phase-*/lesson-*') if p.is_dir())
    for d in lesson_dirs:
        rel = d.relative_to(ROOT / 'curriculum')
        sol = ROOT / 'solutions' / rel
        for f in ['lesson.md', 'problems.md', 'starter/main.py', 'tests/test_main.py']:
            if not (d / f).exists():
                errors.append(f'{rel}: missing {f}')
        for f in ['hints.md', 'method.md', 'worked-solutions.md', 'solution.py']:
            if not (sol / f).exists():
                errors.append(f'{rel}: missing solutions/{rel}/{f}')
        if not (d / 'lesson.md').exists():
            continue
        text = (d / 'lesson.md').read_text()
        for row in META_ROWS:
            if f'**{row}**' not in text:
                errors.append(f'{rel}: metadata row "{row}" missing')
        pos = -1
        for part in PARTS:
            i = text.find(part)
            if i < 0:
                errors.append(f'{rel}: heading "{part}" missing')
            elif i < pos:
                errors.append(f'{rel}: heading "{part}" out of order')
            pos = max(pos, i)
        if not re.search(r'#+\s*Adversarial self-checks', text):
            errors.append(f'{rel}: "Adversarial self-checks" section missing')
        if not re.search(r'## Check yourself', text):
            errors.append(f'{rel}: "## Check yourself" section missing')
        if not re.search(r'\*\*Time\*\*\s*\|\s*unmeasured', text):
            errors.append(f'{rel}: Time row must read "unmeasured" until pilots measure it')
        m = re.search(r'\*\*Practical activity\*\*\s*\|([^\n]*)', text)
        if m and not any(mode in m.group(1).lower() for mode in MODES):
            errors.append(f'{rel}: practical-activity row names none of {MODES}')
        for sid in sorted(set(re.findall(r'(?<![\w.])@([a-z0-9][a-z0-9-]+)', text))):
            if sid not in tiers:
                errors.append(f'{rel}: unknown source @{sid}')
        for line in text.splitlines():
            for sid in re.findall(r'(?<![\w.])@([a-z0-9][a-z0-9-]+)', line):
                if tiers.get(sid) == 'optional-paid' and 'required' in line.lower():
                    errors.append(f'{rel}: optional-paid source @{sid} appears on a line marked required')
        for mm in CLAIM_WORDS.finditer(text):
            ctx = text[max(0, mm.start() - 60):mm.end() + 60].replace('\n', ' ')
            if not re.search(r'\b(not|no|never|does not|is not|without)\b', ctx, re.I):
                warnings.append(f'{rel}: credential-like word "{mm.group(0)}" near: ...{ctx}...')
    print(f'{len(lesson_dirs)} lessons checked')
    for w in warnings:
        print('WARNING', w)
    for e in errors:
        print('ERROR', e)
    if errors:
        sys.exit(1)
    print('lesson structure OK')


if __name__ == '__main__':
    main()
