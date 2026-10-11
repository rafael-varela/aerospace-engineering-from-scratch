#!/usr/bin/env python3
"""Verify that every learning path in ROADMAP.md contains the prerequisites of every lesson it includes.

Parses the lesson tables (the `Needs` column) and the learning-paths table directly from the
Markdown, so it needs no separate metadata file. Exit status 1 on any gap.

Rules it checks:
  * Every `Pxx.yy` prerequisite of every lesson in a path is in that path (Core 00-09 is in every path).
    A `Needs` cell of the form `A or B` is satisfied by either.
  * Every `Needs` reference points at a lesson that exists.
  * No non-★ lesson depends on a ★ lesson.
  * Each path's capstone rule names at least one design and one build-and-test capstone, and
    the path contains each capstone's prerequisites.
"""
import re
import sys
from pathlib import Path

ROADMAP = Path(__file__).resolve().parent.parent / 'ROADMAP.md'
CORE = {f'{i:02d}' for i in range(10)}
DASH = '[–-]'


def parse(md):
    lessons = {}          # 'P10.05' -> {'title', 'needs': [[alts], ...], 'star'}
    phase = None
    capstones = {}        # '05' -> {'type', 'needs'}
    for line in md.splitlines():
        m = re.match(r'## Phase (\d\d):', line)
        if m:
            phase = m.group(1)
            continue
        if re.match(r'## (Bridge|Learning|Home|Assessment|Outcome|Decisions|Benchmark|Phase overview)', line):
            phase = None
        m = re.match(r'\| (\d\d) \| ([^|]+)\| ([^|]*)\| ([^|]*)\|', line)
        if not (m and phase):
            continue
        num, title, build, needs = m.group(1), m.group(2).strip(), m.group(3), m.group(4).strip()
        if phase == '26':
            capstones[num] = {'type': build.strip(), 'needs': needs}
            continue
        groups = []
        needs_main = needs.split(';')[0]
        for part in re.split(r',\s*', needs_main):
            part = part.strip()
            if not part:
                continue
            alts = [a.strip() for a in re.split(r'\s+or\s+', part)]
            groups.append(alts)
        lessons[f'P{phase}.{num}'] = {'title': title, 'needs': groups, 'star': '★' in title}
    return lessons, capstones


def expand_phase_spec(spec, lessons):
    """'10 (01–05, 09–13), 17, 18–19, 20 (03)' -> set of lesson ids."""
    ids = set()
    spec = spec.strip()
    for m in re.finditer(r'(\d\d)' + DASH + r'(\d\d)(?!\s*\))|(\d\d)(?:\s*\(([^)]*)\))?', spec):
        if m.group(1):
            for ph in range(int(m.group(1)), int(m.group(2)) + 1):
                ids |= {k for k in lessons if k.startswith(f'P{ph:02d}.')}
        else:
            ph = m.group(3)
            if m.group(4):
                for rng in re.split(r',\s*', m.group(4)):
                    a, _, b = rng.partition('–') if '–' in rng else rng.partition('-')
                    a = a.strip(); b = (b or a).strip()
                    for n in range(int(a), int(b) + 1):
                        ids.add(f'P{ph}.{n:02d}')
            else:
                ids |= {k for k in lessons if k.startswith(f'P{ph}.')}
    return ids


def capstone_prereqs_ok(needs, path_ids, lessons):
    """Capstone Needs like 'P15, P22.10' or 'P20 (04–10), P21.17'. Text after ';' is optional."""
    needs = needs.split(';')[0]
    for part in re.split(r',\s*(?![^()]*\))', needs):  # commas outside parentheses only
        part = part.strip()
        if not part or part == 'Any':
            continue
        m = re.match(r'P(\d\d)(?:\s*\(([^)]*)\))?$', part)
        if m:
            req = expand_phase_spec(m.group(1) + (f' ({m.group(2)})' if m.group(2) else ''), lessons)
        else:
            req = {part}
        missing = sorted(req - path_ids)
        if missing:
            return missing
    return []


def main():
    md = ROADMAP.read_text()
    lessons, capstones = parse(md)
    errors = []

    for lid, info in lessons.items():
        for alts in info['needs']:
            for a in alts:
                if a not in lessons:
                    errors.append(f'{lid} needs {a}, which does not exist')
                elif lessons[a]['star'] and not info['star']:
                    errors.append(f'{lid} (mandatory) depends on ★ lesson {a}')

    core_ids = {k for k in lessons if k[1:3] in CORE}
    table = md[md.index('## Learning paths'):]
    table = table[:table.index('\n## ', 5)]  # only the learning-paths section
    rows = re.findall(r'^\| \*\*([^*]+)\*\* \| ([^|]+) \| ([^|]+) \| [^|]+\|$', table, flags=re.M)
    if not rows:
        errors.append('no learning-path rows found')
    for name, spec, caps in rows:
        ids = core_ids | expand_phase_spec(spec, lessons)
        for lid in sorted(ids):
            for alts in lessons[lid]['needs']:
                if not any(a in ids for a in alts):
                    errors.append(f'{name}: {lid} needs {" or ".join(alts)}, missing from path')
        cap_nums = re.findall(r'\b(\d\d)\b', caps)
        cap_ranges = re.findall(r'(\d\d)' + DASH + r'(\d\d)', caps)
        for a, b in cap_ranges:
            cap_nums += [f'{n:02d}' for n in range(int(a), int(b) + 1)]
        types = {capstones[c]['type'] for c in cap_nums if c in capstones}
        if not any(t.startswith('Design') for t in types) or not any(t.startswith('Build') for t in types):
            errors.append(f'{name}: capstone rule "{caps}" lacks a design or a build-and-test option')
        for c in sorted(set(cap_nums)):
            if c not in capstones:
                errors.append(f'{name}: capstone {c} does not exist')
                continue
            missing = capstone_prereqs_ok(capstones[c]['needs'], ids, lessons)
            if missing:
                errors.append(f'{name}: capstone {c} needs {missing[:4]}{"..." if len(missing) > 4 else ""}, missing from path')

    print(f'{len(lessons)} lessons, {len(capstones)} capstones, {len(rows)} paths')
    for e in errors:
        print('ERROR', e)
    if errors:
        sys.exit(1)
    print('all paths are prerequisite-complete')


if __name__ == '__main__':
    main()
