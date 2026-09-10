#!/usr/bin/env python3
"""Summarise run_all.sh output into a markdown report.

Usage: ./sum_all.py [out dir] [report.md]      defaults: out-all, <out dir>/summary.md
Reads <graph>_<scheme>_<metric>_<load>.csv, e.g. USA_FHL-inregion_time_rank.csv.
"""
import datetime, pathlib, re, statistics as stat, sys
OUT = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else 'out-all')
REPORT = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 else OUT / 'summary.md'
GRAPHS = ['NY', 'FLA', 'CAL', 'E', 'W', 'USA']
SCHEMES = [('CCH', 'CCH (baseline)'), ('CTL', 'CTL'), ('CTNR', 'CTNR (full table)'),
           ('FHL', '**FHL**'), ('FHL-inregion', 'FHL + in-region hubs'),
           ('FHL-toplabels', 'FHL + top labels (K=3)'),
           ('FHL-inregion-toplabels', 'FHL + in-region hubs + top labels'),
           ('FHL-overlay', 'FHL overlay')]
MEM, CUST = r'total: (\d+) MB', r'customization in (\d+) micro'
md = []

def rows(run):  # one dict per query row. Columns are looked up BY NAME because the layout
    # differs per mode: CCH writes distance,query_time; CTNR adds a mode column; the rank
    # load prepends dijkstra_rank. Non-numeric values (mode) are dropped here.
    f = OUT / (run + '.csv')
    if not f.exists() or f.stat().st_size == 0:
        return None
    lines = [l for l in f.read_text().splitlines() if l[:1] != '#']
    head = lines[0].split(',')
    return [{k: int(v) for k, v in zip(head, l.split(',')) if v.isdigit()} for l in lines[1:]]

def times(run, rank=None):  # mean query time, microseconds. A rank keeps only that bucket,
    # which also drops the rank-0 warm-up prefix.
    rs = rows(run)
    if not rs:
        return '-'
    v = [r['query_time'] for r in rs if rank in (None, r.get('dijkstra_rank'))]
    return '%.2f' % (sum(v) / len(v) / 1e3)

def grep(name, pat):  # pull one field out of a csv comment header or a log; '?' if missing
    f = OUT / name
    m = re.search(pat, f.read_text()) if f.exists() else None
    return m.group(1) if m else '?'

def memory(g, s):  # total memory of the query-time structures, MB
    v = grep('%s_%s_time_long.csv' % (g, s), MEM)
    return v if v.isdigit() else '-'

def customization(g, s):  # full customization, seconds: median of the -n runs in the query run
    f = OUT / ('%s_%s_time_long.csv' % (g, s))
    m = re.search(r'# Customization: ([\d.]+) ms \(median of (\d+)\)', f.read_text()) if f.exists() else None
    if m and int(m.group(2)) > 1:
        return '%.2f' % (float(m.group(1)) / 1e3)
    # Older batches (no repeated runs recorded): FHL/CTNR from the one logged customization,
    # CCH/CTL from the separate -a CTL-custom run (total_time, which includes the CCH part).
    # Marked with * because that is not the same measurement.
    if s in ('CCH', 'CTL'):
        rs = rows('%s_customization_time' % g)
        us = stat.median(r['cch_customization' if s == 'CCH' else 'total_time']
                         for r in rs) if rs else None
    else:
        v = grep('%s_%s_time_long.log' % (g, s), CUST)
        us = int(v) if v.isdigit() else None
    return '%.2f*' % (us / 1e6) if us else '-'

def table(title, cell, note=None):
    md.extend(['', '### ' + title, ''] + ([note, ''] if note else []) +
              ['| Scheme | ' + ' | '.join(GRAPHS) + ' |', '|---|' + '---|' * len(GRAPHS)])
    md.extend('| %s | ' % label + ' | '.join(cell(g, s) for g in GRAPHS) + ' |'
              for s, label in SCHEMES)

# ---- header: what produced these numbers ---------------------------------------------------
md += ['# Benchmark results', '',
       'Generated %s from `%s`.' % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M'), OUT), '']
cfg = OUT / 'run-config.txt'
if cfg.exists():
    md += ['Run configuration:', '', '```', cfg.read_text().rstrip(), '```']
else:  # batches from before run_all.sh recorded its settings: recover what the logs show
    th = sorted({m for f in OUT.glob('*.log')
                 for m in re.findall(r'Using (\d+) threads', f.read_text())})
    md += ['Run configuration: `run-config.txt` not found (batch predates it). Threads seen in '
           'the logs: %s. NUMA policy: not recorded.' % (', '.join(th) or 'unknown')]
md += ['', 'Times are in microseconds per query. Every scheme is checked query by query '
       'against CCH-tree; see the correctness section at the end.']

table('Long-distance queries (10,000 uniform random pairs; mean)',
      lambda g, s: times('%s_%s_time_long' % (g, s)))
table('Distance metric, long-distance queries (mean)',
      lambda g, s: times('%s_%s_dist_long' % (g, s)))
table('Memory (MB)', memory)
table('Full customization (s; median of 3 customizations in each scheme\'s own run)', customization,
      note='Unlike preprocessing, nothing is added here: every scheme\'s customization starts by '
           're-customizing the CCH, so its number already includes the CCH row. (Cells marked * '
           'come from an older batch that recorded a single run, or for CCH/CTL a separate '
           '-a CTL-custom run.)')

# ---- preprocessing -----------------------------------------------------------------------------
# Read from the csv header the launcher writes (# Preprocessing CCH / scheme). Customization has
# its own table above. Batches from before these headers existed show '-'.
PHASE = r'# Preprocessing CCH: ([\d.]+) ms \((built|loaded from cache)\)', r'# Preprocessing scheme: ([\d.]+) ms'

def phases(g, s):  # (CCH ms, 'built' or 'loaded from cache', scheme ms), or None
    f = OUT / ('%s_%s_time_long.csv' % (g, s))
    t = f.read_text() if f.exists() else ''
    m = [re.search(p, t) for p in PHASE]
    return None if not all(m) else (float(m[0].group(1)), m[0].group(2), float(m[1].group(1)))

md.extend(['', '### Preprocessing (s; metric-independent, one run each)', '',
           'All schemes are built on the CCH: a scheme\'s total preprocessing time is its row '
           'plus the CCH row.', '',
           '| Scheme | ' + ' | '.join(GRAPHS) + ' |', '|---|' + '---|' * len(GRAPHS)])
for s, label in SCHEMES:
    cells = []
    for g in GRAPHS:
        p = phases(g, s)
        if s == 'CCH':  # CCH's preprocessing IS the CCH construction; nothing sits on top of it
            cells.append('%.2f' % (p[0] / 1e3) if p and p[1] == 'built' else ('(cached)' if p else '-'))
        else:
            cells.append('%.2f' % (p[2] / 1e3) if p else '-')
    md.append('| %s | ' % label + ' | '.join(cells) + ' |')
md.append('')

md += ['', '## Query time by Dijkstra rank (mean)']
for g in GRAPHS:
    ref = rows('%s_CCH_time_rank' % g)
    if not ref:
        continue
    md += ['', '### ' + g, '', '| Rank | ' + ' | '.join(l.replace('**', '') for _, l in SCHEMES) + ' |',
           '|---|' + '---|' * len(SCHEMES)]
    for r in sorted({x['dijkstra_rank'] for x in ref} - {0}):
        md.append('| 2^%d | ' % (r.bit_length() - 1) + ' | '.join(
            times('%s_%s_time_rank' % (g, s), r) for s, _ in SCHEMES) + ' |')

md += ['', '## Correctness', '', 'Mismatching distances / queries compared, per scheme, against '
       'CCH-tree over all three loads.', '', '| Graph | ' + ' | '.join(s for s, _ in SCHEMES[1:]) + ' |',
       '|---|' + '---|' * (len(SCHEMES) - 1)]
total = 0
for g in GRAPHS:
    cells = []
    for s, _ in SCHEMES[1:]:                             # everything except CCH itself
        n = k = 0
        for tag in ('time_long', 'time_rank', 'dist_long'):
            mine, ref = rows('%s_%s_%s' % (g, s, tag)), rows('%s_CCH_%s' % (g, tag))
            if mine and ref:
                k += min(len(mine), len(ref))
                n += sum(a['distance'] != b['distance'] for a, b in zip(mine, ref))
        total += n
        cells.append('%d / %d' % (n, k) if k else '-')
    md.append('| %s | ' % g + ' | '.join(cells) + ' |')
md += ['', '**Total mismatches: %d**' % total, '']

REPORT.write_text('\n'.join(md))
print('wrote %s  (total mismatches: %d)' % (REPORT, total))
