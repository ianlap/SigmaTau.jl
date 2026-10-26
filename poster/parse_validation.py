"""Parse the per-tau tables in docs/src/validation.md into poster/data/validation.json."""
import re, json, os
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
md = open(os.path.join(ROOT, 'docs', 'src', 'validation.md')).read()
sections = re.split(r'\n### ', md)
out = {}
names = {'Overlapping Allan':'adev','Modified Allan':'mdev','Time':'tdev','Overlapping Hadamard':'hdev','Total':'totdev','Modified Total':'mtotdev','Hadamard Total':'htotdev'}
for sec in sections:
    title = sec.split('\n',1)[0].strip()
    if title not in names: continue
    rows = []
    for line in sec.split('\n'):
        if not line.startswith('|'): continue
        cells = [c.strip() for c in line.strip('|').split('|')]
        if not re.match(r'^\d+', cells[0]): continue
        tau = float(cells[0].rstrip('¹')); flip = cells[0].endswith('¹')
        s32 = float(cells[1]); at = float(cells[2])
        if len(cells) == 6:   # core: tau, s32, at, st, s32ci, stci
            st = float(cells[3]); s32ci = cells[4]; stci = cells[5]; raw = st
        else:                 # total: tau, s32, at(raw), st raw, st corrected, s32ci, stci
            raw = float(cells[3]); st = float(cells[4]); s32ci = cells[5]; stci = cells[6]
        lo, hi = [float(v) for v in re.findall(r'[-+0-9.e]+', stci)]
        rows.append(dict(tau=tau, stable32=s32, allantools=at, sigmatau=st, sigmatau_raw=raw, lo=lo, hi=hi, noise_flip=flip))
    out[names[title]] = rows
os.makedirs(os.path.join(HERE, 'data'), exist_ok=True)
json.dump(out, open(os.path.join(HERE, 'data', 'validation.json'), 'w'), indent=1)
for k,v in out.items():
    worst = max(abs(r['sigmatau_raw' if k in ('mtotdev',) else 'sigmatau']-r['stable32'])/r['stable32'] for r in v if not r['noise_flip'])
    print(k, len(v), 'rows; worst rel diff (non-flip rows) vs Stable32: %.2e' % worst)
