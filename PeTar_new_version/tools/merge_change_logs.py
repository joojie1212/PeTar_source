#!/usr/bin/env python3
"""Merge a stopped/quiescent run's data.change.<rank>.<thread> shards.
Output is sorted by time_new, time_old, rank, thread, row; source files are retained.
Use --include-legacy when continuing a run that previously wrote data.change.
"""
import argparse
from pathlib import Path
import re
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('prefix',type=Path,help='e.g. /run/data.change')
p.add_argument('--output',type=Path,required=True)
p.add_argument('--include-legacy',action='store_true')
a=p.parse_args()
sources=sorted(f for f in a.prefix.parent.glob(a.prefix.name+'.*')
               if re.fullmatch(re.escape(a.prefix.name)+r'\.\d+\.\d+',f.name))
if a.include_legacy and a.prefix.exists(): sources.insert(0,a.prefix)
if not sources:p.error('no change log sources found')
if any(f.resolve()==a.output.resolve() for f in sources):p.error('output must not overwrite an input')
rows=[];header=None
for f in sources:
    with f.open() as stream:
        h=next(stream).strip()
        if header is None:header=h
        if h!=header:raise ValueError(f'inconsistent header: {f}')
        for i,line in enumerate(stream):
            fields=line.split()
            if not line.endswith('\n') or len(fields)!=19:raise ValueError(f'incomplete row: {f}:{i+2}')
            rows.append((float(fields[1]),float(fields[0]),f.name,i,line))
with a.output.open('x') as out:
    out.write(header+'\n')
    for *_,line in sorted(rows):out.write(line)
print(f'{len(rows)} rows from {len(sources)} files -> {a.output}')
