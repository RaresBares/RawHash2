#!/usr/bin/env python3
"""Extract all read IDs from a FAST5 directory (single or multi-read)."""
import sys, os
import h5py

SRC = sys.argv[1]
OUT = sys.argv[2]

def decode(v):
    return v.decode() if isinstance(v, bytes) else str(v)

ids = []
files = [f for f in sorted(os.listdir(SRC)) if f.endswith('.fast5')]
for f in files:
    path = os.path.join(SRC, f)
    try:
        with h5py.File(path, 'r') as h:
            keys = list(h.keys())
            # multi-read: top-level groups start with read_
            multi = any(k.startswith('read_') for k in keys)
            if multi:
                for k in keys:
                    if k.startswith('read_'):
                        ids.append(k[5:])  # strip "read_"
            else:
                # single-read: /Raw/Reads/Read_<N>/read_id attr
                if 'Raw' in h and 'Reads' in h['Raw']:
                    for rk in h['Raw/Reads']:
                        try:
                            rid = h[f'Raw/Reads/{rk}'].attrs.get('read_id')
                            if rid is not None:
                                ids.append(decode(rid))
                        except Exception:
                            pass
    except Exception as e:
        print(f"Warn: {f}: {e}", file=sys.stderr)

# Dedupe
seen = set(); unique = []
for i in ids:
    if i and i not in seen:
        seen.add(i); unique.append(i)

with open(OUT, 'w') as f:
    for i in unique:
        f.write(i + '\n')
print(f"Extracted {len(unique)} read IDs from {len(files)} FAST5 files to {OUT}")
