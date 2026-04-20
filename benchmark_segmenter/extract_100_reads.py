#!/usr/bin/env python3
"""Extract the first N reads from a dataset's FAST5 files into a new directory.

Works for both single-read and multi-read FAST5. For multi-read,
creates new multi-read FAST5 with just the first N reads (in 1 file).
For single-read, copies the first N files.

Usage: extract_100_reads.py <input_fast5_dir> <output_dir> [n_reads=100]
"""
import sys
import os
import shutil
import h5py

SRC_DIR = sys.argv[1]
OUT_DIR = sys.argv[2]
N_READS = int(sys.argv[3]) if len(sys.argv) > 3 else 100

os.makedirs(OUT_DIR, exist_ok=True)

# Find all fast5 files in source
all_f5 = sorted([os.path.join(SRC_DIR, f) for f in os.listdir(SRC_DIR) if f.endswith('.fast5')])
if not all_f5:
    print(f"No FAST5 in {SRC_DIR}")
    sys.exit(1)

# Detect: single-read (one read per file, no /read_xxx group) or multi-read
with h5py.File(all_f5[0], 'r') as h:
    keys = list(h.keys())
    is_multi = any(k.startswith('read_') for k in keys)

print(f"Source: {SRC_DIR}")
print(f"FAST5 files: {len(all_f5)}")
print(f"Format: {'multi-read' if is_multi else 'single-read'}")
print(f"Target reads: {N_READS}")

if not is_multi:
    # Single-read: just copy first N files
    for i, src in enumerate(all_f5[:N_READS]):
        dst = os.path.join(OUT_DIR, os.path.basename(src))
        if not os.path.exists(dst):
            shutil.copy2(src, dst)
    n_out = len([f for f in os.listdir(OUT_DIR) if f.endswith('.fast5')])
    print(f"Copied {n_out} single-read FAST5 files to {OUT_DIR}")

else:
    # Multi-read: iterate through files, collect N reads into a single output file
    out_path = os.path.join(OUT_DIR, 'subset_100reads.fast5')
    if os.path.exists(out_path):
        os.remove(out_path)

    n_written = 0
    with h5py.File(out_path, 'w') as out_h:
        for src in all_f5:
            if n_written >= N_READS:
                break
            with h5py.File(src, 'r') as in_h:
                # Copy top-level attrs if new file
                if n_written == 0:
                    for ak, av in in_h.attrs.items():
                        try:
                            out_h.attrs[ak] = av
                        except Exception:
                            pass
                for k in in_h.keys():
                    if k.startswith('read_'):
                        if n_written >= N_READS:
                            break
                        in_h.copy(k, out_h)
                        n_written += 1

    print(f"Wrote {n_written} reads to {out_path}")
