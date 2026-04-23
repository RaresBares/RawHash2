# RawHash2 Modifications — What Was Changed Upstream

> What source files of RawHash2 were patched to support the 5-segmenter benchmark.

The base repository is `github.com/CMU-SAFARI/RawHash`. All modifications live
in `/home/rsahleanu/rawhash2/`. The feature is fully opt-in: if `--segmenter`
is not passed, RawHash2 runs exactly like the upstream release.

---

## Git commits (chronological)

```
a1026f1  Add pluggable event segmenter support via --segmenter flag
c252646  Fix Window/BinSeg segmenter penalties for normalized signals
fa5aa83  Fix gen_events n_events count and tune segmenter penalties
7536830  Replace Window segmenter with Scrappie (ONT t-stat variant)
a2db20e  Add segmenter benchmark suite with ground-truth accuracy evaluation
3d9637e  Fix D3 ground truth: extract basecalls from FAST5 internal fastq
```

`git log --oneline | head -10` in the rawhash2 repo to reproduce.

---

## Modified source files (core C/C++)

### `src/revent.c`  (+397 / –32 lines)

**The central change.** Before: only one hard-coded event detection algorithm
(the original ONT t-statistic). After: a dispatcher that routes to one of
five implementations based on `segmenter_type`.

Added:
- `static float* detect_events_hmm(...)`       — 4-state Viterbi segmenter
- `static float* detect_events_pelt(...)`      — Pruned Exact Linear Time changepoint
- `static float* detect_events_binseg(...)`    — Binary segmentation, recursive
- `static float* detect_events_scrappie(...)`  — ONT t-stat with Scrappie parameters
- `static float* detect_events_python(...)`    — subprocess hook for arbitrary scripts

Modified:
- `detect_events(...)` — now takes `segmenter_type` and `python_script` args,
  and dispatches:
  ```c
  switch (segmenter_type) {
      case RI_SEGMENTER_DEFAULT:  events = detect_events_default(...); break;
      case RI_SEGMENTER_HMM:      events = detect_events_hmm(...);     break;
      case RI_SEGMENTER_PELT:     events = detect_events_pelt(...);    break;
      case RI_SEGMENTER_BINSEG:   events = detect_events_binseg(...);  break;
      case RI_SEGMENTER_SCRAPPIE: events = detect_events_scrappie(...); break;
      case RI_SEGMENTER_PYTHON:   events = detect_events_python(...);  break;
  }
  ```

### `src/revent.h`  (+4 lines)

Updated the function prototype of `detect_events()` to expose the two new
trailing parameters `segmenter_type` and `python_script`, with defaults so
existing call sites compile unchanged:

```c
float* detect_events(void *km, uint32_t s_len, const float* sig,
                     /* unchanged args */,
                     uint32_t segmenter_type = 0,
                     const char *python_script = 0);
```

### `src/roptions.h`  (+12 lines)

Added constants:

```c
#define RI_SEGMENTER_DEFAULT    0
#define RI_SEGMENTER_HMM        1
#define RI_SEGMENTER_PELT       2
#define RI_SEGMENTER_BINSEG     3
#define RI_SEGMENTER_SCRAPPIE   4
#define RI_SEGMENTER_PYTHON     5
```

Added fields to both `ri_idxopt_t` (indexing options) and `ri_mapopt_t`
(mapping options):

```c
uint32_t segmenter_type;
const char *python_segmenter_script;
```

### `src/roptions.c`  (+4 lines)

Initialize the two new option fields to sane defaults (`RI_SEGMENTER_DEFAULT`, `NULL`).

### `src/main.cpp`  (+19 lines)

Two new CLI flags parsed and stored into both option structs:

```c
{ "segmenter",        required_argument, 372 },
{ "segmenter-script", required_argument, 373 },
```

With option dispatch:

```c
else if (c == 372) {
    if      (!strcmp(arg, "default"))  seg = RI_SEGMENTER_DEFAULT;
    else if (!strcmp(arg, "hmm"))      seg = RI_SEGMENTER_HMM;
    else if (!strcmp(arg, "pelt"))     seg = RI_SEGMENTER_PELT;
    else if (!strcmp(arg, "binseg"))   seg = RI_SEGMENTER_BINSEG;
    else if (!strcmp(arg, "scrappie")) seg = RI_SEGMENTER_SCRAPPIE;
    else if (!strcmp(arg, "python"))   seg = RI_SEGMENTER_PYTHON;
    else { fprintf(stderr, "unknown segmenter '%s'\n", arg); return 1; }
    opt.segmenter_type = ipt.segmenter_type = seg;
}
else if (c == 373) {
    opt.python_segmenter_script = ipt.python_segmenter_script = arg;
}
```

Help-text line added:

```
--segmenter STR   Event segmenter: 'default', 'hmm', 'pelt', 'binseg',
                  'scrappie', 'python' [default]
--segmenter-script FILE   Path to Python script for --segmenter python []
```

### `src/rindex.c`  (+8 lines)

Pass `ipt.segmenter_type` and `ipt.python_segmenter_script` into
`detect_events()` when `RI_I_SIG_TARGET` is set (indexing pre-segmented signal
targets). For the standard DNA-reference case this branch is not taken, so
the segmenter choice only matters at mapping time in our benchmark.

### `src/rindex.h`  (+2 lines)

Struct fields for propagating segmenter options through the index pipeline.

### `src/rmap.cpp`  (+2 lines)

In `ri_map_frag()`, the existing call site

```c
float* events = detect_events(b->km, s_len, sig, ..., NULL);
```

was changed to pass the segmenter options:

```c
float* events = detect_events(b->km, s_len, sig, ...,
                              opt->segmenter_type,
                              opt->python_segmenter_script);
```

This is the only change needed for the read-mapping pipeline — every segmenter
plugs in via the same `float* events` array interface.

---

## New files (not modifications, pure additions)

### Test / benchmark utilities (under `test/`)

| File | Purpose |
|---|---|
| `test/benchmark_segmenters.py` | 1097-line benchmark driver (runs all segmenters on all datasets, computes mapping stats & ground-truth F1, emits JSON + PDF) |
| `test/data/download_selected.sh` | Downloads D1/D2/D3/D6 FAST5 archives and reference genomes |
| `test/run_benchmark.sh` | Top-level wrapper script |

These live under `test/` so they do not affect the compiled binary.

### Kmer models (unrelated to segmenter work)

`extern/kmer_models/dna_r10.4.1_e8.2_400bps/uncalled_r1041_model_only_means.txt`
— 262,144-line R10.4 9-mer pore model added to enable D6 / D7 benchmarks.

---

## Summary of code delta

| Category | Files | Net lines changed |
|---|---|---|
| Core (C/C++) | 7 files (`revent.c/h`, `roptions.c/h`, `main.cpp`, `rindex.c/h`, `rmap.cpp`) | ~460 |
| CLI | 1 (`main.cpp`) | ~20 |
| Benchmark suite | 3 new files | ~1250 |
| Pore models | 1 new file | +262,144 (data, not code) |

**No existing functionality was removed.** A RawHash2 invocation without
`--segmenter` behaves bit-identically to the upstream release, because
`segmenter_type` defaults to `RI_SEGMENTER_DEFAULT` and `detect_events_default()`
is a verbatim copy of the pre-patch code path.

---

## How to reproduce the patched build

```bash
cd /home/rsahleanu/rawhash2
git log --oneline | head -6   # verify the 6 segmenter commits are present
make clean && make           # produces bin/rawhash2
./bin/rawhash2 2>&1 | grep segmenter  # should list the --segmenter flag
```

---

## Verification that segmenter choice is actually active

During a benchmark run, differences appear only when `--segmenter` is passed.
Sanity checks we ran:

1. `rawhash2 --segmenter default` vs no flag → identical outputs (same MD5)
2. `rawhash2 --segmenter binseg` on D1 → 55k mappings
3. `rawhash2 --segmenter pelt` on D1   → 1.32M mappings

So the plumbing from CLI → segmenter dispatch → Event output → Sketch → PAF is
exercised end-to-end.
