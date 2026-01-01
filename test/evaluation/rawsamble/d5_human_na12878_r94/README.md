# Rawsamble - d5_human_na12878_r94

We assume your current directory points to this directory ([`d5_human_na12878_r94`](./)) and your binary containing the Rawsamble functionality is in your `$PATH` and called `rawhash2`.

## Single-liner to Run Everything

The script below is a single liner to run: 1) Rawsamble for all-vs-all read overlapping (for raw nanopore signals), 2) minimap2 for all-vs-all read overlapping (for basecalled reads), 3) miniasm to generate assembly graphs from the output of (1) and (2), 4) evalation output used in our manuscript. 

The following command will use 64 threads. You can change the maximum threads to use by providing a different value than 64 below. You may also need to change some paths signal files (or to other files) if you do not have access them (e.g., if you cannot convert your FAST5 files to POD5 files).

```bash
bash run_rawhash2.sh 64
```
