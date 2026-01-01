#!/usr/bin/env python3
"""
compute_aun.py – Calculate the auN statistic from a GFA v1 file.

Definition
----------
auN = (Σ Li²) / (Σ Li), where Li are the lengths of all contigs.

Usage
-----
    python compute_aun.py path/to/assembly.gfa
"""

import argparse
import sys
from typing import Iterable, List, Tuple


def extract_segment_lengths(gfa_stream) -> Iterable[int]:
    """
    Yield the length of every segment in a GFA v1 stream.

    A segment (contig) line starts with 'S'.

    Field layout:
        1. 'S'
        2. Segment name
        3. Sequence (may be '*' if the sequence is external)
        4+. Optional tags, e.g. 'LN:i:<length>' which stores length.

    If the sequence is '*', the function searches for the LN:i tag.
    Lines without a determinable length are ignored with a warning.
    """
    for line in gfa_stream:
        if not line or line[0] != 'S':
            continue  # Skip non-segment records

        parts = line.rstrip('\n').split('\t')
        if parts[0] != 'S':                       # Malformed line
            continue

        seq_field = parts[2]
        if seq_field != '*':
            yield len(seq_field)
            continue

        # Sequence is absent; fall back to LN:i tag
        length = None
        for tag in parts[3:]:
            if tag.startswith('LN:i:'):
                length = int(tag.split(':')[2])
                break

        if length is None:
            sys.stderr.write(
                f"Warning: segment without length skipped:\n{line}"
            )
        else:
            yield length


def compute_aun(lengths: List[int]) -> Tuple[float, int, int]:
    """
    Compute auN, total length, and segment count from a list of lengths.
    """
    if not lengths:
        raise ValueError("No segment lengths were found in the input.")

    total_length = sum(lengths)
    squared_sum = sum(l * l for l in lengths)
    aun_value = squared_sum / total_length
    return aun_value, total_length, len(lengths)


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(
        description="Compute the auN statistic from a GFA v1 assembly."
    )
    parser.add_argument(
        "gfa",
        help="Path to the input GFA file",
    )
    args = parser.parse_args(argv)

    with open(args.gfa, "r", encoding="utf-8") as fh:
        lengths = list(extract_segment_lengths(fh))
        aun, total_bp, segment_count = compute_aun(lengths)

    print(f"File: {args.gfa}")
    print(f"Segments: {segment_count}")
    print(f"Total length: {total_bp:,} bp")
    print(f"auN: {aun:,.2f} bp")


if __name__ == "__main__":
    main()

