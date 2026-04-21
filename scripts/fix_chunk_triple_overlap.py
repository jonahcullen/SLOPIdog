#!/usr/bin/env python3
"""
Fix triple input-region overlaps in GLIMPSE2 chunk TSV files.

GLIMPSE2_ligate reads the input_region stored in each imputed BCF's header
and requires that no genomic position is covered by more than 2 chunks'
input regions simultaneously. When three consecutive chunks all overlap in
their input regions (which can happen when a middle chunk has an unusually
large output window, causing the flanking chunks' buffers to reach each
other), GLIMPSE2_ligate fails.

Fix strategy: for each consecutive triple (i, i+1, i+2) where chunk i's
input_end >= chunk i+2's input_start, trim chunk i's input_end to
chunk i+2's input_start - 1. If that would cut into chunk i's output_end,
trim chunk i+2's input_start instead. Repeat until no triple overlaps remain.
"""

import sys
import argparse


def parse_region(region_str):
    """Parse 'chr:start-end' -> (chrom, start, end) with integer coords."""
    chrom, coords = region_str.rsplit(":", 1)
    start, end = coords.split("-")
    return chrom, int(start), int(end)


def format_region(chrom, start, end):
    return f"{chrom}:{start}-{end}"


def read_chunks(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            rows.append(line.split("\t"))
    return rows


def write_chunks(rows, path):
    with open(path, "w") as fh:
        for row in rows:
            fh.write("\t".join(row) + "\n")


def fix_triple_overlaps(rows):
    """
    Iteratively trim input regions until no three consecutive chunks share
    any genomic position in their input regions.

    Returns (fixed_rows, n_fixes).
    """
    rows = [list(r) for r in rows]  # deep copy
    n_fixes = 0
    changed = True

    while changed:
        changed = False
        for i in range(len(rows) - 2):
            c0_in_chrom, c0_in_s, c0_in_e = parse_region(rows[i][2])
            _,           c2_in_s, _       = parse_region(rows[i + 2][2])

            # Triple overlap: chunk i's input end reaches into chunk i+2's input
            if c0_in_e >= c2_in_s:
                _, c0_out_s, c0_out_e = parse_region(rows[i][3])

                new_end = c2_in_s - 1

                if new_end >= c0_out_e:
                    # Safe to trim chunk i's input end
                    rows[i][2] = format_region(c0_in_chrom, c0_in_s, new_end)
                    n_fixes += 1
                    changed = True
                else:
                    # Trimming chunk i would cut into its output region.
                    # Instead, trim chunk i+2's input start.
                    c2_in_chrom, c2_in_s2, c2_in_e = parse_region(rows[i + 2][2])
                    _, c2_out_s, c2_out_e = parse_region(rows[i + 2][3])

                    new_start = c0_in_e + 1

                    if new_start <= c2_out_s:
                        # Safe to trim chunk i+2's input start
                        rows[i + 2][2] = format_region(c2_in_chrom, new_start, c2_in_e)
                        n_fixes += 1
                        changed = True
                    else:
                        # Neither trim is safe — degenerate case, warn and skip
                        print(
                            f"WARNING: cannot safely fix triple input overlap for chunks "
                            f"{rows[i][0]}, {rows[i+1][0]}, {rows[i+2][0]} "
                            f"(trimming would cut into an output region). "
                            f"Consider adjusting --window-mb / --buffer-mb.",
                            file=sys.stderr,
                        )

    return rows, n_fixes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input",  required=True, help="Input chunk TSV (from GLIMPSE2_chunk)")
    parser.add_argument("-o", "--output", required=True, help="Output chunk TSV with overlaps fixed")
    args = parser.parse_args()

    rows = read_chunks(args.input)
    fixed, n_fixes = fix_triple_overlaps(rows)

    if n_fixes:
        print(
            f"Fixed {n_fixes} triple input-region overlap(s) in {args.input}",
            file=sys.stderr,
        )
    else:
        print(f"No triple input-region overlaps found in {args.input}", file=sys.stderr)

    write_chunks(fixed, args.output)


if __name__ == "__main__":
    main()
