#!/usr/bin/env python3
# naive_2mm.py
# Naïve approximate matcher allowing up to 2 mismatches.
# Requirements satisfied:
# - Manual FASTA/FASTQ parsing (no Biopython)
# - Phred+33 quality handling; keep reads with all Q > 5
# - Report read id and all start positions in reference (<=2 mismatches)
# - Note how many reads do NOT make the cut (quality filter)

import sys
import argparse

def phred33_to_q(ch):
    """Convert a single FASTQ quality char to integer Phred score (Sanger / +33)."""
    return ord(ch) - 33

def parse_fasta(path):
    """Read a (possibly multi-record) FASTA and return the concatenated uppercase sequence.
       Also returns the first header for reference."""
    header = None
    seq_chunks = []
    with open(path, 'r', encoding='utf-8') as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith('>'):
                if header is None:
                    header = line[1:].strip()
                # allow multiple records; we just concatenate
                continue
            # accumulate sequence lines, stripping whitespace and uppercasing
            line = line.strip()
            if line:
                # manual uppercase without locale issues
                up = []
                for c in line:
                    # uppercase via ASCII logic (avoid .upper())
                    if 'a' <= c <= 'z':
                        up.append(chr(ord(c) - 32))
                    else:
                        up.append(c)
                seq_chunks.append(''.join(up))
    return (header if header is not None else "reference"), ''.join(seq_chunks)

def parse_fastq(path):
    """Generator over FASTQ records: yields (read_id, seq, qual_scores_list)."""
    with open(path, 'r', encoding='utf-8') as fh:
        while True:
            header = fh.readline()
            if not header:
                break  # EOF
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()

            if not qual:
                # malformed/truncated record
                break

            # parse read id (up to first whitespace)
            read_id = header.strip()
            if read_id.startswith('@'):
                read_id = read_id[1:]
            if read_id == '':
                read_id = 'unnamed_read'

            # clean newline
            seq = seq.strip()
            qual = qual.strip()

            # uppercase sequence manually
            up_seq_chars = []
            for c in seq:
                if 'a' <= c <= 'z':
                    up_seq_chars.append(chr(ord(c) - 32))
                else:
                    up_seq_chars.append(c)
            up_seq = ''.join(up_seq_chars)

            # convert qualities to ints
            q_scores = [phred33_to_q(c) for c in qual]

            yield (read_id, up_seq, q_scores)

def passes_quality_all_gt(q_scores, threshold):
    """Return True iff every base quality is strictly greater than threshold."""
    # manual 'all' check (avoid built-ins)
    for q in q_scores:
        if not (q > threshold):
            return False
    return True

def naive_approx_match(ref, pat, max_mismatches=2):
    """Return list of 0-based positions where pat matches ref with <= max_mismatches."""
    positions = []
    n = len(ref)
    m = len(pat)
    if m == 0 or m > n:
        return positions
    last_start = n - m
    # slide window
    i = 0
    while i <= last_start:
        mism = 0
        j = 0
        # compare pat[j] vs ref[i+j]
        while j < m:
            if ref[i + j] != pat[j]:
                mism += 1
                if mism > max_mismatches:
                    break  # early exit for this i
            j += 1
        if mism <= max_mismatches and j == m:
            positions.append(i)
        i += 1
    return positions

def main():
    p = argparse.ArgumentParser(description="Naïve matcher allowing up to two mismatches with quality filtering.")
    p.add_argument("fasta", help="Reference genome in FASTA (e.g., lambda_virus.fa)")
    p.add_argument("fastq", help="Reads in FASTQ (e.g., myReads.fastq)")
    p.add_argument("--qmin", type=int, default=5, help="Keep reads only if ALL base Q > qmin (default: 5)")
    p.add_argument("--k", type=int, default=2, help="Maximum mismatches allowed (default: 2)")
    p.add_argument("--one-based", action="store_true", help="Report positions as 1-based instead of 0-based")
    args = p.parse_args()
 
    ref_name, ref_seq = parse_fasta(args.fasta)
    ref_len = len(ref_seq)
    if ref_len == 0:
        print("ERROR: Reference sequence is empty.", file=sys.stderr)
        sys.exit(1)

    total_reads = 0
    quality_filtered_out = 0
    no_match_count = 0
    matched_reads = 0

    # Process reads one-by-one
    for read_id, read_seq, q_scores in parse_fastq(args.fastq):
        total_reads += 1

        # Check lengths consistency with quality string
        if len(read_seq) != len(q_scores):
            # Skip malformed read gracefully
            quality_filtered_out += 1
            continue

        if not passes_quality_all_gt(q_scores, args.qmin):
            quality_filtered_out += 1
            continue

        # Find approximate matches (<= k mismatches)
        positions = naive_approx_match(ref_seq, read_seq, max_mismatches=args.k)

        if positions:
            matched_reads += 1
            # convert to 1-based if requested
            if args.one_based:
                # convert without list comprehensions
                out_positions = []
                for pos in positions:
                    out_positions.append(str(pos + 1))
            else:
                out_positions = [str(pos) for pos in positions]
            print(f"{read_id}\t{','.join(out_positions)}")
        else:
            no_match_count += 1
            print(f"{read_id}\tno_match")

    # Summary
    print(f"# total_reads: {total_reads}", file=sys.stderr)
    print(f"# quality_filtered_out (Q <= {args.qmin} present): {quality_filtered_out}", file=sys.stderr)
    print(f"# matched_reads: {matched_reads}", file=sys.stderr)
    print(f"# no_match_after_filter: {no_match_count}", file=sys.stderr)

if __name__ == "__main__":
    main()