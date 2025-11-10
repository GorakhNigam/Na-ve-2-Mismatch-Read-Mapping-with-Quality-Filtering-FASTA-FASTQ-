# Naive<=2-Mismatch-Read-Mapping-with-Quality-Filtering-FASTA-FASTQ-
A pure-Python mapper that parses FASTA/FASTQ, filters reads whose all base qualities > 5, and reports all reference start positions allowing up to 2 mismatches.
This demonstrates fundamentals of read mapping without indexes—clean I/O, Phred+33 handling, and a tight early-break mismatch loop.

**Features**
	•	Manual FASTA & FASTQ parsing (no external libraries).
	•	Phred+33 quality decoding; strict all-bases > qmin filter (default qmin=5).
	•	Naïve approximate matching with ≤2 mismatches (early exit when mismatches > k).
	•	0-based or 1-based coordinates via flag.
	•	Summary stats to stderr (totals, filtered-out reads, matches, no-matches).

**Data**
	•	Reference: lambda_virus.fa (FASTA).
	•	Reads: myReads.fastq (FASTQ with qualities).

**How it works**
	•	Load reference into a single uppercase string.
	•	Stream FASTQ in 4-line chunks; convert quality chars (ord(ch)-33) and filter reads.
	•	Slide each passing read across the reference; count mismatches character-by-character; record positions with mismatches ≤ k.


