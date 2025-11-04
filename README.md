# kfilt - K-mer Based Read Filtering

A fast, multi-threaded tool for filtering FASTA/FASTQ reads based on k-mer matching using BK-trees for efficient approximate matching with support for Hamming distance tolerance.

## Features

- **Hybrid index architecture**: Bloom filter + hash table + BK-tree
- **Multiple input formats**: Paired-end, single-end, or interleaved (mixed paired/single)
- **Auto-detection**: Automatically detects FASTA or FASTQ input
- **Flexible output**: Output as FASTA or FASTQ (optionally gzip-compressed)
- **Approximate matching**: Configurable Hamming distance tolerance (0-2+)
- **Multi-threaded**: Scales across multiple CPU cores

## Installation

### Prerequisites

- Go 1.21 or later
- Linux/macOS/Windows

### Build

```bash
git clone https://github.com/yourusername/kfilt.git
cd kfilt
go mod init kfilt
go mod tidy
go build
```

## Quick Start

### Step 1: Generate unique k-mers using meryl

```bash
#Count k-mers from reference
meryl count k=31 memory=30 threads=20 reference.fa.gz output reference.meryl

#Count k-mers from alleles/variants
meryl count k=31 memory=30 threads=20 alleles.fa.gz output alleles.meryl

#Get difference (kmers in the alleles not in the reference)
meryl difference alleles.meryl reference.meryl output unique.meryl

#Export to text
meryl print unique.meryl > unique_kmers.txt
```

### Step 2: Build BK-tree index

```bash
./kfilt build -k unique_kmers.txt -o unique.bktree.idx -K 31
```

### Step 3: Extract unmapped reads (if starting from BAM/CRAM)

```bash
#Extract unmapped reads from CRAM and convert to interleaved FASTA
samtools view -u -f 4 -@ 5 -T reference.fa.gz input.cram |samtools sort -@ 5 -n | samtools fasta -0 /dev/null unmapped_sorted.bam | gzip > unmapped.fa.gz
```

### Step 4: Filter reads

```bash
# Filter interleaved FASTA, output as FASTA (compressed)
./kfilt filter \
    -I unmapped.fa.gz \
    -i unique.bktree.idx \
    -o filtered.fa.gz \
    -f fasta \
    -z \
    -v stats.tsv \
    -n 1 \
    -m 0 \
    -t 8
```

## Usage

### kfilt build

Build a BK-tree index from k-mer list.

```bash
kfilt build -k KMERS_FILE -o INDEX_FILE [-K KMER_SIZE]

Options:
-k, --kmers string Input k-mer file (meryl print output) [required]
-o, --output string Output BK-tree index file (default: bktree.idx)
-K, --kmer-size int K-mer size (default: 31)
```

### kfilt filter

Filter reads based on k-mer matching.

```bash
kfilt filter [options]

Input modes:
-1, --input1 string Input FASTA/FASTQ R1 or single-end
-2, --input2 string Input FASTA/FASTQ R2 for paired-end
-I, --interleaved string Interleaved FASTA/FASTQ (mixed paired/single)

Output options:
-o, --output string Output file [required]
-f, --output-format Output format: fasta or fastq (default: fastq)
-z, --compress Compress output with gzip
-Z, --compress-level int Gzip compression level 1-9 (default: 6)

Filtering options:
-i, --index string BK-tree index file (default: bktree.idx)
-n, --min-matches int Minimum matching k-mers required (default: 5)
-m, --hamming-dist int Maximum Hamming distance (default: 1)
-B, --both-match For pairs, require BOTH reads meet threshold

Performance:
-t, --threads int Number of threads (default: CPU count)
-v, --verbose string Verbose per-read output file (optional)
```

## Examples

### Example 1: Single-end reads

```bash
kfilt filter \
    -1 reads.fq \
    -i alleles.idx \
    -o filtered.fq \ #uncompressed
    -n 1 \ #at least 1 matching k-mer in a read 
    -m 0  #perfect matches
```

### Example 2: Paired-end reads (separate files)

```bash
kfilt filter \
    -1 R1.fq.gz \
    -2 R2.fq.gz \
    -i alleles.idx \
    -o filtered.fa.gz \
    -f fasta \ #output fasta
    -z \ #compress
    -n 5 \ #at least 5 matching k-mers in a read
    -m 1 #allows up to hamming-distance 1
```

### Example 3: Mixed interleaved format

```bash
kfilt filter \
    -I mixed_interleaved.fa.gz \ #mixed paired end and single-end in interleaved fasta
    -i alleles.idx \
    -o filtered.fa.gz \
    -f fasta \
    -z \
    -v stats.tsv \ #output per-read statistics on filtering
    -n 1 \
    -m 0
```

## Output

### Filtered reads output

By default, reads that pass the filter are written to the output file in the specified format (FASTA or FASTQ).

For paired-end reads:
- Default behavior: Keeps pair if **combined** k-mer matches ≥ threshold
- With `--both-match`: Requires **each read independently** to meet threshold

### Verbose output (TSV format)

When using `-v`, produces tab-separated output with per-read statistics:

```txt
read_id read_type total_kmers matched_kmers best_hamming_dist status
READ001 R1 240 15 0 KEPT
READ001 R2 240 8 1 KEPT
READ002 R1 240 2 1 FILTERED
READ002 R2 240 1 1 FILTERED
```

Columns:
- `read_id`: Read identifier (without @ or >)
- `read_type`: R1 (first of pair), R2 (second of pair), or SE (single-end)
- `total_kmers`: Total k-mers extracted from read (forward + reverse complement)
- `matched_kmers`: Number of k-mers matching the index
- `best_hamming_dist`: Lowest Hamming distance found (-1 if no match)
- `status`: KEPT or FILTERED

## Performance Tips

1. **Use appropriate Hamming distance**: Higher values slow down search. Start with `-m 0` for exact matching or `-m 1` for near-exact matching.

2. **Thread count**: Use `-t` to match your CPU count for best performance.

3. **Compression**: Gzip compression (`-z`) adds overhead but reduces output size for interleaved data.

4. **Large k-mer sets**: BK-tree performance depends on k-mer distribution. For very large sets (>10M k-mers), expect longer build times.

## Algorithm

### BK-Tree (Burkhard-Kelling Tree)

Metric tree data structure that enables fast approximate matching. Each node stores:
- A k-mer string
- Edges to child nodes labeled by Hamming distance

At query time, uses triangle inequality to prune search space, dramatically reducing comparisons needed for approximate matching.

## File Format Support

**Input:**
- FASTA: Headers start with `>`
- FASTQ: Headers start with `@`
- Gzip-compressed versions (`.gz`) automatically detected

**Output:**
- FASTA: `>header\nsequence\n`
- FASTQ: `@header\nsequence\n+\nquality\n`
- Both formats support gzip compression when `-z` is specified

## License

MIT License - see LICENSE file

## Contact

For issues, suggestions, or contributions, please open an issue on GitHub.
