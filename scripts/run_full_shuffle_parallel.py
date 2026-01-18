#!/usr/bin/env python3
"""
Run full genome-wide 3'UTR shuffled control analysis with 10 parallel replicates.

IMPORTANT: Output is NOT deduplicated. Raw BLAST results with aligned sequences.

This script:
1. Loads ALL 3'UTR sequences (full genome, not sampled)
2. Creates 10 dinucleotide-shuffled versions in parallel
3. Runs BLAST on each with full output including qseq/sseq
4. Saves all results for downstream comparison

Output format: 16 columns
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qseq sseq
"""

import subprocess
import sys
import random
import time
import os
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime

# Configuration
N_REPLICATES = 10
BASE_SEED = 12345
BLAST_THREADS_PER_JOB = 2  # 2 threads × 10 jobs = 20 threads max
INPUT_FASTA = Path('data/references/dmel_3utr.fasta')
TE_DB = Path('data/blastdb/dmel_te_flybase')
OUTPUT_DIR = Path('results/shuffled_full')

# BLAST output format - 16 columns WITH sequences
BLAST_OUTFMT = '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qseq sseq'


def count_dinucs(seq):
    """Count dinucleotide frequencies."""
    counts = defaultdict(int)
    seq = str(seq).upper()
    for i in range(len(seq) - 1):
        dinuc = seq[i:i+2]
        if 'N' not in dinuc and all(c in 'ACGT' for c in dinuc):
            counts[dinuc] += 1
    return dict(counts)


def dinucleotide_shuffle(sequence, rng):
    """Shuffle sequence preserving dinucleotide frequencies using Euler path."""
    seq = str(sequence).upper()

    if len(seq) < 2:
        return seq

    # Build edge list
    edges = []
    for i in range(len(seq) - 1):
        if seq[i] in 'ACGT' and seq[i+1] in 'ACGT':
            edges.append((seq[i], seq[i+1]))

    if not edges:
        return seq

    rng.shuffle(edges)

    # Build adjacency list
    graph = defaultdict(list)
    for src, dst in edges:
        graph[src].append(dst)

    # Find Eulerian path
    start = seq[0] if seq[0] in graph and graph[seq[0]] else next((n for n in graph if graph[n]), seq[0])

    stack = [start]
    path = []

    while stack:
        current = stack[-1]
        if graph[current]:
            stack.append(graph[current].pop())
        else:
            path.append(stack.pop())

    return ''.join(reversed(path))


def parse_fasta(filepath):
    """Simple FASTA parser."""
    sequences = {}
    current_id = None
    current_seq = []

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            sequences[current_id] = ''.join(current_seq)

    return sequences


def write_fasta(sequences, filepath):
    """Write sequences to FASTA file."""
    with open(filepath, 'w') as f:
        for seq_id, seq in sequences.items():
            f.write(f">{seq_id}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')


def run_single_replicate(args):
    """Run a single shuffle replicate. Called in parallel."""
    rep_num, sequences, output_dir, te_db, seed = args

    start_time = time.time()

    # Create shuffled sequences
    rng = random.Random(seed)
    shuffled_seqs = {}
    for seq_id, seq in sequences.items():
        shuffled_seqs[f"{seq_id}_shuf{rep_num}"] = dinucleotide_shuffle(seq, rng)

    # Write FASTA
    fasta_path = output_dir / f'replicate_{rep_num:02d}.fasta'
    write_fasta(shuffled_seqs, fasta_path)

    # Run BLAST
    blast_output = output_dir / f'replicate_{rep_num:02d}_blast.tsv'

    cmd = [
        'blastn',
        '-query', str(fasta_path),
        '-db', str(te_db),
        '-out', str(blast_output),
        '-outfmt', BLAST_OUTFMT,
        '-word_size', '7',
        '-gapopen', '2',
        '-gapextend', '1',
        '-penalty', '-1',
        '-reward', '1',
        '-dust', 'yes',
        '-evalue', '10',
        '-num_threads', str(BLAST_THREADS_PER_JOB),
        '-max_target_seqs', '500'
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    elapsed = time.time() - start_time

    # Count hits
    hit_count = 0
    if blast_output.exists():
        with open(blast_output) as f:
            hit_count = sum(1 for _ in f)

    # Get file size
    fasta_size = fasta_path.stat().st_size / (1024*1024)  # MB
    blast_size = blast_output.stat().st_size / (1024*1024) if blast_output.exists() else 0

    return {
        'replicate': rep_num,
        'hits': hit_count,
        'fasta_mb': fasta_size,
        'blast_mb': blast_size,
        'elapsed_sec': elapsed,
        'success': result.returncode == 0,
        'error': result.stderr if result.returncode != 0 else None
    }


def main():
    print("=" * 70)
    print("FULL GENOME-WIDE 3'UTR SHUFFLED CONTROL ANALYSIS")
    print("=" * 70)
    print(f"\nStarted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"\nConfiguration:")
    print(f"  Input: {INPUT_FASTA}")
    print(f"  TE database: {TE_DB}")
    print(f"  Replicates: {N_REPLICATES}")
    print(f"  Base seed: {BASE_SEED}")
    print(f"  Threads per BLAST: {BLAST_THREADS_PER_JOB}")
    print(f"  Output: {OUTPUT_DIR}")
    print(f"\n  *** OUTPUT IS NOT DEDUPLICATED ***")
    print(f"  *** BLAST output includes qseq/sseq (16 columns) ***")

    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Save questions/notes for later
    questions_file = OUTPUT_DIR / 'QUESTIONS_FOR_LATER.md'
    with open(questions_file, 'w') as f:
        f.write("# Questions to Address After Shuffle Analysis\n\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write("## Saved Questions\n\n")
        f.write("1. Is 10 replicates enough? Need to check saturation after runs complete.\n")
        f.write("2. Should we deduplicate before or after comparison?\n")
        f.write("3. How to best compare alignment patterns (qseq/sseq) between real and shuffled?\n")
        f.write("4. Storage format for long-term: keep TSV or convert to database?\n")
        f.write("5. What statistical test for comparing hit distributions?\n\n")
        f.write("## Notes\n\n")
        f.write("- Using dinucleotide shuffling (preserves GC and dinuc frequencies)\n")
        f.write("- Each replicate has independent seed (BASE_SEED + replicate_number)\n")
        f.write("- FASTA files are kept (not deleted) for reproducibility\n")

    # Write metadata
    metadata_file = OUTPUT_DIR / 'METADATA.md'
    with open(metadata_file, 'w') as f:
        f.write("# Shuffled Control Analysis Metadata\n\n")
        f.write(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write("## Configuration\n\n")
        f.write(f"- Input FASTA: `{INPUT_FASTA}`\n")
        f.write(f"- TE Database: `{TE_DB}`\n")
        f.write(f"- Number of replicates: {N_REPLICATES}\n")
        f.write(f"- Base random seed: {BASE_SEED}\n")
        f.write(f"- BLAST threads per job: {BLAST_THREADS_PER_JOB}\n\n")
        f.write("## BLAST Parameters\n\n")
        f.write("```\n")
        f.write("word_size=7, gapopen=2, gapextend=1, penalty=-1, reward=1\n")
        f.write("dust=yes, evalue=10, max_target_seqs=500\n")
        f.write("```\n\n")
        f.write("## Output Format\n\n")
        f.write("**16 columns (BLAST outfmt 6):**\n")
        f.write("```\n")
        f.write("qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qseq sseq\n")
        f.write("```\n\n")
        f.write("## ⚠️ IMPORTANT NOTES\n\n")
        f.write("1. **NOT DEDUPLICATED** - Raw BLAST output, may contain duplicate hits\n")
        f.write("2. **Includes aligned sequences** - qseq and sseq columns for pattern comparison\n")
        f.write("3. **Full genome** - All 30,324 transcripts, not sampled\n")
        f.write("4. **Shuffling method** - Dinucleotide shuffle (preserves dinuc frequencies)\n\n")
        f.write("## Files\n\n")
        f.write("- `replicate_NN.fasta` - Shuffled sequences for replicate NN\n")
        f.write("- `replicate_NN_blast.tsv` - BLAST results for replicate NN\n")
        f.write("- `run_summary.tsv` - Per-replicate statistics\n")
        f.write("- `METADATA.md` - This file\n")
        f.write("- `QUESTIONS_FOR_LATER.md` - Saved questions for post-analysis\n")

    # Load sequences
    print("\n" + "-" * 50)
    print("Loading sequences...")
    sequences = parse_fasta(INPUT_FASTA)
    print(f"  Loaded {len(sequences):,} sequences")

    total_bp = sum(len(s) for s in sequences.values())
    print(f"  Total base pairs: {total_bp:,}")

    # Prepare arguments for parallel execution
    print("\n" + "-" * 50)
    print(f"Starting {N_REPLICATES} parallel shuffle+BLAST jobs...")
    print(f"  (This may take 30-60 minutes per replicate)")

    job_args = [
        (i, sequences, OUTPUT_DIR, TE_DB, BASE_SEED + i)
        for i in range(1, N_REPLICATES + 1)
    ]

    # Run in parallel
    results = []
    start_total = time.time()

    # Use up to N_REPLICATES workers (each job uses BLAST_THREADS_PER_JOB threads)
    max_workers = min(N_REPLICATES, os.cpu_count() or 4)
    print(f"  Using {max_workers} parallel workers")

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(run_single_replicate, args): args[0] for args in job_args}

        for future in as_completed(futures):
            rep_num = futures[future]
            try:
                result = future.result()
                results.append(result)
                status = "✓" if result['success'] else "✗"
                print(f"  {status} Replicate {result['replicate']:2d}: {result['hits']:,} hits, "
                      f"{result['blast_mb']:.1f}MB, {result['elapsed_sec']/60:.1f}min")
            except Exception as e:
                print(f"  ✗ Replicate {rep_num}: FAILED - {e}")
                results.append({
                    'replicate': rep_num,
                    'hits': 0,
                    'fasta_mb': 0,
                    'blast_mb': 0,
                    'elapsed_sec': 0,
                    'success': False,
                    'error': str(e)
                })

    total_elapsed = time.time() - start_total

    # Sort results by replicate number
    results.sort(key=lambda x: x['replicate'])

    # Write summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    summary_file = OUTPUT_DIR / 'run_summary.tsv'
    with open(summary_file, 'w') as f:
        f.write("replicate\thits\tfasta_mb\tblast_mb\telapsed_sec\tsuccess\n")
        for r in results:
            f.write(f"{r['replicate']}\t{r['hits']}\t{r['fasta_mb']:.2f}\t{r['blast_mb']:.2f}\t{r['elapsed_sec']:.1f}\t{r['success']}\n")

    total_hits = sum(r['hits'] for r in results)
    total_storage = sum(r['blast_mb'] + r['fasta_mb'] for r in results)
    successful = sum(1 for r in results if r['success'])

    print(f"\nCompleted: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Total time: {total_elapsed/60:.1f} minutes")
    print(f"Successful replicates: {successful}/{N_REPLICATES}")
    print(f"Total hits across all replicates: {total_hits:,}")
    print(f"Total storage: {total_storage:.1f} MB")
    print(f"\nResults saved to: {OUTPUT_DIR}")
    print(f"Summary: {summary_file}")

    # Update metadata with results
    with open(metadata_file, 'a') as f:
        f.write("\n## Run Results\n\n")
        f.write(f"- Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"- Total time: {total_elapsed/60:.1f} minutes\n")
        f.write(f"- Successful replicates: {successful}/{N_REPLICATES}\n")
        f.write(f"- Total hits: {total_hits:,}\n")
        f.write(f"- Total storage: {total_storage:.1f} MB\n")

    print("\n" + "=" * 70)
    print("NEXT STEPS")
    print("=" * 70)
    print("\n1. Check saturation: Do hits per replicate stabilize?")
    print("2. Compare to real hits: results/3utr_deduplicated/3utr_deduplicated_hits.tsv")
    print("3. Analyze TE nucleotide patterns using qseq/sseq columns")
    print(f"4. Review saved questions: {questions_file}")

    return 0 if successful == N_REPLICATES else 1


if __name__ == '__main__':
    sys.exit(main())
