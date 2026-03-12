#!/usr/bin/env python3
"""
Full parameter sweep for TE BLAST analysis.

Runs BLAST with multiple parameter combinations on a sample of genome-wide UTRs
with paired dinucleotide-shuffled controls.

Parameter space:
- word_size: 7, 11
- gapopen: 2, 5, 10
- gapextend: 1 (only with go=2), 2
- penalty: -1, -2
- dust: yes, no

Output columns (17 total):
qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend,
sstart, send, evalue, bitscore, qlen, slen, qseq, sseq, stitle

Usage:
    python scripts/full_param_sweep.py [--sample-frac 0.1] [--seed 42] [--threads 4]
"""

import argparse
import sys
import random
import subprocess
import json
from pathlib import Path
from collections import defaultdict
from datetime import datetime
from typing import Dict, List, Tuple

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).parent))

from utils.paths import get_blastdb_dir, get_results_dir, get_data_dir


def dinucleotide_shuffle(sequence: str, rng: random.Random) -> str:
    """
    Shuffle sequence preserving dinucleotide frequencies.
    Uses Altschul-Erickson algorithm with Eulerian path.
    """
    if len(sequence) < 2:
        return sequence

    sequence = sequence.upper()

    # Build graph of dinucleotide transitions
    graph = defaultdict(list)
    for i in range(len(sequence) - 1):
        graph[sequence[i]].append(sequence[i+1])

    # Shuffle edges for each node
    for node in graph:
        rng.shuffle(graph[node])

    # Find Eulerian path starting from first nucleotide
    result = [sequence[0]]
    current = sequence[0]

    while graph[current]:
        next_nuc = graph[current].pop()
        result.append(next_nuc)
        current = next_nuc

    shuffled = ''.join(result)

    # Verify length preserved
    if len(shuffled) != len(sequence):
        # Fallback to simple shuffle if Eulerian path fails
        seq_list = list(sequence)
        rng.shuffle(seq_list)
        return ''.join(seq_list)

    return shuffled


def parse_fasta(fasta_path: Path) -> Dict[str, str]:
    """Parse FASTA file into dict of id -> sequence."""
    sequences = {}
    current_id = None
    current_seq = []

    with open(fasta_path) as f:
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


def write_fasta(sequences: Dict[str, str], output_path: Path):
    """Write sequences dict to FASTA file."""
    with open(output_path, 'w') as f:
        for seq_id, seq in sequences.items():
            f.write(f">{seq_id}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')


# BLAST path
BLASTN_PATH = "/Users/jacobboysen/miniconda3/envs/bioinformatics-program/bin/blastn"

# Output format: 17 columns
BLAST_OUTFMT = '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qseq sseq stitle'


def run_blast(query_fasta: Path, db_path: Path, output_file: Path,
              word_size: int, gapopen: int, gapextend: int,
              penalty: int, reward: int, dust: str,
              threads: int = 4) -> Tuple[bool, str]:
    """Run BLAST with specified parameters."""
    cmd = [
        BLASTN_PATH,
        '-query', str(query_fasta),
        '-db', str(db_path),
        '-out', str(output_file),
        '-outfmt', BLAST_OUTFMT,
        '-word_size', str(word_size),
        '-gapopen', str(gapopen),
        '-gapextend', str(gapextend),
        '-penalty', str(penalty),
        '-reward', str(reward),
        '-dust', dust,
        '-evalue', '10',
        '-num_threads', str(threads),
        '-max_target_seqs', '500'
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        return False, result.stderr
    return True, ""


def generate_param_combos() -> List[Dict]:
    """Generate all valid parameter combinations."""
    combos = []

    for word_size in [7, 11]:
        for dust in ['yes', 'no']:
            for penalty in [-1, -2]:
                # gapextend=1 only valid with gapopen=2
                for gapopen, gapextend in [(2, 1), (2, 2), (5, 2), (10, 2)]:
                    combo = {
                        'word_size': word_size,
                        'gapopen': gapopen,
                        'gapextend': gapextend,
                        'penalty': penalty,
                        'reward': 1,
                        'dust': dust,
                    }
                    # Create readable name
                    combo['name'] = f"ws{word_size}_go{gapopen}_ge{gapextend}_p{penalty}_dust{dust}"
                    combos.append(combo)

    return combos


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sample-frac', type=float, default=0.1,
                        help='Fraction of UTRs to sample (default: 0.1 = 10%%)')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed for reproducibility')
    parser.add_argument('--threads', type=int, default=4,
                        help='BLAST threads per run')
    parser.add_argument('--output-dir', type=str, default=None,
                        help='Output directory (default: results/param_sweep_full)')

    args = parser.parse_args()

    # Paths
    utr_fasta = get_data_dir() / 'references' / 'dmel_3utr.fasta'
    db_path = get_blastdb_dir() / 'dmel_te_combined'

    if args.output_dir:
        output_dir = Path(args.output_dir)
    else:
        output_dir = get_results_dir() / 'param_sweep_full'

    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("FULL PARAMETER SWEEP")
    print(f"Sample fraction: {args.sample_frac*100:.0f}%")
    print(f"Random seed: {args.seed}")
    print(f"Output directory: {output_dir}")
    print("=" * 70)

    # Load UTRs
    print(f"\nLoading UTRs from {utr_fasta}...")
    all_sequences = parse_fasta(utr_fasta)
    print(f"  Total UTRs: {len(all_sequences):,}")

    # Sample
    rng = random.Random(args.seed)
    seq_ids = list(all_sequences.keys())
    sample_size = int(len(seq_ids) * args.sample_frac)
    sampled_ids = rng.sample(seq_ids, sample_size)

    sampled_real = {sid: all_sequences[sid] for sid in sampled_ids}
    print(f"  Sampled: {len(sampled_real):,} UTRs ({args.sample_frac*100:.0f}%)")

    # Create paired shuffled sequences
    print("\nCreating paired dinucleotide-shuffled sequences...")
    shuffle_rng = random.Random(args.seed + 1)
    sampled_shuffled = {}

    for sid, seq in sampled_real.items():
        shuffled_seq = dinucleotide_shuffle(seq, shuffle_rng)
        sampled_shuffled[f"{sid}_shuffled"] = shuffled_seq

    print(f"  Created {len(sampled_shuffled):,} paired shuffled sequences")

    # Save sample info
    sample_info = {
        'sample_frac': args.sample_frac,
        'seed': args.seed,
        'n_real': len(sampled_real),
        'n_shuffled': len(sampled_shuffled),
        'total_utrs': len(all_sequences),
        'timestamp': datetime.now().isoformat(),
    }

    # UTR length distribution for the sample
    utr_lengths = [len(seq) for seq in sampled_real.values()]
    sample_info['utr_length_stats'] = {
        'min': min(utr_lengths),
        'max': max(utr_lengths),
        'mean': sum(utr_lengths) / len(utr_lengths),
        'median': sorted(utr_lengths)[len(utr_lengths)//2],
    }

    with open(output_dir / 'sample_info.json', 'w') as f:
        json.dump(sample_info, f, indent=2)

    # Write FASTA files
    real_fasta = output_dir / 'sampled_real.fasta'
    shuffled_fasta = output_dir / 'sampled_shuffled.fasta'

    write_fasta(sampled_real, real_fasta)
    write_fasta(sampled_shuffled, shuffled_fasta)
    print(f"\nWrote FASTAs to {output_dir}")

    # Generate parameter combinations
    param_combos = generate_param_combos()
    print(f"\nParameter combinations: {len(param_combos)}")

    # Save parameter info
    with open(output_dir / 'param_combos.json', 'w') as f:
        json.dump(param_combos, f, indent=2)

    # Run BLAST for each combination
    total_runs = len(param_combos) * 2  # real + shuffled
    run_num = 0
    results_summary = []

    print("\n" + "=" * 70)
    print("RUNNING BLAST")
    print("=" * 70)

    for combo in param_combos:
        for seq_type, fasta_file in [('real', real_fasta), ('shuffled', shuffled_fasta)]:
            run_num += 1
            output_file = output_dir / f"{combo['name']}_{seq_type}.tsv"

            print(f"\n[{run_num}/{total_runs}] {combo['name']} ({seq_type})")

            # Skip if already exists
            if output_file.exists() and output_file.stat().st_size > 0:
                print(f"  -> Already exists, skipping")
                # Count lines for summary
                with open(output_file) as f:
                    n_hits = sum(1 for _ in f)
                results_summary.append({
                    'combo': combo['name'],
                    'seq_type': seq_type,
                    'n_hits': n_hits,
                    'status': 'skipped',
                    **combo
                })
                continue

            success, error_msg = run_blast(
                fasta_file, db_path, output_file,
                word_size=combo['word_size'],
                gapopen=combo['gapopen'],
                gapextend=combo['gapextend'],
                penalty=combo['penalty'],
                reward=combo['reward'],
                dust=combo['dust'],
                threads=args.threads
            )

            if success:
                # Count hits
                with open(output_file) as f:
                    n_hits = sum(1 for _ in f)
                print(f"  -> {n_hits:,} hits")
                results_summary.append({
                    'combo': combo['name'],
                    'seq_type': seq_type,
                    'n_hits': n_hits,
                    'status': 'success',
                    **combo
                })
            else:
                print(f"  -> FAILED: {error_msg[:100]}")
                results_summary.append({
                    'combo': combo['name'],
                    'seq_type': seq_type,
                    'n_hits': 0,
                    'status': 'failed',
                    'error': error_msg[:200],
                    **combo
                })

    # Save summary
    import pandas as pd
    summary_df = pd.DataFrame(results_summary)
    summary_df.to_csv(output_dir / 'sweep_summary.tsv', sep='\t', index=False)

    print("\n" + "=" * 70)
    print("SWEEP COMPLETE")
    print("=" * 70)
    print(f"\nResults saved to: {output_dir}")
    print(f"Summary: {output_dir / 'sweep_summary.tsv'}")

    # Quick summary stats
    print("\n--- Quick Summary ---")
    for dust in ['yes', 'no']:
        for ws in [7, 11]:
            subset = summary_df[(summary_df['dust'] == dust) &
                               (summary_df['word_size'] == ws) &
                               (summary_df['seq_type'] == 'real')]
            if len(subset) > 0:
                mean_hits = subset['n_hits'].mean()
                print(f"  dust={dust}, ws={ws}: {mean_hits:,.0f} avg hits (real)")

    return 0


if __name__ == '__main__':
    sys.exit(main())
