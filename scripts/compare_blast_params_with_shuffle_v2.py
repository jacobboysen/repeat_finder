#!/usr/bin/env python3
"""
Compare BLAST parameters with paired shuffled controls.

Properly samples from genome-wide 3'UTRs (30k sequences) and creates
paired dinucleotide-shuffled versions for each sampled sequence.

Tests combinations of:
- gapopen: 2, 5, 10
- reward/penalty: (1,-1), (1,-2)
- dust: yes, no
"""

import sys
import random
import subprocess
import tempfile
from pathlib import Path
from collections import defaultdict
from typing import Dict, Tuple

import pandas as pd
import numpy as np

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).parent))

from utils.paths import get_blastdb_dir, get_results_dir
from utils.blast_io import BLAST_COLUMNS


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


BLASTN_PATH = "/Users/jacobboysen/miniconda3/envs/bioinformatics-program/bin/blastn"


def run_blast(query_fasta: Path, db_path: Path, output_file: Path,
              gapopen: int, gapextend: int, reward: int, penalty: int,
              dust: str) -> Tuple[bool, str]:
    """Run BLAST with specified parameters."""
    cmd = [
        BLASTN_PATH,
        '-query', str(query_fasta),
        '-db', str(db_path),
        '-out', str(output_file),
        '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qseq sseq',
        '-word_size', '7',
        '-gapopen', str(gapopen),
        '-gapextend', str(gapextend),
        '-penalty', str(penalty),
        '-reward', str(reward),
        '-dust', dust,
        '-evalue', '10',
        '-num_threads', '4',
        '-max_target_seqs', '500'
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        return False, result.stderr
    return True, ""


def load_blast_results(tsv_path: Path) -> pd.DataFrame:
    """Load BLAST results from TSV."""
    if not tsv_path.exists() or tsv_path.stat().st_size == 0:
        return pd.DataFrame(columns=BLAST_COLUMNS[:16])

    df = pd.read_csv(tsv_path, sep='\t', header=None, names=BLAST_COLUMNS[:16])
    return df


def analyze_hits(df: pd.DataFrame, label: str = "") -> Dict:
    """Compute statistics for BLAST hits."""
    if len(df) == 0:
        return {
            'label': label,
            'total_hits': 0,
            'unique_queries': 0,
            'queries_with_hits': 0,
            'mean_length': 0,
            'median_length': 0,
            'mean_pident': 0,
            'mean_evalue': 0,
            'hits_with_gaps': 0,
            'pct_with_gaps': 0,
            'hits_e_lt_1': 0,
            'hits_e_lt_01': 0,
            'hits_e_lt_001': 0,
            'length_p25': 0,
            'length_p75': 0,
        }

    return {
        'label': label,
        'total_hits': len(df),
        'unique_queries': df['qseqid'].nunique(),
        'queries_with_hits': df['qseqid'].nunique(),
        'mean_length': df['length'].mean(),
        'median_length': df['length'].median(),
        'mean_pident': df['pident'].mean(),
        'mean_evalue': df['evalue'].mean(),
        'hits_with_gaps': (df['gapopen'] > 0).sum(),
        'pct_with_gaps': 100 * (df['gapopen'] > 0).mean(),
        'hits_e_lt_1': (df['evalue'] < 1).sum(),
        'hits_e_lt_01': (df['evalue'] < 0.1).sum(),
        'hits_e_lt_001': (df['evalue'] < 0.01).sum(),
        'length_p25': df['length'].quantile(0.25),
        'length_p75': df['length'].quantile(0.75),
    }


def main():
    # Configuration
    SAMPLE_FRACTION = 0.01  # 1% = ~300 sequences for faster testing
    SEED = 42

    # Source: genome-wide 3'UTRs
    GENOME_WIDE_3UTR = Path("/Users/jacobboysen/git_repos/repeat_finder/data/references/dmel_3utr.fasta")

    # Valid parameter combinations (gapextend=1 only works with gapopen=2)
    PARAM_COMBOS = [
        # gapopen, gapextend, reward, penalty
        (2, 1, 1, -1),
        (2, 1, 1, -2),
        (2, 2, 1, -1),
        (2, 2, 1, -2),
        (5, 2, 1, -1),
        (5, 2, 1, -2),
        (10, 2, 1, -1),
        (10, 2, 1, -2),
    ]
    DUST_VALUES = ['yes', 'no']

    # Paths
    db_path = get_blastdb_dir() / 'dmel_te_combined'
    results_dir = get_results_dir() / 'param_comparison_shuffle_v2'
    results_dir.mkdir(parents=True, exist_ok=True)

    # Load genome-wide 3'UTRs
    print(f"Loading genome-wide 3'UTRs from {GENOME_WIDE_3UTR}...")
    all_sequences = parse_fasta(GENOME_WIDE_3UTR)
    print(f"Total sequences: {len(all_sequences)}")

    # Sample 10%
    rng = random.Random(SEED)
    seq_ids = list(all_sequences.keys())
    sample_size = int(len(seq_ids) * SAMPLE_FRACTION)
    sampled_ids = rng.sample(seq_ids, sample_size)

    sampled_real = {sid: all_sequences[sid] for sid in sampled_ids}
    print(f"Sampled {len(sampled_real)} sequences ({SAMPLE_FRACTION*100:.0f}%)")

    # Create PAIRED shuffled versions (same seed for reproducibility)
    print("\nCreating paired dinucleotide-shuffled sequences...")
    shuffle_rng = random.Random(SEED + 1)
    sampled_shuffled = {}

    for sid, seq in sampled_real.items():
        # Each real sequence gets its own shuffled partner
        shuffled_seq = dinucleotide_shuffle(seq, shuffle_rng)
        sampled_shuffled[f"{sid}_shuffled"] = shuffled_seq

    print(f"Created {len(sampled_shuffled)} paired shuffled sequences")

    # Verify pairing: check a few sequences
    print("\nVerifying shuffle pairing (first 3 sequences):")
    for i, sid in enumerate(list(sampled_real.keys())[:3]):
        real_seq = sampled_real[sid]
        shuf_seq = sampled_shuffled[f"{sid}_shuffled"]
        real_gc = (real_seq.count('G') + real_seq.count('C')) / len(real_seq) * 100
        shuf_gc = (shuf_seq.count('G') + shuf_seq.count('C')) / len(shuf_seq) * 100
        print(f"  {sid}: len={len(real_seq)}, real_GC={real_gc:.1f}%, shuf_GC={shuf_gc:.1f}%")

    # Write FASTA files (persistent for inspection)
    real_fasta = results_dir / 'sampled_real.fasta'
    shuffled_fasta = results_dir / 'sampled_shuffled.fasta'

    write_fasta(sampled_real, real_fasta)
    write_fasta(sampled_shuffled, shuffled_fasta)
    print(f"\nWrote FASTAs to {results_dir}")

    # Run BLAST for each parameter combination
    results = []

    total_combos = len(PARAM_COMBOS) * len(DUST_VALUES) * 2
    combo_num = 0

    for gapopen, gapextend, reward, penalty in PARAM_COMBOS:
        for dust in DUST_VALUES:
            combo_name = f"go{gapopen}_ge{gapextend}_r{reward}p{penalty}_dust{dust}"

            for seq_type, fasta_file in [('real', real_fasta), ('shuffled', shuffled_fasta)]:
                combo_num += 1
                print(f"\n[{combo_num}/{total_combos}] Running: {combo_name} ({seq_type})")

                output_file = results_dir / f"{combo_name}_{seq_type}.tsv"

                success, error_msg = run_blast(
                    fasta_file, db_path, output_file,
                    gapopen=gapopen, gapextend=gapextend,
                    reward=reward, penalty=penalty, dust=dust
                )

                if success:
                    df = load_blast_results(output_file)
                    stats = analyze_hits(df, f"{combo_name}_{seq_type}")
                    stats['gapopen'] = gapopen
                    stats['gapextend'] = gapextend
                    stats['reward'] = reward
                    stats['penalty'] = penalty
                    stats['dust'] = dust
                    stats['seq_type'] = seq_type
                    stats['combo_name'] = combo_name
                    stats['n_input_seqs'] = len(sampled_real)
                    results.append(stats)
                    pct_with_hits = 100 * stats['queries_with_hits'] / len(sampled_real)
                    print(f"  -> {stats['total_hits']} hits from {stats['queries_with_hits']}/{len(sampled_real)} seqs ({pct_with_hits:.1f}%)")
                else:
                    print(f"  -> BLAST failed: {error_msg[:100]}")

    # Create summary DataFrame
    results_df = pd.DataFrame(results)
    results_df.to_csv(results_dir / 'summary.tsv', sep='\t', index=False)
    print(f"\n\nSaved summary to {results_dir / 'summary.tsv'}")

    # Print detailed comparison
    print("\n" + "="*130)
    print("PARAMETER COMPARISON: REAL vs PAIRED SHUFFLED CONTROLS")
    print(f"Sample size: {len(sampled_real)} sequences (10% of {len(all_sequences)} genome-wide 3'UTRs)")
    print("="*130)

    for dust in DUST_VALUES:
        print(f"\n### DUST = {dust} ###\n")

        subset = results_df[results_df['dust'] == dust]

        header = f"{'Parameters':<30} {'Real':>8} {'Shuf':>8} {'Ratio':>7} {'R/S Seqs':>10} {'RealLen':>8} {'ShufLen':>8} {'R Gaps%':>8} {'S Gaps%':>8}"
        print(header)
        print("-" * len(header))

        for gapopen, gapextend, reward, penalty in PARAM_COMBOS:
            real_row = subset[(subset['gapopen'] == gapopen) &
                              (subset['gapextend'] == gapextend) &
                              (subset['reward'] == reward) &
                              (subset['seq_type'] == 'real')]
            shuf_row = subset[(subset['gapopen'] == gapopen) &
                              (subset['gapextend'] == gapextend) &
                              (subset['reward'] == reward) &
                              (subset['seq_type'] == 'shuffled')]

            if len(real_row) > 0 and len(shuf_row) > 0:
                r = real_row.iloc[0]
                s = shuf_row.iloc[0]

                ratio = r['total_hits'] / s['total_hits'] if s['total_hits'] > 0 else float('inf')
                seq_ratio = f"{r['queries_with_hits']}/{s['queries_with_hits']}"

                param_str = f"go={gapopen} ge={gapextend} r={reward}/p={penalty}"
                print(f"{param_str:<30} {r['total_hits']:>8} {s['total_hits']:>8} {ratio:>7.2f} {seq_ratio:>10} {r['mean_length']:>8.1f} {s['mean_length']:>8.1f} {r['pct_with_gaps']:>7.1f}% {s['pct_with_gaps']:>7.1f}%")

    # Deeper analysis: what's different between real and shuffled?
    print("\n" + "="*130)
    print("DUST EFFECT ANALYSIS")
    print("="*130)

    for gapopen, gapextend, reward, penalty in [(2, 1, 1, -1)]:  # Focus on main params
        print(f"\nParameters: gapopen={gapopen}, gapextend={gapextend}, reward={reward}, penalty={penalty}")

        for seq_type in ['real', 'shuffled']:
            dust_no = results_df[(results_df['gapopen'] == gapopen) &
                                  (results_df['gapextend'] == gapextend) &
                                  (results_df['reward'] == reward) &
                                  (results_df['dust'] == 'no') &
                                  (results_df['seq_type'] == seq_type)]
            dust_yes = results_df[(results_df['gapopen'] == gapopen) &
                                   (results_df['gapextend'] == gapextend) &
                                   (results_df['reward'] == reward) &
                                   (results_df['dust'] == 'yes') &
                                   (results_df['seq_type'] == seq_type)]

            if len(dust_no) > 0 and len(dust_yes) > 0:
                hits_no = dust_no.iloc[0]['total_hits']
                hits_yes = dust_yes.iloc[0]['total_hits']
                removed = hits_no - hits_yes
                pct_removed = 100 * removed / hits_no if hits_no > 0 else 0
                print(f"  {seq_type.upper()}: dust=no: {hits_no} hits -> dust=yes: {hits_yes} hits (removed {removed}, {pct_removed:.1f}%)")


if __name__ == '__main__':
    main()
