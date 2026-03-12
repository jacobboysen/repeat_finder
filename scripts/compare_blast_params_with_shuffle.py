#!/usr/bin/env python3
"""
Compare BLAST parameters with shuffled controls.

Tests combinations of:
- gapopen: 2, 5, 10
- reward/penalty: (1,-1), (1,-2)
- dust: yes, no

Uses 10% sample of real 3'UTRs and their dinucleotide-shuffled versions.
"""

import sys
import random
import subprocess
import tempfile
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple

import pandas as pd

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).parent))

from utils.paths import get_blastdb_dir, get_results_dir, get_query_fasta_path
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
        dinuc = sequence[i:i+2]
        graph[dinuc[0]].append(dinuc[1])

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
            # Write in 80-char lines
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')


def run_blast(query_fasta: Path, db_path: Path, output_file: Path,
              gapopen: int, gapextend: int, reward: int, penalty: int, dust: str) -> Tuple[bool, str]:
    """Run BLAST with specified parameters."""
    cmd = [
        'blastn',
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
        return pd.DataFrame(columns=BLAST_COLUMNS[:16])  # Exclude 'strand'

    df = pd.read_csv(tsv_path, sep='\t', header=None,
                     names=BLAST_COLUMNS[:16])
    return df


def analyze_hits(df: pd.DataFrame) -> Dict:
    """Compute statistics for BLAST hits."""
    if len(df) == 0:
        return {
            'total_hits': 0,
            'unique_queries': 0,
            'mean_length': 0,
            'mean_pident': 0,
            'mean_evalue': 0,
            'hits_with_gaps': 0,
            'pct_with_gaps': 0,
            'hits_e_lt_1': 0,
            'hits_e_lt_01': 0,
            'hits_e_lt_001': 0,
        }

    return {
        'total_hits': len(df),
        'unique_queries': df['qseqid'].nunique(),
        'mean_length': df['length'].mean(),
        'mean_pident': df['pident'].mean(),
        'mean_evalue': df['evalue'].mean(),
        'hits_with_gaps': (df['gapopen'] > 0).sum(),
        'pct_with_gaps': 100 * (df['gapopen'] > 0).mean(),
        'hits_e_lt_1': (df['evalue'] < 1).sum(),
        'hits_e_lt_01': (df['evalue'] < 0.1).sum(),
        'hits_e_lt_001': (df['evalue'] < 0.01).sum(),
    }


def main():
    # Configuration
    SAMPLE_FRACTION = 0.10
    SEED = 42

    # Parameter combinations to test
    GAPOPEN_VALUES = [2, 5, 10]
    GAPEXTEND_VALUES = [1, 2]  # Need gapextend >= 1 for valid scoring
    REWARD_PENALTY_PAIRS = [(1, -1), (1, -2)]
    DUST_VALUES = ['yes', 'no']

    # Paths
    db_path = get_blastdb_dir() / 'dmel_te_combined'
    results_dir = get_results_dir() / 'param_comparison_shuffle'
    results_dir.mkdir(parents=True, exist_ok=True)

    # Collect all 3'UTRs from multiple gene groups
    print("Loading 3'UTR sequences...")
    all_sequences = {}
    gene_groups = ['germ_plasm', 'housekeeping', 'somatic', 'cleared', 'adult']

    for group in gene_groups:
        try:
            fasta_path = get_query_fasta_path(group, strand='sense')
            seqs = parse_fasta(fasta_path)
            print(f"  {group}: {len(seqs)} sequences")
            all_sequences.update(seqs)
        except Exception as e:
            print(f"  {group}: skipped ({e})")

    print(f"\nTotal sequences: {len(all_sequences)}")

    # Sample 10%
    rng = random.Random(SEED)
    seq_ids = list(all_sequences.keys())
    sample_size = max(1, int(len(seq_ids) * SAMPLE_FRACTION))
    sampled_ids = rng.sample(seq_ids, sample_size)

    sampled_real = {sid: all_sequences[sid] for sid in sampled_ids}
    print(f"Sampled {len(sampled_real)} sequences ({SAMPLE_FRACTION*100:.0f}%)")

    # Create shuffled versions (using same seed for reproducibility)
    print("\nCreating dinucleotide-shuffled sequences...")
    shuffle_rng = random.Random(SEED + 1)  # Different seed but deterministic
    sampled_shuffled = {}
    for sid, seq in sampled_real.items():
        sampled_shuffled[f"{sid}_shuffled"] = dinucleotide_shuffle(seq, shuffle_rng)

    # Write temporary FASTA files
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        real_fasta = tmpdir / 'real_sample.fasta'
        shuffled_fasta = tmpdir / 'shuffled_sample.fasta'

        write_fasta(sampled_real, real_fasta)
        write_fasta(sampled_shuffled, shuffled_fasta)

        print(f"\nWrote {len(sampled_real)} real sequences to temp file")
        print(f"Wrote {len(sampled_shuffled)} shuffled sequences to temp file")

        # Run BLAST for each parameter combination
        results = []

        total_combos = len(GAPOPEN_VALUES) * len(GAPEXTEND_VALUES) * len(REWARD_PENALTY_PAIRS) * len(DUST_VALUES) * 2
        combo_num = 0

        for gapopen in GAPOPEN_VALUES:
            for gapextend in GAPEXTEND_VALUES:
                for reward, penalty in REWARD_PENALTY_PAIRS:
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
                                stats = analyze_hits(df)
                                stats['gapopen'] = gapopen
                                stats['gapextend'] = gapextend
                                stats['reward'] = reward
                                stats['penalty'] = penalty
                                stats['dust'] = dust
                                stats['seq_type'] = seq_type
                                stats['combo_name'] = combo_name
                                results.append(stats)
                                print(f"  -> {stats['total_hits']} hits, {stats['unique_queries']} unique queries")
                            else:
                                print(f"  -> BLAST failed: {error_msg[:100]}")

    # Create summary DataFrame
    results_df = pd.DataFrame(results)

    # Save raw results
    results_df.to_csv(results_dir / 'summary.tsv', sep='\t', index=False)
    print(f"\n\nSaved summary to {results_dir / 'summary.tsv'}")

    # Print comparison table
    print("\n" + "="*120)
    print("PARAMETER COMPARISON: REAL vs SHUFFLED")
    print("="*120)

    # Pivot for easy comparison
    for dust in DUST_VALUES:
        print(f"\n### DUST = {dust} ###\n")

        subset = results_df[results_df['dust'] == dust]

        print(f"{'Parameters':<35} {'Real Hits':>10} {'Shuf Hits':>10} {'Ratio':>8} {'Real Gaps%':>10} {'Shuf Gaps%':>10} {'Real E<0.01':>12}")
        print("-" * 105)

        for gapopen in GAPOPEN_VALUES:
            for gapextend in GAPEXTEND_VALUES:
                for reward, penalty in REWARD_PENALTY_PAIRS:
                    real_row = subset[(subset['gapopen'] == gapopen) &
                                      (subset['gapextend'] == gapextend) &
                                      (subset['reward'] == reward) &
                                      (subset['seq_type'] == 'real')]
                    shuf_row = subset[(subset['gapopen'] == gapopen) &
                                      (subset['gapextend'] == gapextend) &
                                      (subset['reward'] == reward) &
                                      (subset['seq_type'] == 'shuffled')]

                    if len(real_row) > 0 and len(shuf_row) > 0:
                        real_hits = real_row['total_hits'].values[0]
                        shuf_hits = shuf_row['total_hits'].values[0]
                        ratio = real_hits / shuf_hits if shuf_hits > 0 else float('inf')
                        real_gaps = real_row['pct_with_gaps'].values[0]
                        shuf_gaps = shuf_row['pct_with_gaps'].values[0]
                        real_strong = real_row['hits_e_lt_001'].values[0]

                        param_str = f"go={gapopen} ge={gapextend} r={reward}/p={penalty}"
                        print(f"{param_str:<35} {real_hits:>10} {shuf_hits:>10} {ratio:>8.2f} {real_gaps:>9.1f}% {shuf_gaps:>9.1f}% {real_strong:>12}")

    # Summary statistics
    print("\n" + "="*100)
    print("KEY INSIGHTS")
    print("="*100)

    # Compare dust on vs off
    dust_yes = results_df[results_df['dust'] == 'yes']
    dust_no = results_df[results_df['dust'] == 'no']

    real_dust_yes = dust_yes[dust_yes['seq_type'] == 'real']['total_hits'].mean()
    real_dust_no = dust_no[dust_no['seq_type'] == 'real']['total_hits'].mean()
    shuf_dust_yes = dust_yes[dust_yes['seq_type'] == 'shuffled']['total_hits'].mean()
    shuf_dust_no = dust_no[dust_no['seq_type'] == 'shuffled']['total_hits'].mean()

    print(f"\nDUST effect on hit counts (averaged across params):")
    print(f"  Real UTRs:     dust=no: {real_dust_no:.0f} hits  ->  dust=yes: {real_dust_yes:.0f} hits")
    print(f"  Shuffled UTRs: dust=no: {shuf_dust_no:.0f} hits  ->  dust=yes: {shuf_dust_yes:.0f} hits")

    # Best signal-to-noise ratio
    print(f"\nSignal-to-noise (Real/Shuffled ratio) by parameter combo:")

    for _, row in results_df[results_df['seq_type'] == 'real'].iterrows():
        shuf_row = results_df[(results_df['combo_name'] == row['combo_name']) &
                              (results_df['seq_type'] == 'shuffled')]
        if len(shuf_row) > 0:
            shuf_hits = shuf_row['total_hits'].values[0]
            ratio = row['total_hits'] / shuf_hits if shuf_hits > 0 else float('inf')
            print(f"  {row['combo_name']}: {ratio:.2f}x")


if __name__ == '__main__':
    main()
