#!/usr/bin/env python3
"""
Systematic BLAST parameter sweep for TE detection sensitivity.

Runs BLAST with multiple parameter combinations to find optimal settings
for detecting ancient/degraded TE sequences in 3'UTRs.

Parameter grid (324 combinations):
- word_size: [7, 9, 11]
- gapopen: [2, 5, 10]
- gapextend: [1, 2, 4]
- penalty: [-1, -2, -3]
- reward: [1, 2]
- task: ["blastn", "blastn-short"]
"""

import argparse
import itertools
import os
import subprocess
import sys
import yaml
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import pandas as pd


# Parameter grid - using only BLAST-supported reward/penalty combinations
# Valid combinations: reward=1 with penalty=-1,-2,-3,-4; reward=2 with penalty=-3,-5,-7
PARAM_GRID = {
    'word_size': [7, 9, 11],
    'gapopen': [2, 5, 10],
    'gapextend': [1, 2, 4],
    'reward_penalty': [(1, -1), (1, -2), (1, -3), (2, -3)],  # Valid BLAST combinations
    'task': ['blastn', 'blastn-short']
}

# Fixed parameters for sensitive TE detection
FIXED_PARAMS = {
    'evalue': 10,
    'dust': 'no',
    'soft_masking': 'false',
    'max_target_seqs': 1000,
    'max_hsps': 50,
    'outfmt': '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qseq sseq'
}

# Column names for BLAST output
BLAST_COLUMNS = [
    'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
    'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
    'qlen', 'slen', 'qseq', 'sseq'
]


def classify_strand(sstart, send):
    """
    Classify hit strand based on subject coordinates.

    BLAST reports subject strand via coordinate order:
    - Plus strand (sense): sstart < send
    - Minus strand (antisense): sstart > send

    Returns:
        str: 'plus' or 'minus'
    """
    return 'plus' if sstart < send else 'minus'


def add_strand_column(df):
    """Add strand classification column to BLAST results DataFrame."""
    if len(df) > 0:
        df['strand'] = df.apply(lambda row: classify_strand(row['sstart'], row['send']), axis=1)
    else:
        df['strand'] = []
    return df


def generate_param_combinations():
    """Generate all parameter combinations."""
    keys = list(PARAM_GRID.keys())
    values = list(PARAM_GRID.values())

    combinations = []
    for combo in itertools.product(*values):
        param_dict = dict(zip(keys, combo))

        # Unpack reward_penalty tuple into separate keys
        if 'reward_penalty' in param_dict:
            reward, penalty = param_dict.pop('reward_penalty')
            param_dict['reward'] = reward
            param_dict['penalty'] = penalty

        combinations.append(param_dict)

    return combinations


def build_blast_command(query, db, output, params, num_threads=1):
    """Build BLAST command string."""
    cmd = [
        'blastn',
        '-query', str(query),
        '-db', str(db),
        '-out', str(output),
        '-num_threads', str(num_threads)
    ]

    # Add fixed parameters
    for key, value in FIXED_PARAMS.items():
        cmd.extend([f'-{key}', str(value)])

    # Add variable parameters
    for key, value in params.items():
        cmd.extend([f'-{key}', str(value)])

    return cmd


def run_single_blast(args):
    """Run BLAST with a single parameter combination."""
    combo_id, params, query, db, output_dir, num_threads = args

    # Create combo directory
    combo_dir = output_dir / f"combo_{combo_id:04d}"
    combo_dir.mkdir(parents=True, exist_ok=True)

    # Output files
    results_file = combo_dir / 'blast_results.tsv'
    params_file = combo_dir / 'params.yaml'

    # Save parameters
    with open(params_file, 'w') as f:
        yaml.dump({**params, **FIXED_PARAMS}, f, default_flow_style=False)

    # Build and run BLAST command
    cmd = build_blast_command(query, db, results_file, params, num_threads)

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
            timeout=600  # 10 minute timeout
        )

        # Parse results
        if results_file.exists() and results_file.stat().st_size > 0:
            df = pd.read_csv(results_file, sep='\t', names=BLAST_COLUMNS)

            # Add strand classification
            df = add_strand_column(df)

            # Calculate strand breakdown
            plus_hits = (df['strand'] == 'plus').sum() if len(df) > 0 else 0
            minus_hits = (df['strand'] == 'minus').sum() if len(df) > 0 else 0

            stats = {
                'combo_id': combo_id,
                'success': True,
                'total_hits': len(df),
                'plus_strand_hits': plus_hits,
                'minus_strand_hits': minus_hits,
                'strand_ratio': plus_hits / minus_hits if minus_hits > 0 else float('inf'),
                'unique_queries': df['qseqid'].nunique(),
                'unique_subjects': df['sseqid'].nunique(),
                'mean_evalue': df['evalue'].mean() if len(df) > 0 else None,
                'median_evalue': df['evalue'].median() if len(df) > 0 else None,
                'min_evalue': df['evalue'].min() if len(df) > 0 else None,
                'mean_pident': df['pident'].mean() if len(df) > 0 else None,
                'mean_length': df['length'].mean() if len(df) > 0 else None,
                'hits_lt_1': (df['evalue'] < 1).sum(),
                'hits_lt_0.1': (df['evalue'] < 0.1).sum(),
                'hits_lt_0.01': (df['evalue'] < 0.01).sum(),
            }

            # Save updated results with strand column
            df.to_csv(results_file, sep='\t', index=False, header=False)
        else:
            stats = {
                'combo_id': combo_id,
                'success': True,
                'total_hits': 0,
                'plus_strand_hits': 0,
                'minus_strand_hits': 0,
                'strand_ratio': None,
                'unique_queries': 0,
                'unique_subjects': 0,
            }

        return {**params, **stats}

    except subprocess.TimeoutExpired:
        return {**params, 'combo_id': combo_id, 'success': False, 'error': 'timeout'}
    except subprocess.CalledProcessError as e:
        return {**params, 'combo_id': combo_id, 'success': False, 'error': str(e.stderr)}


def run_parameter_sweep(query, db, output_dir, num_threads=4, max_workers=4,
                        subset=None, verbose=False):
    """Run full parameter sweep."""
    combinations = generate_param_combinations()
    total = len(combinations)

    if subset:
        combinations = combinations[:subset]
        print(f"Running subset: {len(combinations)} of {total} combinations")
    else:
        print(f"Running {total} parameter combinations")

    # Prepare arguments for parallel execution
    args_list = [
        (i, params, query, db, output_dir, num_threads)
        for i, params in enumerate(combinations)
    ]

    results = []
    completed = 0

    # Run in parallel
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(run_single_blast, args): args[0]
                   for args in args_list}

        for future in as_completed(futures):
            combo_id = futures[future]
            try:
                result = future.result()
                results.append(result)
                completed += 1

                if verbose or completed % 10 == 0:
                    print(f"  Progress: {completed}/{len(combinations)} "
                          f"(combo {combo_id}: {result.get('total_hits', 0)} hits)")
            except Exception as e:
                print(f"  Error in combo {combo_id}: {e}")
                results.append({
                    'combo_id': combo_id,
                    'success': False,
                    'error': str(e)
                })
                completed += 1

    return results


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--query',
        type=Path,
        help='Query FASTA file (default: auto-detect based on --tier and --strand)'
    )
    parser.add_argument(
        '--tier',
        choices=['tier1', 'all'],
        default='tier1',
        help='Gene tier to analyze (default: tier1)'
    )
    parser.add_argument(
        '--strand',
        choices=['sense', 'antisense', 'shuffled'],
        default='sense',
        help='Sequence strand/type (default: sense)'
    )
    parser.add_argument(
        '--database',
        type=Path,
        default=Path('data/blastdb/dmel_te_flybase'),
        help='BLAST database (default: data/blastdb/dmel_te_flybase)'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('results/parameter_sweep'),
        help='Output directory (default: results/parameter_sweep)'
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=4,
        help='Threads per BLAST run (default: 4)'
    )
    parser.add_argument(
        '--workers',
        type=int,
        default=2,
        help='Parallel workers (default: 2)'
    )
    parser.add_argument(
        '--subset',
        type=int,
        help='Run only first N combinations (for testing)'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )
    parser.add_argument(
        '--list-params',
        action='store_true',
        help='List parameter grid and exit'
    )

    args = parser.parse_args()

    # List parameters if requested
    if args.list_params:
        print("Parameter Grid:")
        print("=" * 60)
        for key, values in PARAM_GRID.items():
            if key == 'reward_penalty':
                print(f"  reward/penalty: {values}")
            else:
                print(f"  {key}: {values}")
        combos = generate_param_combinations()
        print(f"\nTotal combinations: {len(combos)}")
        print("\nFixed Parameters:")
        for key, value in FIXED_PARAMS.items():
            print(f"  {key}: {value}")
        return 0

    # Auto-detect query file
    if args.query is None:
        queries_dir = Path('data/queries/germ_plasm')
        if args.strand == 'sense':
            if args.tier == 'tier1':
                args.query = queries_dir / '3UTR_sense_tier1.fasta'
            else:
                args.query = queries_dir / '3UTR_sense.fasta'
        elif args.strand == 'antisense':
            if args.tier == 'tier1':
                args.query = queries_dir / '3UTR_antisense_tier1.fasta'
            else:
                args.query = queries_dir / '3UTR_antisense.fasta'
        else:  # shuffled
            if args.tier == 'tier1':
                args.query = queries_dir / '3UTR_shuffled_tier1.fasta'
            else:
                args.query = queries_dir / '3UTR_shuffled.fasta'

    # Validate inputs
    if not args.query.exists():
        print(f"Error: Query file not found: {args.query}", file=sys.stderr)
        return 1

    # Check database exists
    db_file = Path(str(args.database) + '.nhr')
    if not db_file.exists():
        print(f"Error: Database not found: {args.database}", file=sys.stderr)
        return 1

    # Create timestamped output directory
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    run_dir = args.output_dir / f"{timestamp}_{args.tier}_{args.strand}"
    run_dir.mkdir(parents=True, exist_ok=True)

    print("TE Parameter Sweep")
    print("=" * 60)
    print(f"Query: {args.query}")
    print(f"Database: {args.database}")
    print(f"Output: {run_dir}")
    print(f"Threads per BLAST: {args.threads}")
    print(f"Parallel workers: {args.workers}")
    print()

    # Save sweep configuration
    config_file = run_dir / 'sweep_config.yaml'
    config = {
        'query': str(args.query),
        'database': str(args.database),
        'tier': args.tier,
        'strand': args.strand,
        'timestamp': timestamp,
        'param_grid': PARAM_GRID,
        'fixed_params': FIXED_PARAMS,
        'threads': args.threads,
        'workers': args.workers
    }
    with open(config_file, 'w') as f:
        yaml.dump(config, f, default_flow_style=False)

    # Run sweep
    results = run_parameter_sweep(
        args.query,
        args.database,
        run_dir,
        num_threads=args.threads,
        max_workers=args.workers,
        subset=args.subset,
        verbose=args.verbose
    )

    # Create summary dataframe
    df = pd.DataFrame(results)

    # Add parameter columns to front (use expanded param names)
    param_cols = ['word_size', 'gapopen', 'gapextend', 'reward', 'penalty', 'task']
    stat_cols = [c for c in df.columns if c not in param_cols]
    df = df[[c for c in param_cols if c in df.columns] + stat_cols]

    # Save summary
    summary_file = run_dir / 'sweep_summary.tsv'
    df.to_csv(summary_file, sep='\t', index=False)
    print(f"\nSaved summary to: {summary_file}")

    # Save parameter combinations list
    combos_file = run_dir / 'parameter_combinations.tsv'
    combo_df = pd.DataFrame(generate_param_combinations())
    combo_df.index.name = 'combo_id'
    combo_df.to_csv(combos_file, sep='\t')

    # Print summary statistics
    print("\n" + "=" * 60)
    print("Sweep Summary")
    print("-" * 60)

    successful = df[df['success'] == True] if 'success' in df.columns else df
    print(f"Successful runs: {len(successful)}/{len(df)}")

    if 'total_hits' in successful.columns and len(successful) > 0:
        print(f"\nHit Statistics:")
        print(f"  Total hits range: {successful['total_hits'].min():.0f} - {successful['total_hits'].max():.0f}")
        print(f"  Mean hits: {successful['total_hits'].mean():.1f}")

        # Strand breakdown
        if 'plus_strand_hits' in successful.columns:
            total_plus = successful['plus_strand_hits'].sum()
            total_minus = successful['minus_strand_hits'].sum()
            print(f"\nStrand Breakdown (across all runs):")
            print(f"  Plus strand (sense) hits:  {total_plus:.0f}")
            print(f"  Minus strand (antisense) hits: {total_minus:.0f}")
            if total_minus > 0:
                print(f"  Overall strand ratio (+/-): {total_plus/total_minus:.2f}")

        # Top 5 parameter combinations by hits
        top5 = successful.nlargest(5, 'total_hits')
        print(f"\nTop 5 combinations by hit count:")
        for _, row in top5.iterrows():
            strand_info = ""
            if 'plus_strand_hits' in row:
                strand_info = f", +:{row['plus_strand_hits']:.0f}/-:{row['minus_strand_hits']:.0f}"
            print(f"  combo_{row['combo_id']:04d}: {row['total_hits']:.0f} hits "
                  f"(ws={row['word_size']}, task={row['task']}{strand_info})")

    print("\n" + "=" * 60)
    print("Sweep complete!")

    return 0


if __name__ == '__main__':
    sys.exit(main())
