#!/usr/bin/env python3
"""
RNA secondary structure analysis of TE hits vs non-TE regions.

Compares:
1. Structural properties of TE-hit regions vs non-TE regions
2. Whether specific motifs fall in structured vs unstructured regions
3. Real vs shuffled sequences to test for non-random structure

Uses ViennaRNA RNAfold for structure prediction.
"""

import argparse
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from collections import defaultdict
import re

import numpy as np
import pandas as pd
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_project_root, get_results_dir
from utils.blast_io import iter_blast_results


def run_rnafold(sequences: dict, timeout_per_seq: int = 30) -> dict:
    """
    Run RNAfold on a set of sequences.

    Args:
        sequences: {name: sequence}
        timeout_per_seq: timeout in seconds per sequence

    Returns:
        {name: {mfe, structure, ensemble_diversity, ...}}
    """
    results = {}

    # Create temp FASTA file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        fasta_path = Path(f.name)
        for name, seq in sequences.items():
            # Clean sequence - only ACGU
            clean_seq = ''.join(c for c in seq.upper() if c in 'ACGTU')
            clean_seq = clean_seq.replace('T', 'U')  # Convert to RNA
            if len(clean_seq) >= 10:  # Minimum length
                f.write(f">{name}\n{clean_seq}\n")

    try:
        # Run RNAfold with partition function for ensemble diversity
        cmd = ['RNAfold', '-p', '--noPS', '-i', str(fasta_path)]
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout_per_seq * len(sequences)
        )

        # Parse output
        lines = result.stdout.strip().split('\n')
        i = 0
        while i < len(lines):
            line = lines[i]
            if line.startswith('>'):
                name = line[1:].strip()
                seq_line = lines[i+1] if i+1 < len(lines) else ''
                struct_line = lines[i+2] if i+2 < len(lines) else ''

                # Parse MFE structure line: "...(((...)))... (-5.60)"
                mfe_match = re.search(r'\(([+-]?\d+\.?\d*)\)', struct_line)
                if mfe_match:
                    mfe = float(mfe_match.group(1))
                    structure = struct_line.split()[0] if struct_line else ''

                    # Count base pairs
                    n_paired = structure.count('(')
                    n_unpaired = structure.count('.')
                    total = len(structure)

                    results[name] = {
                        'length': len(seq_line),
                        'mfe': mfe,
                        'mfe_per_nt': mfe / len(seq_line) if len(seq_line) > 0 else 0,
                        'structure': structure,
                        'n_paired': n_paired,
                        'n_unpaired': n_unpaired,
                        'pct_paired': 100 * n_paired / total if total > 0 else 0,
                    }
                i += 3
            else:
                i += 1

    except subprocess.TimeoutExpired:
        print("    RNAfold timed out")
    except Exception as e:
        print(f"    RNAfold error: {e}")
    finally:
        fasta_path.unlink(missing_ok=True)
        # Clean up RNAfold output files
        for suffix in ['_ss.ps', '_dp.ps']:
            for f in Path('.').glob(f'*{suffix}'):
                f.unlink(missing_ok=True)

    return results


def extract_te_regions(blast_path: Path, utr_seqs: dict, buffer: int = 20) -> tuple:
    """
    Extract TE-hit regions and non-TE regions from UTRs.

    Returns:
        (te_regions, non_te_regions) - both are {name: sequence}
    """
    # Build TE intervals per transcript
    te_intervals = defaultdict(list)

    for hit in iter_blast_results(blast_path):
        qseqid = hit.get('qseqid', '')
        qstart = hit.get('qstart', 0)
        qend = hit.get('qend', 0)

        if qseqid and qstart and qend:
            start = min(qstart, qend) - 1  # 0-based
            end = max(qstart, qend)
            te_intervals[qseqid].append((start, end))

    # Merge overlapping intervals with buffer
    for transcript_id in te_intervals:
        intervals = sorted(te_intervals[transcript_id])
        merged = []
        for start, end in intervals:
            # Add buffer
            start = max(0, start - buffer)
            end = end + buffer
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        te_intervals[transcript_id] = merged

    te_regions = {}
    non_te_regions = {}

    for transcript_id, sequence in utr_seqs.items():
        seq_len = len(sequence)
        intervals = te_intervals.get(transcript_id, [])

        if not intervals:
            # No TE hits - entire sequence is non-TE
            non_te_regions[f"{transcript_id}_full"] = sequence
            continue

        # Extract TE regions
        for i, (start, end) in enumerate(intervals):
            end = min(end, seq_len)
            if end > start:
                region_seq = sequence[start:end]
                if len(region_seq) >= 30:  # Minimum length for structure
                    te_regions[f"{transcript_id}_te{i}_{start}_{end}"] = region_seq

        # Extract non-TE regions (gaps between TE hits)
        prev_end = 0
        for i, (start, end) in enumerate(intervals):
            if start > prev_end:
                gap_seq = sequence[prev_end:start]
                if len(gap_seq) >= 30:
                    non_te_regions[f"{transcript_id}_gap{i}_{prev_end}_{start}"] = gap_seq
            prev_end = max(prev_end, min(end, seq_len))

        # Final gap after last TE
        if prev_end < seq_len:
            gap_seq = sequence[prev_end:seq_len]
            if len(gap_seq) >= 30:
                non_te_regions[f"{transcript_id}_gap_final_{prev_end}_{seq_len}"] = gap_seq

    return te_regions, non_te_regions


def find_motif_in_structure(sequence: str, structure: str, motif: str) -> list:
    """
    Find motif positions and their structural context.

    Returns list of dicts with position and structure info.
    """
    results = []
    seq = sequence.upper().replace('U', 'T')
    motif = motif.upper()

    start = 0
    while True:
        pos = seq.find(motif, start)
        if pos == -1:
            break

        # Get structure at motif position
        motif_struct = structure[pos:pos+len(motif)] if pos+len(motif) <= len(structure) else ''
        n_paired = motif_struct.count('(') + motif_struct.count(')')
        n_unpaired = motif_struct.count('.')

        results.append({
            'position': pos,
            'motif_structure': motif_struct,
            'n_paired': n_paired,
            'n_unpaired': n_unpaired,
            'pct_paired': 100 * n_paired / len(motif_struct) if motif_struct else 0,
        })

        start = pos + 1

    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sample-size', type=int, default=100,
                        help='Number of sequences to sample for timing estimate')
    parser.add_argument('--full-run', action='store_true',
                        help='Run full analysis (skip timing estimate)')
    parser.add_argument('--max-seqs', type=int, default=5000,
                        help='Maximum sequences per category for full run')
    parser.add_argument('--buffer', type=int, default=20,
                        help='Buffer around TE hits (default: 20bp)')
    parser.add_argument('--motifs', type=str, default='AATAAA,ATTAAA,TGTATA,CAGCAG',
                        help='Motifs to analyze')
    parser.add_argument('--max-length', type=int, default=300,
                        help='Maximum sequence length for folding (default: 300bp)')
    parser.add_argument('--output', type=Path, default=None)
    args = parser.parse_args()

    project_root = get_project_root()
    results_dir = get_results_dir()
    struct_dir = results_dir / "structure_analysis"
    struct_dir.mkdir(exist_ok=True)

    output_path = args.output or struct_dir / "rna_structure_analysis.tsv"

    # File paths
    utr_fasta = project_root / "data" / "references" / "dmel_3utr.fasta"
    blast_path = results_dir / "genome_wide_all_3utrs.tsv"
    shuffled_dir = results_dir / "shuffled_full"

    motifs = [m.strip().upper() for m in args.motifs.split(',')]

    print("=" * 70, flush=True)
    print("RNA SECONDARY STRUCTURE ANALYSIS", flush=True)
    print("=" * 70, flush=True)
    print(f"Buffer around TE hits: {args.buffer}bp", flush=True)
    print(f"Max sequence length: {args.max_length}bp", flush=True)
    print(f"Motifs to analyze: {', '.join(motifs)}", flush=True)

    # Load UTR sequences
    print("\n1. Loading UTR sequences...")
    utr_seqs = {}
    for record in SeqIO.parse(utr_fasta, "fasta"):
        utr_seqs[record.id] = str(record.seq).upper()
    print(f"   Loaded {len(utr_seqs):,} UTR sequences")

    # Extract TE and non-TE regions
    print(f"\n2. Extracting TE and non-TE regions (buffer={args.buffer}bp)...")
    te_regions, non_te_regions = extract_te_regions(blast_path, utr_seqs, args.buffer)
    print(f"   TE regions: {len(te_regions):,}")
    print(f"   Non-TE regions: {len(non_te_regions):,}")

    # Load shuffled sequences (first replicate)
    print("\n3. Loading shuffled sequences...")
    shuf_fasta = shuffled_dir / "replicate_01.fasta"
    shuf_seqs = {}
    for record in SeqIO.parse(shuf_fasta, "fasta"):
        shuf_seqs[record.id] = str(record.seq).upper()
    print(f"   Loaded {len(shuf_seqs):,} shuffled sequences")

    # Extract TE regions from shuffled
    shuf_blast = shuffled_dir / "replicate_01_blast.tsv"
    shuf_te_regions, shuf_non_te_regions = extract_te_regions(shuf_blast, shuf_seqs, args.buffer)
    print(f"   Shuffled TE regions: {len(shuf_te_regions):,}")
    print(f"   Shuffled non-TE regions: {len(shuf_non_te_regions):,}")

    # Timing estimate
    if not args.full_run:
        print(f"\n4. Running timing estimate with {args.sample_size} sequences...")

        # Sample from each category
        sample_te = dict(list(te_regions.items())[:args.sample_size])
        sample_non_te = dict(list(non_te_regions.items())[:args.sample_size])

        start_time = time.time()
        results_te = run_rnafold(sample_te)
        results_non_te = run_rnafold(sample_non_te)
        elapsed = time.time() - start_time

        n_processed = len(results_te) + len(results_non_te)
        time_per_seq = elapsed / n_processed if n_processed > 0 else 1

        total_seqs = min(len(te_regions), args.max_seqs) + min(len(non_te_regions), args.max_seqs)
        total_seqs += min(len(shuf_te_regions), args.max_seqs) + min(len(shuf_non_te_regions), args.max_seqs)

        est_time = time_per_seq * total_seqs
        est_hours = est_time / 3600

        print(f"\n   Timing Results:")
        print(f"   - Processed {n_processed} sequences in {elapsed:.1f}s")
        print(f"   - Time per sequence: {time_per_seq:.3f}s")
        print(f"   - Total sequences for full run: {total_seqs:,}")
        print(f"   - Estimated full run time: {est_time:.0f}s ({est_hours:.2f} hours)")

        if est_hours > 2:
            print(f"\n   ⚠️  Estimated time > 2 hours. Consider reducing --max-seqs or running --full-run manually.")
            print(f"   Current --max-seqs: {args.max_seqs}")
            suggested = int(args.max_seqs * (2 * 3600) / est_time)
            print(f"   Suggested --max-seqs for 2hr limit: {suggested}")
            return 0
        else:
            print(f"\n   ✓ Estimated time < 2 hours. Proceeding with full analysis...")
            args.full_run = True

    # Full analysis
    if args.full_run:
        print(f"\n5. Running full structure analysis (max {args.max_seqs} per category)...")

        # Filter by max length and limit sample sizes
        def filter_by_length(regions, max_len, max_n):
            filtered = {k: v for k, v in regions.items() if len(v) <= max_len}
            return dict(list(filtered.items())[:max_n])

        te_sample = filter_by_length(te_regions, args.max_length, args.max_seqs)
        non_te_sample = filter_by_length(non_te_regions, args.max_length, args.max_seqs)
        shuf_te_sample = filter_by_length(shuf_te_regions, args.max_length, args.max_seqs)
        shuf_non_te_sample = filter_by_length(shuf_non_te_regions, args.max_length, args.max_seqs)

        print(f"\n   After length filter (≤{args.max_length}bp):", flush=True)
        print(f"   Real TE: {len(te_sample)}, Non-TE: {len(non_te_sample)}", flush=True)
        print(f"   Shuf TE: {len(shuf_te_sample)}, Non-TE: {len(shuf_non_te_sample)}", flush=True)

        print(f"\n   Folding real TE regions ({len(te_sample)})...")
        start = time.time()
        results_te = run_rnafold(te_sample)
        print(f"   Done in {time.time()-start:.1f}s")

        print(f"\n   Folding real non-TE regions ({len(non_te_sample)})...")
        start = time.time()
        results_non_te = run_rnafold(non_te_sample)
        print(f"   Done in {time.time()-start:.1f}s")

        print(f"\n   Folding shuffled TE regions ({len(shuf_te_sample)})...")
        start = time.time()
        results_shuf_te = run_rnafold(shuf_te_sample)
        print(f"   Done in {time.time()-start:.1f}s")

        print(f"\n   Folding shuffled non-TE regions ({len(shuf_non_te_sample)})...")
        start = time.time()
        results_shuf_non_te = run_rnafold(shuf_non_te_sample)
        print(f"   Done in {time.time()-start:.1f}s")

        # Compile results
        print("\n6. Analyzing results...")

        def summarize_structures(results: dict, label: str) -> dict:
            if not results:
                return {}
            mfes = [r['mfe_per_nt'] for r in results.values()]
            pct_paired = [r['pct_paired'] for r in results.values()]
            return {
                'category': label,
                'n_seqs': len(results),
                'mean_mfe_per_nt': np.mean(mfes),
                'std_mfe_per_nt': np.std(mfes),
                'median_mfe_per_nt': np.median(mfes),
                'mean_pct_paired': np.mean(pct_paired),
                'std_pct_paired': np.std(pct_paired),
            }

        summaries = [
            summarize_structures(results_te, 'real_te'),
            summarize_structures(results_non_te, 'real_non_te'),
            summarize_structures(results_shuf_te, 'shuffled_te'),
            summarize_structures(results_shuf_non_te, 'shuffled_non_te'),
        ]

        summary_df = pd.DataFrame([s for s in summaries if s])

        # Print comparison
        print("\n" + "=" * 70)
        print("RESULTS: Structural Properties by Region Type")
        print("=" * 70)

        print(f"\n{'Category':<20} {'N seqs':>10} {'MFE/nt':>12} {'± std':>10} {'% paired':>12}")
        print("-" * 70)
        for _, row in summary_df.iterrows():
            print(f"{row['category']:<20} {row['n_seqs']:>10,} "
                  f"{row['mean_mfe_per_nt']:>12.4f} {row['std_mfe_per_nt']:>10.4f} "
                  f"{row['mean_pct_paired']:>11.1f}%")

        # Statistical comparison
        print("\n" + "=" * 70)
        print("COMPARISON: TE vs Non-TE Regions")
        print("=" * 70)

        if results_te and results_non_te:
            te_mfe = [r['mfe_per_nt'] for r in results_te.values()]
            non_te_mfe = [r['mfe_per_nt'] for r in results_non_te.values()]

            diff_mfe = np.mean(te_mfe) - np.mean(non_te_mfe)
            print(f"\nMFE per nucleotide:")
            print(f"  TE regions:     {np.mean(te_mfe):.4f} ± {np.std(te_mfe):.4f}")
            print(f"  Non-TE regions: {np.mean(non_te_mfe):.4f} ± {np.std(non_te_mfe):.4f}")
            print(f"  Difference:     {diff_mfe:+.4f} (negative = more stable)")

            te_paired = [r['pct_paired'] for r in results_te.values()]
            non_te_paired = [r['pct_paired'] for r in results_non_te.values()]

            diff_paired = np.mean(te_paired) - np.mean(non_te_paired)
            print(f"\nPercent base-paired:")
            print(f"  TE regions:     {np.mean(te_paired):.1f}% ± {np.std(te_paired):.1f}%")
            print(f"  Non-TE regions: {np.mean(non_te_paired):.1f}% ± {np.std(non_te_paired):.1f}%")
            print(f"  Difference:     {diff_paired:+.1f}%")

        # Real vs Shuffled comparison
        print("\n" + "=" * 70)
        print("COMPARISON: Real vs Shuffled")
        print("=" * 70)

        if results_te and results_shuf_te:
            real_mfe = [r['mfe_per_nt'] for r in results_te.values()]
            shuf_mfe = [r['mfe_per_nt'] for r in results_shuf_te.values()]

            diff = np.mean(real_mfe) - np.mean(shuf_mfe)
            print(f"\nTE regions - MFE per nucleotide:")
            print(f"  Real:     {np.mean(real_mfe):.4f} ± {np.std(real_mfe):.4f}")
            print(f"  Shuffled: {np.mean(shuf_mfe):.4f} ± {np.std(shuf_mfe):.4f}")
            print(f"  Difference: {diff:+.4f}")

        # Motif structural context
        print("\n" + "=" * 70)
        print("MOTIF STRUCTURAL CONTEXT")
        print("=" * 70)

        motif_results = []

        for motif in motifs:
            te_motif_data = []
            non_te_motif_data = []

            # Analyze in TE regions
            for name, result in results_te.items():
                seq = te_sample.get(name, '')
                struct = result.get('structure', '')
                if seq and struct:
                    hits = find_motif_in_structure(seq, struct, motif)
                    te_motif_data.extend(hits)

            # Analyze in non-TE regions
            for name, result in results_non_te.items():
                seq = non_te_sample.get(name, '')
                struct = result.get('structure', '')
                if seq and struct:
                    hits = find_motif_in_structure(seq, struct, motif)
                    non_te_motif_data.extend(hits)

            if te_motif_data:
                te_pct_paired = np.mean([h['pct_paired'] for h in te_motif_data])
            else:
                te_pct_paired = np.nan

            if non_te_motif_data:
                non_te_pct_paired = np.mean([h['pct_paired'] for h in non_te_motif_data])
            else:
                non_te_pct_paired = np.nan

            motif_results.append({
                'motif': motif,
                'n_in_te': len(te_motif_data),
                'n_in_non_te': len(non_te_motif_data),
                'pct_paired_in_te': te_pct_paired,
                'pct_paired_in_non_te': non_te_pct_paired,
            })

            print(f"\n{motif}:")
            print(f"  In TE regions:     {len(te_motif_data):,} occurrences, {te_pct_paired:.1f}% base-paired")
            print(f"  In non-TE regions: {len(non_te_motif_data):,} occurrences, {non_te_pct_paired:.1f}% base-paired")

        # Save results
        summary_df.to_csv(output_path, sep='\t', index=False)

        motif_df = pd.DataFrame(motif_results)
        motif_output = struct_dir / "motif_structure_context.tsv"
        motif_df.to_csv(motif_output, sep='\t', index=False)

        print(f"\n\nResults saved to:")
        print(f"  {output_path}")
        print(f"  {motif_output}")

        # Interpretation
        print("\n" + "=" * 70)
        print("INTERPRETATION")
        print("=" * 70)
        print("""
MFE (Minimum Free Energy):
  - More negative = more stable structure
  - MFE per nucleotide normalizes for length

Percent base-paired:
  - Higher = more structured
  - Lower = more single-stranded / accessible

If TE-derived sequences have functional RNA structures:
  - Expect different MFE between TE and non-TE regions
  - Expect real to differ from shuffled (non-random structure)
  - Regulatory motifs may be in structured or unstructured regions
    depending on their function
""")

    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)

    return 0


if __name__ == '__main__':
    sys.exit(main())
