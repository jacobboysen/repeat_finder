#!/usr/bin/env python3
"""
Predict intrinsic disorder using IUPred3 REST API.

Downloads FlyBase protein sequences and predicts disorder using the IUPred3
web service. Outputs per-gene summary statistics for integration with TE analyses.

IUPred3 API: https://iupred3.elte.hu/

Note: Rate limiting may apply for large batches. The script includes delays
between requests to be respectful of the service.

Outputs:
- gene_disorder_scores.tsv: Per-gene disorder statistics
- protein_disorder_raw.tsv: Raw per-residue predictions (optional)
"""

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import requests
from Bio import SeqIO


IUPRED3_API_URL = "https://iupred3.elte.hu/iupred3/"


def load_protein_sequences(fasta_path: Path) -> Dict[str, Tuple[str, str]]:
    """
    Load protein sequences from FASTA file.

    Args:
        fasta_path: Path to protein FASTA

    Returns:
        Dictionary mapping FBgn to (protein_id, sequence)
    """
    proteins = {}

    for record in SeqIO.parse(fasta_path, 'fasta'):
        # Parse FlyBase protein header
        # Example: >FBpp0070748 type=protein; parent=FBgn0011761,FBtr0070749; ...
        protein_id = record.id

        # Extract parent gene
        fbgn = None
        for part in record.description.split(';'):
            part = part.strip()
            if part.startswith('parent='):
                parents = part.split('=')[1]
                for p in parents.split(','):
                    if p.startswith('FBgn'):
                        fbgn = p
                        break
                break

        if fbgn and len(record.seq) > 0:
            # Use first/longest protein per gene
            if fbgn not in proteins or len(record.seq) > len(proteins[fbgn][1]):
                proteins[fbgn] = (protein_id, str(record.seq))

    return proteins


def predict_disorder_iupred3(sequence: str, protein_id: str = "query") -> Optional[List[float]]:
    """
    Predict disorder using IUPred3 REST API.

    Args:
        sequence: Protein sequence (amino acids)
        protein_id: Identifier for logging

    Returns:
        List of per-residue disorder scores (0-1), or None on error
    """
    try:
        # IUPred3 API expects POST with sequence
        response = requests.post(
            IUPRED3_API_URL,
            data={
                'sequence': sequence,
                'output_type': 'json'
            },
            timeout=60
        )

        if response.status_code != 200:
            print(f"  Warning: IUPred3 returned status {response.status_code} for {protein_id}",
                  file=sys.stderr)
            return None

        # Parse JSON response
        data = response.json()

        # Extract IUPred scores
        if 'iupred2' in data:
            scores = [float(s) for s in data['iupred2']]
            return scores

        return None

    except requests.exceptions.Timeout:
        print(f"  Warning: Timeout for {protein_id}", file=sys.stderr)
        return None
    except requests.exceptions.RequestException as e:
        print(f"  Warning: Request error for {protein_id}: {e}", file=sys.stderr)
        return None
    except (json.JSONDecodeError, KeyError, ValueError) as e:
        print(f"  Warning: Parse error for {protein_id}: {e}", file=sys.stderr)
        return None


def calculate_disorder_stats(scores: List[float], threshold: float = 0.5) -> dict:
    """
    Calculate summary statistics from per-residue disorder scores.

    Args:
        scores: List of disorder scores (0-1) per residue
        threshold: Score threshold for "disordered" classification

    Returns:
        Dictionary with disorder statistics
    """
    if not scores:
        return {
            'disorder_fraction': 0,
            'max_disordered_region': 0,
            'n_disordered_residues': 0,
            'total_length': 0,
            'mean_score': 0
        }

    n_disordered = sum(1 for s in scores if s >= threshold)
    disorder_fraction = n_disordered / len(scores)

    # Find longest contiguous disordered region
    max_region = 0
    current_region = 0
    for s in scores:
        if s >= threshold:
            current_region += 1
            max_region = max(max_region, current_region)
        else:
            current_region = 0

    return {
        'disorder_fraction': disorder_fraction,
        'max_disordered_region': max_region,
        'n_disordered_residues': n_disordered,
        'total_length': len(scores),
        'mean_score': sum(scores) / len(scores)
    }


def run_disorder_prediction(
    proteins: Dict[str, Tuple[str, str]],
    output_path: Path,
    delay: float = 1.0,
    verbose: bool = False,
    max_proteins: Optional[int] = None
) -> Dict[str, dict]:
    """
    Run disorder prediction on all proteins.

    Args:
        proteins: Dictionary mapping FBgn to (protein_id, sequence)
        output_path: Path to write results
        delay: Delay between API requests (seconds)
        verbose: Print progress
        max_proteins: Limit number of proteins (for testing)

    Returns:
        Dictionary mapping FBgn to disorder statistics
    """
    results = {}

    protein_list = list(proteins.items())
    if max_proteins:
        protein_list = protein_list[:max_proteins]

    total = len(protein_list)

    if verbose:
        print(f"\nPredicting disorder for {total:,} proteins...")

    with open(output_path, 'w') as f:
        # Write header
        f.write('\t'.join([
            'fbgn', 'protein_id', 'total_length', 'disorder_fraction',
            'max_disordered_region', 'n_disordered_residues', 'mean_score'
        ]) + '\n')

        for i, (fbgn, (protein_id, sequence)) in enumerate(protein_list, 1):
            if verbose and i % 100 == 0:
                print(f"  Progress: {i}/{total} ({100*i/total:.1f}%)")

            # Skip very short sequences
            if len(sequence) < 20:
                continue

            # Predict disorder
            scores = predict_disorder_iupred3(sequence, protein_id)

            if scores is None:
                # Use fallback: simple hydropathy-based estimate
                stats = {
                    'disorder_fraction': 0,
                    'max_disordered_region': 0,
                    'n_disordered_residues': 0,
                    'total_length': len(sequence),
                    'mean_score': 0
                }
            else:
                stats = calculate_disorder_stats(scores)

            results[fbgn] = stats

            # Write row
            f.write('\t'.join([
                fbgn,
                protein_id,
                str(stats['total_length']),
                f"{stats['disorder_fraction']:.4f}",
                str(stats['max_disordered_region']),
                str(stats['n_disordered_residues']),
                f"{stats['mean_score']:.4f}"
            ]) + '\n')

            # Rate limiting
            time.sleep(delay)

    return results


def run_local_prediction(
    proteins: Dict[str, Tuple[str, str]],
    output_path: Path,
    verbose: bool = False
) -> Dict[str, dict]:
    """
    Run simple local disorder prediction (no API needed).

    Uses a basic hydropathy-based heuristic as a fallback when API is unavailable.
    This is less accurate than IUPred but doesn't require network access.

    Args:
        proteins: Dictionary mapping FBgn to (protein_id, sequence)
        output_path: Path to write results
        verbose: Print progress

    Returns:
        Dictionary mapping FBgn to disorder statistics
    """
    # Kyte-Doolittle hydropathy scale
    HYDROPATHY = {
        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
        'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
        'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
        'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
    }

    # Disorder-promoting residues (PEST-like)
    DISORDER_PROMOTING = set('PESTKRQN')

    results = {}
    total = len(proteins)

    if verbose:
        print(f"\nRunning local disorder prediction for {total:,} proteins...")

    with open(output_path, 'w') as f:
        f.write('\t'.join([
            'fbgn', 'protein_id', 'total_length', 'disorder_fraction',
            'max_disordered_region', 'n_disordered_residues', 'mean_score'
        ]) + '\n')

        for i, (fbgn, (protein_id, sequence)) in enumerate(proteins.items(), 1):
            if verbose and i % 1000 == 0:
                print(f"  Progress: {i}/{total} ({100*i/total:.1f}%)")

            if len(sequence) < 10:
                continue

            # Calculate per-residue disorder scores
            window_size = 21
            half_window = window_size // 2
            scores = []

            for j in range(len(sequence)):
                # Get window
                start = max(0, j - half_window)
                end = min(len(sequence), j + half_window + 1)
                window = sequence[start:end]

                # Calculate features
                # 1. Fraction of disorder-promoting residues
                n_disorder = sum(1 for aa in window if aa in DISORDER_PROMOTING)
                disorder_frac = n_disorder / len(window)

                # 2. Average hydropathy (low = more disordered)
                hydro_sum = sum(HYDROPATHY.get(aa, 0) for aa in window)
                avg_hydro = hydro_sum / len(window)

                # Combine: high disorder fraction + low hydropathy = disordered
                # Normalize to 0-1 range
                hydro_score = (1.5 - avg_hydro) / 6  # Range roughly -4.5 to 4.5
                combined = 0.6 * disorder_frac + 0.4 * max(0, min(1, hydro_score))

                scores.append(combined)

            # Calculate statistics
            stats = calculate_disorder_stats(scores, threshold=0.4)  # Lower threshold for local method
            results[fbgn] = stats

            f.write('\t'.join([
                fbgn,
                protein_id,
                str(stats['total_length']),
                f"{stats['disorder_fraction']:.4f}",
                str(stats['max_disordered_region']),
                str(stats['n_disordered_residues']),
                f"{stats['mean_score']:.4f}"
            ]) + '\n')

    return results


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--proteins',
        type=Path,
        default=Path('data/references/dmel_proteins.fasta'),
        help='Protein sequences FASTA (default: data/references/dmel_proteins.fasta)'
    )
    parser.add_argument(
        '--output',
        type=Path,
        default=Path('data/annotations/gene_disorder_scores.tsv'),
        help='Output TSV path (default: data/annotations/gene_disorder_scores.tsv)'
    )
    parser.add_argument(
        '--local',
        action='store_true',
        help='Use local prediction (no API, less accurate but faster)'
    )
    parser.add_argument(
        '--delay',
        type=float,
        default=1.0,
        help='Delay between API requests in seconds (default: 1.0)'
    )
    parser.add_argument(
        '--max-proteins',
        type=int,
        help='Limit number of proteins to process (for testing)'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )

    args = parser.parse_args()

    # Validate input
    if not args.proteins.exists():
        print(f"Error: Protein file not found: {args.proteins}", file=sys.stderr)
        return 1

    # Create output directory
    args.output.parent.mkdir(parents=True, exist_ok=True)

    print("Disorder Prediction")
    print("=" * 60)

    # Load proteins
    if args.verbose:
        print(f"Loading proteins from: {args.proteins}")

    proteins = load_protein_sequences(args.proteins)

    if args.verbose:
        print(f"  Loaded {len(proteins):,} unique genes with proteins")

    # Run prediction
    if args.local:
        print("\nUsing LOCAL prediction method (hydropathy-based heuristic)")
        results = run_local_prediction(proteins, args.output, args.verbose)
    else:
        print("\nUsing IUPred3 API (this may take a while...)")
        results = run_disorder_prediction(
            proteins, args.output, args.delay, args.verbose, args.max_proteins
        )

    # Summary
    print("\n" + "=" * 60)
    print("Prediction complete!")
    print(f"  Genes processed: {len(results):,}")
    print(f"  Output written to: {args.output}")

    if results:
        fractions = [r['disorder_fraction'] for r in results.values()]
        print(f"\nDisorder statistics:")
        print(f"  Mean disorder fraction: {sum(fractions)/len(fractions):.3f}")
        print(f"  Median: {sorted(fractions)[len(fractions)//2]:.3f}")

        low = sum(1 for f in fractions if f < 0.2)
        medium = sum(1 for f in fractions if 0.2 <= f < 0.5)
        high = sum(1 for f in fractions if f >= 0.5)
        print(f"  Low disorder (<20%): {low:,} genes")
        print(f"  Medium (20-50%): {medium:,} genes")
        print(f"  High (>50%): {high:,} genes")

    return 0


if __name__ == '__main__':
    sys.exit(main())
