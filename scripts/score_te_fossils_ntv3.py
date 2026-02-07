#!/usr/bin/env python3
"""
Score TE fossils in regulatory regions using Nucleotide Transformer v3.

This script takes the LM-ready dataset (from prepare_te_fossils_for_lm.py) and
runs three analyses through NTv3's post-trained Drosophila functional track heads:

1. **Baseline scoring**: Predict regulatory activity (ATAC, ChIP, RNA-seq tracks)
   across each region. Compare signal at TE-derived vs non-TE positions.

2. **Ablation**: Shuffle or mask the TE-derived subsequence, re-run NTv3, and
   measure the delta in predicted regulatory activity. A large delta indicates
   the TE fossil contributes to regulatory function.

3. **In-silico saturation mutagenesis**: For top candidates, compute every
   possible SNV within the TE-derived region and score through NTv3. Generates
   per-position functional constraint maps.

Requires:
  - GPU with >=40GB VRAM (A100 recommended)
  - transformers >= 4.55.0
  - torch >= 2.0
  - pyfaidx

Usage:
    # Score all TE-bearing regions (baseline + ablation)
    python scripts/score_te_fossils_ntv3.py -v

    # Also run saturation mutagenesis on top N candidates
    python scripts/score_te_fossils_ntv3.py -v --satmut --top-n 100

    # Quick test on a small subset
    python scripts/score_te_fossils_ntv3.py -v --max-regions 50 --satmut --top-n 10
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_project_root, get_results_dir


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
# NTv3 requires input length divisible by 2^7 = 128
DIVISOR = 128
# Post-trained model crops output to center 37.5%
KEEP_FRACTION = 0.375
# Nucleotide encoding for saturation mutagenesis
NUC_TOKENS = {'A': 6, 'T': 7, 'C': 8, 'G': 9}  # will be set from tokenizer
NUCS = ['A', 'T', 'C', 'G']


# ---------------------------------------------------------------------------
# Genome sequence provider
# ---------------------------------------------------------------------------
class GenomeProvider:
    """Lazy-loading genome FASTA provider using pyfaidx."""

    def __init__(self, fasta_path):
        from pyfaidx import Fasta
        fai_path = Path(str(fasta_path) + '.fai')
        self.fasta = Fasta(str(fasta_path), rebuild=not fai_path.exists())
        # Map chromosome name variants
        self.chrom_map = {}
        for name in self.fasta.keys():
            self.chrom_map[name] = name
            # Also map with/without 'chr' prefix
            if name.startswith('chr'):
                self.chrom_map[name[3:]] = name
            else:
                self.chrom_map['chr' + name] = name

    def fetch(self, chrom, start, end):
        """Fetch sequence (0-based half-open coordinates)."""
        key = self.chrom_map.get(chrom)
        if key is None:
            return None
        chrom_len = len(self.fasta[key])
        # Clamp to chromosome bounds
        s = max(0, start)
        e = min(chrom_len, end)
        seq = str(self.fasta[key][s:e]).upper()
        # Pad if we went past chromosome edges
        left_pad = s - start
        right_pad = end - e
        if left_pad > 0:
            seq = 'N' * left_pad + seq
        if right_pad > 0:
            seq = seq + 'N' * right_pad
        return seq

    def chrom_length(self, chrom):
        key = self.chrom_map.get(chrom)
        if key is None:
            return 0
        return len(self.fasta[key])


# ---------------------------------------------------------------------------
# NTv3 model wrapper
# ---------------------------------------------------------------------------
class NTv3Scorer:
    """Wrapper for NTv3 post-trained model inference."""

    def __init__(self, model_name="InstaDeepAI/NTv3_100M_post", device=None,
                 use_bf16=True):
        from transformers import AutoTokenizer, AutoModel

        if device is None:
            if torch.cuda.is_available():
                device = 'cuda'
            elif hasattr(torch.backends, 'mps') and torch.backends.mps.is_available():
                device = 'mps'
            else:
                device = 'cpu'
        self.device = torch.device(device)

        print(f"Loading NTv3 model: {model_name}")
        print(f"Device: {self.device}")

        self.tokenizer = AutoTokenizer.from_pretrained(
            model_name, trust_remote_code=True
        )

        dtype_kwargs = {}
        if use_bf16 and self.device.type == 'cuda':
            dtype_kwargs = {
                'stem_compute_dtype': 'bfloat16',
                'down_convolution_compute_dtype': 'bfloat16',
                'transformer_qkvo_compute_dtype': 'bfloat16',
                'transformer_ffn_compute_dtype': 'bfloat16',
                'up_convolution_compute_dtype': 'bfloat16',
                'modulation_compute_dtype': 'bfloat16',
            }

        self.model = AutoModel.from_pretrained(
            model_name, trust_remote_code=True, **dtype_kwargs
        ).to(self.device).eval()

        # Discover species key for Drosophila
        species_map = self.model.config.species_to_token_id
        self.species_key = None
        for candidate in ['fruit_fly', 'drosophila', 'drosophila_melanogaster',
                          'dm6', 'dm', 'Drosophila_melanogaster']:
            if candidate in species_map:
                self.species_key = candidate
                break
        if self.species_key is None:
            # Try case-insensitive search
            for key in species_map:
                if 'droso' in key.lower() or 'fly' in key.lower() or 'dm' in key.lower():
                    self.species_key = key
                    break
        if self.species_key is None:
            print(f"WARNING: Could not find Drosophila in species map.")
            print(f"  Available species: {list(species_map.keys())}")
            print(f"  Falling back to first available species.")
            self.species_key = list(species_map.keys())[0]

        print(f"Using species key: '{self.species_key}'")

        # Get track names
        self.bigwig_names = list(
            self.model.config.bigwigs_per_species.get(self.species_key, [])
        )
        self.bed_element_names = list(self.model.config.bed_elements_names)
        self.keep_fraction = getattr(
            self.model.config, 'keep_target_center_fraction', KEEP_FRACTION
        )

        print(f"BigWig tracks for {self.species_key}: {len(self.bigwig_names)}")
        print(f"BED element types: {len(self.bed_element_names)}")
        if self.bigwig_names:
            print(f"  Example tracks: {self.bigwig_names[:5]}")
        if self.bed_element_names:
            print(f"  Example elements: {self.bed_element_names[:5]}")

        # Build nucleotide token ID map from tokenizer
        self.nuc_token_ids = {}
        for nuc in NUCS:
            ids = self.tokenizer.encode(nuc, add_special_tokens=False)
            if ids:
                self.nuc_token_ids[nuc] = ids[0]

    def tokenize(self, sequence):
        """Tokenize a single sequence, padded to multiple of 128."""
        batch = self.tokenizer(
            [sequence],
            add_special_tokens=False,
            padding=True,
            pad_to_multiple_of=DIVISOR,
            return_tensors='pt',
        )
        return batch['input_ids'].to(self.device)

    @torch.no_grad()
    def predict(self, input_ids):
        """Run forward pass, return bigwig and bed predictions."""
        species_ids = self.model.encode_species([self.species_key])
        species_ids = species_ids.to(self.device)

        out = self.model(
            input_ids=input_ids,
            species_ids=species_ids,
        )

        result = {}
        if out.bigwig_tracks_logits is not None:
            result['bigwig'] = out.bigwig_tracks_logits[0].cpu().float().numpy()
        if out.bed_tracks_logits is not None:
            result['bed'] = out.bed_tracks_logits[0].cpu().float().numpy()
        if out.logits is not None:
            result['logits'] = out.logits[0].cpu().float().numpy()

        return result

    def score_sequence(self, sequence):
        """Score a single sequence, return predictions dict."""
        input_ids = self.tokenize(sequence)
        return self.predict(input_ids)

    def score_mutant_batch(self, sequences):
        """Score a batch of sequences (for saturation mutagenesis)."""
        batch = self.tokenizer(
            sequences,
            add_special_tokens=False,
            padding=True,
            pad_to_multiple_of=DIVISOR,
            return_tensors='pt',
        )
        input_ids = batch['input_ids'].to(self.device)
        species_ids = self.model.encode_species(
            [self.species_key] * len(sequences)
        )
        species_ids = species_ids.to(self.device)

        with torch.no_grad():
            out = self.model(
                input_ids=input_ids,
                species_ids=species_ids,
            )

        result = {}
        if out.bigwig_tracks_logits is not None:
            result['bigwig'] = out.bigwig_tracks_logits.cpu().float().numpy()
        del out, input_ids, species_ids
        torch.cuda.empty_cache()
        return result


# ---------------------------------------------------------------------------
# Window construction
# ---------------------------------------------------------------------------
def build_genomic_window(genome, chrom, region_start, region_end, region_strand,
                         target_window=32768):
    """
    Build a genomic window centered on the region, sized to target_window.

    Returns:
        sequence: str of length target_window (padded with N if at chrom edge)
        region_offset: int, 0-based offset of region_start within the window
        pred_start: int, 0-based start of the prediction window within the sequence
        pred_end: int, 0-based end of the prediction window
    """
    region_len = region_end - region_start
    region_center = (region_start + region_end) // 2

    # Center the window on the region
    window_start = region_center - target_window // 2
    window_end = window_start + target_window

    seq = genome.fetch(chrom, window_start, window_end)
    if seq is None:
        return None, None, None, None

    region_offset = region_start - window_start

    # Prediction window = center KEEP_FRACTION of the input
    margin = int(target_window * (1 - KEEP_FRACTION) / 2)
    pred_start = margin
    pred_end = target_window - margin

    return seq, region_offset, pred_start, pred_end


def shuffle_region(sequence, start, end, rng=None):
    """Dinucleotide-preserving shuffle of a subsequence within a larger sequence."""
    if rng is None:
        rng = np.random.default_rng()

    prefix = sequence[:start]
    region = list(sequence[start:end])
    suffix = sequence[end:]

    # Simple shuffle (preserves composition but not dinucleotide)
    rng.shuffle(region)

    return prefix + ''.join(region) + suffix


def mask_region(sequence, start, end):
    """Replace a subsequence with N's."""
    return sequence[:start] + 'N' * (end - start) + sequence[end:]


# ---------------------------------------------------------------------------
# Analysis functions
# ---------------------------------------------------------------------------
def compute_track_summary(bigwig, pred_start, pred_end, region_offset,
                          region_len, mask, window_start_genomic,
                          bigwig_names):
    """
    Compute summary statistics comparing TE vs non-TE positions.

    Args:
        bigwig: (pred_len, n_tracks) array of predicted track values
        pred_start/pred_end: prediction window offsets within full sequence
        region_offset: offset of region within full sequence
        region_len: length of the regulatory region
        mask: uint8 array of TE mask for the region
        window_start_genomic: genomic start of the full window
        bigwig_names: list of track names

    Returns:
        dict with summary statistics
    """
    pred_len = bigwig.shape[0]
    n_tracks = bigwig.shape[1]

    # Map region positions to prediction window positions
    # Region covers [region_offset, region_offset + region_len) in full seq
    # Prediction covers [pred_start, pred_end) in full seq
    reg_in_pred_start = max(0, region_offset - pred_start)
    reg_in_pred_end = min(pred_len, region_offset + region_len - pred_start)

    if reg_in_pred_end <= reg_in_pred_start:
        return None  # Region is outside prediction window

    # Extract region predictions
    region_tracks = bigwig[reg_in_pred_start:reg_in_pred_end]  # (reg_len_in_pred, n_tracks)

    # Map mask to prediction window coordinates
    mask_start_in_region = max(0, pred_start - region_offset)
    mask_end_in_region = min(region_len, pred_end - region_offset)
    region_mask = mask[mask_start_in_region:mask_end_in_region]

    # Ensure mask and region_tracks have same length
    actual_len = min(len(region_mask), region_tracks.shape[0])
    region_mask = region_mask[:actual_len]
    region_tracks = region_tracks[:actual_len]

    if actual_len == 0:
        return None

    # Compute per-track mean signal for TE vs non-TE positions
    te_positions = region_mask >= 3  # strict tier
    non_te_positions = region_mask == 0
    all_positions = np.ones(actual_len, dtype=bool)

    result = {
        'n_positions_scored': int(actual_len),
        'n_te_positions': int(te_positions.sum()),
        'n_non_te_positions': int(non_te_positions.sum()),
    }

    # Mean signal across all tracks
    result['mean_signal_all'] = float(region_tracks.mean())

    if te_positions.any():
        te_signal = region_tracks[te_positions].mean(axis=0)  # (n_tracks,)
        result['mean_signal_te'] = float(te_signal.mean())
        result['max_signal_te'] = float(te_signal.max())
        # Per-track TE signal
        result['te_track_means'] = te_signal.tolist()
    else:
        result['mean_signal_te'] = 0.0
        result['max_signal_te'] = 0.0

    if non_te_positions.any():
        non_te_signal = region_tracks[non_te_positions].mean(axis=0)
        result['mean_signal_non_te'] = float(non_te_signal.mean())
    else:
        result['mean_signal_non_te'] = 0.0

    # TE vs non-TE enrichment ratio
    if result['mean_signal_non_te'] > 0 and te_positions.any():
        result['te_enrichment'] = result['mean_signal_te'] / result['mean_signal_non_te']
    else:
        result['te_enrichment'] = float('nan')

    return result


def compute_ablation_delta(ref_bigwig, abl_bigwig, pred_start, pred_end,
                           region_offset, region_len, mask):
    """
    Compute the change in predicted signal when TE region is ablated.

    Returns dict with delta statistics.
    """
    pred_len = ref_bigwig.shape[0]

    reg_in_pred_start = max(0, region_offset - pred_start)
    reg_in_pred_end = min(pred_len, region_offset + region_len - pred_start)

    if reg_in_pred_end <= reg_in_pred_start:
        return None

    ref_region = ref_bigwig[reg_in_pred_start:reg_in_pred_end]
    abl_region = abl_bigwig[reg_in_pred_start:reg_in_pred_end]

    # Map mask
    mask_start = max(0, pred_start - region_offset)
    mask_end = min(region_len, pred_end - region_offset)
    region_mask = mask[mask_start:mask_end]

    actual_len = min(len(region_mask), ref_region.shape[0])
    region_mask = region_mask[:actual_len]
    ref_region = ref_region[:actual_len]
    abl_region = abl_region[:actual_len]

    delta = ref_region - abl_region  # positive = signal lost when TE removed

    te_pos = region_mask >= 3
    result = {
        'mean_delta_all': float(delta.mean()),
        'mean_abs_delta_all': float(np.abs(delta).mean()),
        'max_abs_delta': float(np.abs(delta).max()),
    }

    if te_pos.any():
        result['mean_delta_te'] = float(delta[te_pos].mean())
        result['mean_abs_delta_te'] = float(np.abs(delta[te_pos]).mean())

    # Track-wise delta (mean across positions)
    result['track_deltas'] = delta.mean(axis=0).tolist()

    # Fraction of tracks where signal decreased (TE was contributing)
    track_mean_delta = delta.mean(axis=0)
    result['frac_tracks_decreased'] = float(
        (track_mean_delta > 0).sum() / max(1, len(track_mean_delta))
    )

    return result


def run_saturation_mutagenesis(scorer, sequence, region_offset, region_len,
                               mask, pred_start, pred_end, batch_size=32):
    """
    In-silico saturation mutagenesis over TE-derived positions.

    For each TE position in the region, compute the effect of mutating to
    each of the 3 alternative nucleotides. Returns a per-position effect matrix.

    Returns:
        effect_matrix: (n_te_positions, 4) array, mean |delta| for each mutation
        positions: list of 0-based positions within the region
    """
    pred_len = pred_start  # just for offset calculation base
    te_positions = np.where(mask >= 3)[0]  # strict TE positions within region

    if len(te_positions) == 0:
        return None, []

    # Only mutagenize positions that fall within the prediction window
    valid = []
    for pos in te_positions:
        abs_pos = region_offset + pos
        if pred_start <= abs_pos < pred_end:
            valid.append(pos)

    if not valid:
        return None, []

    # Get reference prediction
    ref_preds = scorer.score_sequence(sequence)
    ref_bigwig = ref_preds.get('bigwig')
    if ref_bigwig is None:
        return None, []

    # Build mutant sequences in batches
    effect_matrix = np.zeros((len(valid), 4), dtype=np.float32)
    seq_list = list(sequence)

    mutant_specs = []  # (pos_idx, nuc_idx, mutant_sequence)
    for pi, pos in enumerate(valid):
        abs_pos = region_offset + pos
        ref_nuc = sequence[abs_pos]
        for ni, nuc in enumerate(NUCS):
            if nuc == ref_nuc:
                effect_matrix[pi, ni] = 0.0
            else:
                mut_seq = sequence[:abs_pos] + nuc + sequence[abs_pos + 1:]
                mutant_specs.append((pi, ni, mut_seq))

    # Score mutants in batches
    for batch_start in range(0, len(mutant_specs), batch_size):
        batch = mutant_specs[batch_start:batch_start + batch_size]
        seqs = [spec[2] for spec in batch]

        mut_preds = scorer.score_mutant_batch(seqs)
        mut_bigwig = mut_preds.get('bigwig')
        if mut_bigwig is None:
            continue

        for i, (pi, ni, _) in enumerate(batch):
            delta = np.abs(ref_bigwig - mut_bigwig[i]).mean()
            effect_matrix[pi, ni] = delta

    return effect_matrix, valid


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '--dataset-dir', type=Path,
        default=get_results_dir() / 'te_fossil_lm_dataset',
        help='Input dataset directory (from prepare_te_fossils_for_lm.py)',
    )
    parser.add_argument(
        '--genome', type=Path,
        default=get_project_root() / 'data' / 'references' / 'dmel_genome.fasta',
        help='Reference genome FASTA (indexed with samtools faidx)',
    )
    parser.add_argument(
        '--output-dir', type=Path,
        default=get_results_dir() / 'ntv3_te_scoring',
        help='Output directory for results',
    )
    parser.add_argument(
        '--model', type=str,
        default='InstaDeepAI/NTv3_100M_post',
        help='HuggingFace model ID',
    )
    parser.add_argument(
        '--window-size', type=int, default=32768,
        help='Genomic context window size (must be multiple of 128)',
    )
    parser.add_argument(
        '--device', type=str, default=None,
        help='Device (cuda, mps, cpu). Auto-detected if not set.',
    )
    parser.add_argument(
        '--satmut', action='store_true',
        help='Run saturation mutagenesis on top candidates',
    )
    parser.add_argument(
        '--top-n', type=int, default=100,
        help='Number of top candidates for saturation mutagenesis',
    )
    parser.add_argument(
        '--satmut-batch-size', type=int, default=4,
        help='Batch size for saturation mutagenesis forward passes',
    )
    parser.add_argument(
        '--max-regions', type=int, default=None,
        help='Limit number of regions to process (for testing)',
    )
    parser.add_argument(
        '--ablation-reps', type=int, default=5,
        help='Number of shuffle replicates for ablation',
    )
    parser.add_argument(
        '--min-te-bp', type=int, default=20,
        help='Minimum TE bp (strict) to include a region',
    )
    parser.add_argument(
        '-v', '--verbose', action='store_true',
    )

    args = parser.parse_args()

    assert args.window_size % DIVISOR == 0, \
        f"Window size must be divisible by {DIVISOR}"

    print("NTv3 TE Fossil Scoring Pipeline")
    print("=" * 60)

    # ── Load dataset ─────────────────────────────────────────────
    print("\nLoading dataset...")
    regions_df = pd.read_csv(args.dataset_dir / 'regions.tsv', sep='\t')
    masks_npz = np.load(str(args.dataset_dir / 'masks.npz'), allow_pickle=False)

    # Filter to TE-bearing regions above threshold
    te_regions = regions_df[regions_df['te_bp_strict'] >= args.min_te_bp].copy()
    te_regions = te_regions.sort_values('te_bp_strict', ascending=False)

    if args.max_regions:
        te_regions = te_regions.head(args.max_regions)

    print(f"  Total regions: {len(regions_df):,}")
    print(f"  TE-bearing (>={args.min_te_bp}bp strict): {len(te_regions):,}")

    # ── Load genome ──────────────────────────────────────────────
    print(f"\nLoading genome: {args.genome}")
    genome = GenomeProvider(args.genome)

    # ── Load model ───────────────────────────────────────────────
    scorer = NTv3Scorer(
        model_name=args.model,
        device=args.device,
        use_bf16=True,
    )

    # ── Prepare output ───────────────────────────────────────────
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # ── Check for checkpoint (resume support) ────────────────────
    checkpoint_path = args.output_dir / '_checkpoint.json'
    completed_rids = set()
    if checkpoint_path.exists():
        with open(checkpoint_path) as f:
            ckpt = json.load(f)
        completed_rids = set(ckpt.get('completed_region_ids', []))
        print(f"\nResuming from checkpoint: {len(completed_rids):,} regions already done")

    # ── Phase 1: Baseline scoring + ablation ─────────────────────
    print(f"\n{'='*60}")
    print(f"Phase 1: Baseline scoring + ablation ({len(te_regions):,} regions)")
    print(f"{'='*60}")

    scoring_results = []
    ablation_results = []
    rng = np.random.default_rng(42)

    # Open TSV files for streaming writes (append if resuming)
    scores_streaming_path = args.output_dir / '_scores_partial.tsv'
    ablation_streaming_path = args.output_dir / '_ablation_partial.tsv'
    score_columns = [
        'region_id', 'region_type', 'chrom', 'start', 'end', 'fbgn', 'symbol',
        'te_bp_strict', 'te_fraction_strict', 'n_te_hits_strict',
        'n_positions_scored', 'n_te_positions', 'n_non_te_positions',
        'mean_signal_all', 'mean_signal_te', 'max_signal_te',
        'mean_signal_non_te', 'te_enrichment',
    ]
    abl_columns = [
        'region_id', 'ablation_type', 'mean_delta_all', 'mean_abs_delta_all',
        'max_abs_delta', 'frac_tracks_decreased', 'mean_delta_te',
        'mean_abs_delta_te', 'n_reps',
    ]

    write_mode = 'a' if completed_rids else 'w'
    if write_mode == 'w':
        with open(scores_streaming_path, 'w') as f:
            f.write('\t'.join(score_columns) + '\n')
        with open(ablation_streaming_path, 'w') as f:
            f.write('\t'.join(abl_columns) + '\n')

    t0 = time.time()
    for idx, (_, row) in enumerate(te_regions.iterrows()):
        rid = row['region_id']

        # Skip if already completed (resume support)
        if rid in completed_rids:
            continue

        chrom = str(row['chrom'])
        region_start = int(row['start'])
        region_end = int(row['end'])
        region_strand = str(row['strand'])
        region_len = int(row['length'])

        if args.verbose and idx % 100 == 0:
            elapsed = time.time() - t0
            rate = idx / elapsed if elapsed > 0 else 0
            print(f"  [{idx:,}/{len(te_regions):,}] {rid} "
                  f"({rate:.1f} regions/sec)")

        # Get mask
        if rid in masks_npz.files:
            mask = masks_npz[rid]
        else:
            continue

        # Build genomic window
        seq, region_offset, pred_start, pred_end = build_genomic_window(
            genome, chrom, region_start, region_end, region_strand,
            target_window=args.window_size,
        )
        if seq is None:
            continue

        # Ensure region is within prediction window
        # (regions longer than the prediction window are partially scored)
        region_in_pred = (
            region_offset < pred_end and
            region_offset + region_len > pred_start
        )
        if not region_in_pred:
            if args.verbose:
                print(f"    SKIP {rid}: region outside prediction window")
            continue

        # --- Baseline scoring ---
        try:
            ref_preds = scorer.score_sequence(seq)
        except Exception as e:
            if args.verbose:
                print(f"    ERROR scoring {rid}: {e}")
            continue

        ref_bigwig = ref_preds.get('bigwig')
        if ref_bigwig is None:
            continue

        summary = compute_track_summary(
            ref_bigwig, pred_start, pred_end, region_offset, region_len,
            mask, region_start - region_offset, scorer.bigwig_names,
        )
        if summary is None:
            continue

        summary['region_id'] = rid
        summary['region_type'] = row['region_type']
        summary['chrom'] = chrom
        summary['start'] = region_start
        summary['end'] = region_end
        summary['fbgn'] = row['fbgn']
        summary['symbol'] = row['symbol']
        summary['te_bp_strict'] = int(row['te_bp_strict'])
        summary['te_fraction_strict'] = float(row['te_fraction_strict'])
        summary['n_te_hits_strict'] = int(row['n_te_hits_strict'])
        scoring_results.append(summary)

        # --- Ablation (shuffle + mask) ---
        # Find bounding box of all strict TE positions in the region
        te_indices = np.where(mask >= 3)[0]
        if len(te_indices) == 0:
            continue
        te_start_in_seq = region_offset + int(te_indices[0])
        te_end_in_seq = region_offset + int(te_indices[-1]) + 1

        # 1. N-mask ablation
        masked_seq = mask_region(seq, te_start_in_seq, te_end_in_seq)
        try:
            mask_preds = scorer.score_sequence(masked_seq)
        except Exception:
            mask_preds = None

        if mask_preds and mask_preds.get('bigwig') is not None:
            mask_delta = compute_ablation_delta(
                ref_bigwig, mask_preds['bigwig'],
                pred_start, pred_end, region_offset, region_len, mask,
            )
            if mask_delta:
                mask_delta['region_id'] = rid
                mask_delta['ablation_type'] = 'mask'
                ablation_results.append(mask_delta)

        # 2. Shuffle ablation (average over reps)
        shuffle_deltas = []
        for rep in range(args.ablation_reps):
            shuffled_seq = shuffle_region(
                seq, te_start_in_seq, te_end_in_seq, rng=rng
            )
            try:
                shuf_preds = scorer.score_sequence(shuffled_seq)
            except Exception:
                continue
            if shuf_preds.get('bigwig') is not None:
                sd = compute_ablation_delta(
                    ref_bigwig, shuf_preds['bigwig'],
                    pred_start, pred_end, region_offset, region_len, mask,
                )
                if sd:
                    shuffle_deltas.append(sd)

        if shuffle_deltas:
            avg_delta = {
                'region_id': rid,
                'ablation_type': 'shuffle',
                'mean_delta_all': float(np.mean(
                    [d['mean_delta_all'] for d in shuffle_deltas]
                )),
                'mean_abs_delta_all': float(np.mean(
                    [d['mean_abs_delta_all'] for d in shuffle_deltas]
                )),
                'max_abs_delta': float(np.mean(
                    [d['max_abs_delta'] for d in shuffle_deltas]
                )),
                'frac_tracks_decreased': float(np.mean(
                    [d['frac_tracks_decreased'] for d in shuffle_deltas]
                )),
                'n_reps': len(shuffle_deltas),
            }
            if 'mean_delta_te' in shuffle_deltas[0]:
                avg_delta['mean_delta_te'] = float(np.mean(
                    [d['mean_delta_te'] for d in shuffle_deltas]
                ))
                avg_delta['mean_abs_delta_te'] = float(np.mean(
                    [d['mean_abs_delta_te'] for d in shuffle_deltas]
                ))
            ablation_results.append(avg_delta)

        # Stream results to disk and update checkpoint
        completed_rids.add(rid)

        # Append score row to streaming TSV
        if scoring_results and scoring_results[-1]['region_id'] == rid:
            s = scoring_results[-1]
            with open(scores_streaming_path, 'a') as f:
                vals = [str(s.get(c, '')) for c in score_columns]
                f.write('\t'.join(vals) + '\n')

        # Append ablation rows
        new_abl = [a for a in ablation_results if a.get('region_id') == rid]
        if new_abl:
            with open(ablation_streaming_path, 'a') as f:
                for a in new_abl:
                    vals = [str(a.get(c, '')) for c in abl_columns]
                    f.write('\t'.join(vals) + '\n')

        # Save checkpoint every 50 regions
        if len(completed_rids) % 50 == 0:
            with open(checkpoint_path, 'w') as f:
                json.dump({'completed_region_ids': list(completed_rids)}, f)

    # Final checkpoint
    with open(checkpoint_path, 'w') as f:
        json.dump({'completed_region_ids': list(completed_rids)}, f)

    elapsed = time.time() - t0
    print(f"\nPhase 1 complete: {len(scoring_results):,} regions scored "
          f"in {elapsed:.1f}s")

    # ── Write Phase 1 results ────────────────────────────────────
    if scoring_results:
        # Flatten scoring results (drop per-track arrays for TSV)
        score_rows = []
        for s in scoring_results:
            row = {k: v for k, v in s.items() if k != 'te_track_means'}
            score_rows.append(row)
        score_df = pd.DataFrame(score_rows)
        score_path = args.output_dir / 'region_scores.tsv'
        score_df.to_csv(score_path, sep='\t', index=False)
        print(f"  Wrote {len(score_df):,} region scores to {score_path}")

    if ablation_results:
        abl_rows = []
        for a in ablation_results:
            row = {k: v for k, v in a.items() if k != 'track_deltas'}
            abl_rows.append(row)
        abl_df = pd.DataFrame(abl_rows)
        abl_path = args.output_dir / 'ablation_results.tsv'
        abl_df.to_csv(abl_path, sep='\t', index=False)
        print(f"  Wrote {len(abl_df):,} ablation results to {abl_path}")

    # ── Phase 2: Saturation mutagenesis (optional) ───────────────
    if args.satmut and scoring_results:
        print(f"\n{'='*60}")
        print(f"Phase 2: Saturation mutagenesis (top {args.top_n} candidates)")
        print(f"{'='*60}")

        # Rank by ablation delta (mask ablation)
        mask_abl = {a['region_id']: a for a in ablation_results
                    if a.get('ablation_type') == 'mask'}
        ranked = sorted(
            scoring_results,
            key=lambda s: mask_abl.get(s['region_id'], {}).get(
                'mean_abs_delta_all', 0
            ),
            reverse=True,
        )
        top_candidates = ranked[:args.top_n]

        satmut_results = []
        t1 = time.time()

        for ci, candidate in enumerate(top_candidates):
            rid = candidate['region_id']
            row = te_regions[te_regions['region_id'] == rid].iloc[0]

            if args.verbose:
                print(f"  [{ci+1}/{len(top_candidates)}] {rid} "
                      f"(te_bp={candidate['te_bp_strict']})")

            mask = masks_npz[rid] if rid in masks_npz.files else None
            if mask is None:
                continue

            chrom = str(row['chrom'])
            region_start = int(row['start'])
            region_end = int(row['end'])
            region_len = int(row['length'])

            seq, region_offset, pred_start, pred_end = build_genomic_window(
                genome, chrom, region_start, region_end, row['strand'],
                target_window=args.window_size,
            )
            if seq is None:
                continue

            effect_matrix, positions = run_saturation_mutagenesis(
                scorer, seq, region_offset, region_len, mask,
                pred_start, pred_end,
                batch_size=args.satmut_batch_size,
            )
            if effect_matrix is None:
                continue

            # Compute per-position max effect
            max_effects = effect_matrix.max(axis=1)

            satmut_entry = {
                'region_id': rid,
                'region_type': row['region_type'],
                'chrom': chrom,
                'start': region_start,
                'end': region_end,
                'fbgn': row['fbgn'],
                'symbol': row['symbol'],
                'n_positions_mutagenized': len(positions),
                'mean_max_effect': float(max_effects.mean()),
                'median_max_effect': float(np.median(max_effects)),
                'max_effect': float(max_effects.max()),
                'frac_high_effect': float(
                    (max_effects > np.percentile(max_effects, 90)).mean()
                ) if len(max_effects) > 10 else float('nan'),
            }
            satmut_results.append(satmut_entry)

            # Save per-position effect matrix
            np.save(
                str(args.output_dir / f'satmut_{rid}.npy'),
                effect_matrix,
            )
            # Save position map
            with open(args.output_dir / f'satmut_{rid}_positions.json', 'w') as f:
                json.dump({
                    'region_id': rid,
                    'positions_in_region': positions,
                    'nucleotide_order': NUCS,
                }, f)

        elapsed2 = time.time() - t1
        print(f"\nPhase 2 complete: {len(satmut_results)} regions mutagenized "
              f"in {elapsed2:.1f}s")

        if satmut_results:
            satmut_df = pd.DataFrame(satmut_results)
            satmut_path = args.output_dir / 'satmut_summary.tsv'
            satmut_df.to_csv(satmut_path, sep='\t', index=False)
            print(f"  Wrote satmut summary to {satmut_path}")

    # ── Phase 3: Integrated analysis ─────────────────────────────
    print(f"\n{'='*60}")
    print("Phase 3: Integrated analysis")
    print(f"{'='*60}")

    if not scoring_results:
        print("  No results to analyze.")
        return 0

    score_df = pd.DataFrame([
        {k: v for k, v in s.items() if k != 'te_track_means'}
        for s in scoring_results
    ])

    # Quality Paradox test: do regions with weaker BLAST signal (lower identity)
    # but strong regulatory scores show evidence of exaptation?
    # Merge with hits to get mean identity per region
    hits_df = pd.read_csv(args.dataset_dir / 'hits.tsv', sep='\t')
    region_identity = hits_df[hits_df['quality_tier'] == 'strict'].groupby(
        'region_id'
    )['pident'].mean().reset_index()
    region_identity.columns = ['region_id', 'mean_pident']

    merged = score_df.merge(region_identity, on='region_id', how='left')

    # Merge with ablation
    if ablation_results:
        mask_abl_df = pd.DataFrame([
            {k: v for k, v in a.items() if k != 'track_deltas'}
            for a in ablation_results if a.get('ablation_type') == 'mask'
        ])
        if len(mask_abl_df) > 0:
            mask_abl_df = mask_abl_df.rename(columns={
                'mean_delta_all': 'ablation_delta',
                'mean_abs_delta_all': 'ablation_abs_delta',
                'frac_tracks_decreased': 'ablation_frac_decreased',
            })
            merged = merged.merge(
                mask_abl_df[['region_id', 'ablation_delta',
                             'ablation_abs_delta', 'ablation_frac_decreased']],
                on='region_id', how='left',
            )

    # Classify candidates
    # High regulatory + low identity = exapted fossil candidate
    # High regulatory + high identity = young/active TE
    # Low regulatory + any identity = neutral insertion
    if 'ablation_abs_delta' in merged.columns and 'mean_pident' in merged.columns:
        valid = merged.dropna(subset=['ablation_abs_delta', 'mean_pident'])
        if len(valid) > 0:
            abl_median = valid['ablation_abs_delta'].median()
            ident_median = valid['mean_pident'].median()

            merged['candidate_class'] = 'unclassified'
            merged.loc[
                (merged['ablation_abs_delta'] >= abl_median) &
                (merged['mean_pident'] < ident_median),
                'candidate_class'
            ] = 'exapted_fossil'
            merged.loc[
                (merged['ablation_abs_delta'] >= abl_median) &
                (merged['mean_pident'] >= ident_median),
                'candidate_class'
            ] = 'young_active'
            merged.loc[
                merged['ablation_abs_delta'] < abl_median,
                'candidate_class'
            ] = 'neutral'

            print("\nCandidate classification (Quality Paradox test):")
            print(merged['candidate_class'].value_counts().to_string())

    # Write integrated results
    integrated_path = args.output_dir / 'integrated_results.tsv'
    merged.to_csv(integrated_path, sep='\t', index=False)
    print(f"\nWrote integrated results to {integrated_path}")

    # ── Write summary stats ──────────────────────────────────────
    stats = {
        'n_regions_scored': len(scoring_results),
        'n_ablations': len(ablation_results),
        'model': args.model,
        'window_size': args.window_size,
        'min_te_bp': args.min_te_bp,
        'ablation_reps': args.ablation_reps,
    }

    if 'te_enrichment' in score_df.columns:
        enrichment = score_df['te_enrichment'].dropna()
        stats['mean_te_enrichment'] = float(enrichment.mean())
        stats['median_te_enrichment'] = float(enrichment.median())
        stats['frac_te_enriched'] = float((enrichment > 1.0).mean())

    if 'ablation_abs_delta' in merged.columns:
        abd = merged['ablation_abs_delta'].dropna()
        stats['mean_ablation_delta'] = float(abd.mean())
        stats['median_ablation_delta'] = float(abd.median())

    if 'candidate_class' in merged.columns:
        class_counts = merged['candidate_class'].value_counts().to_dict()
        stats['candidate_classes'] = class_counts

    if args.satmut and 'satmut_results' in dir() and satmut_results:
        stats['n_satmut_regions'] = len(satmut_results)
        stats['satmut_top_n'] = args.top_n

    stats_path = args.output_dir / 'scoring_stats.json'
    with open(stats_path, 'w') as f:
        json.dump(stats, f, indent=2)

    # ── Print summary ────────────────────────────────────────────
    print(f"\n{'='*60}")
    print("Summary")
    print(f"{'='*60}")
    print(f"  Regions scored:        {stats['n_regions_scored']:,}")
    print(f"  Ablation tests:        {stats['n_ablations']:,}")
    if 'mean_te_enrichment' in stats:
        print(f"  Mean TE enrichment:    {stats['mean_te_enrichment']:.3f}")
        print(f"  Fraction TE-enriched:  {stats['frac_te_enriched']:.1%}")
    if 'mean_ablation_delta' in stats:
        print(f"  Mean ablation delta:   {stats['mean_ablation_delta']:.4f}")
    if 'candidate_classes' in stats:
        for cls, cnt in stats['candidate_classes'].items():
            print(f"  {cls}: {cnt}")
    print(f"\nOutput directory: {args.output_dir}")

    return 0


if __name__ == '__main__':
    # Tee all output to a persistent log file in the output directory
    import io

    class Tee:
        def __init__(self, *streams):
            self.streams = streams
        def write(self, data):
            for s in self.streams:
                s.write(data)
                s.flush()
        def flush(self):
            for s in self.streams:
                s.flush()

    # Parse args early to get output_dir for log placement
    # Default log goes next to the script if output_dir isn't created yet
    _default_log = Path(__file__).parent.parent / 'results' / 'ntv3_run.log'
    for i, arg in enumerate(sys.argv):
        if arg == '--output-dir' and i + 1 < len(sys.argv):
            _log_dir = Path(sys.argv[i + 1])
            _log_dir.mkdir(parents=True, exist_ok=True)
            _default_log = _log_dir / 'ntv3_run.log'
            break

    _default_log.parent.mkdir(parents=True, exist_ok=True)
    _log_file = open(_default_log, 'w')
    _tee = Tee(sys.stdout, _log_file)
    sys.stdout = _tee
    sys.stderr = Tee(sys.stderr, _log_file)

    try:
        rc = main()
    except Exception as e:
        print(f"\nFATAL: {e}")
        import traceback
        traceback.print_exc()
        rc = 1
    finally:
        _log_file.close()
        print(f"\nLog saved to: {_default_log}", file=sys.__stdout__)

    sys.exit(rc)
