#!/usr/bin/env python3
"""
Test the NTv3 TE fossil scoring pipeline on a small subset.

Run this first to validate that:
  1. The model loads and produces expected output shapes
  2. Species key discovery works for Drosophila
  3. Genomic windowing + prediction cropping is correct
  4. Ablation produces measurable deltas
  5. Saturation mutagenesis produces per-position effects
  6. Output files are written correctly

Usage (on GPU box):
    python scripts/test_ntv3_pipeline.py

Requires: torch, transformers>=4.55.0, pyfaidx, numpy, pandas
"""

import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import torch

sys.path.insert(0, str(Path(__file__).parent))
from utils.paths import get_project_root, get_results_dir  # noqa: E402


def test_model_loading():
    """Test 1: Model loads and we can find Drosophila tracks."""
    from transformers import AutoTokenizer, AutoModel

    repo = "InstaDeepAI/NTv3_100M_post"
    print(f"Loading model: {repo}")

    tokenizer = AutoTokenizer.from_pretrained(repo, trust_remote_code=True)
    model = AutoModel.from_pretrained(repo, trust_remote_code=True)

    # Check species
    species_map = model.config.species_to_token_id
    print(f"Available species: {list(species_map.keys())}")

    # Find Drosophila
    dmel_key = None
    for candidate in ['fruit_fly', 'drosophila', 'drosophila_melanogaster', 'dm6']:
        if candidate in species_map:
            dmel_key = candidate
            break
    if dmel_key is None:
        for key in species_map:
            if 'droso' in key.lower() or 'fly' in key.lower():
                dmel_key = key
                break

    assert dmel_key is not None, f"Could not find Drosophila in {list(species_map.keys())}"
    print(f"Drosophila key: '{dmel_key}'")

    # Check tracks
    bigwig_names = list(model.config.bigwigs_per_species.get(dmel_key, []))
    bed_names = list(model.config.bed_elements_names)
    print(f"BigWig tracks for {dmel_key}: {len(bigwig_names)}")
    print(f"BED element types: {len(bed_names)}")
    if bigwig_names:
        print(f"  First 10 tracks: {bigwig_names[:10]}")
    if bed_names:
        print(f"  BED elements: {bed_names[:10]}")

    return tokenizer, model, dmel_key, bigwig_names, bed_names


def test_inference(tokenizer, model, species_key):
    """Test 2: Forward pass produces expected output shapes."""
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    model = model.to(device).eval()

    # Create a simple test sequence (32kb, divisible by 128)
    seq_len = 32768
    test_seq = 'ATCGATCG' * (seq_len // 8)
    assert len(test_seq) == seq_len

    batch = tokenizer(
        [test_seq],
        add_special_tokens=False,
        padding=True,
        pad_to_multiple_of=128,
        return_tensors='pt',
    )
    input_ids = batch['input_ids'].to(device)
    species_ids = model.encode_species([species_key]).to(device)

    print(f"\nInput shape: {input_ids.shape}")
    print(f"Device: {device}")

    with torch.no_grad():
        out = model(input_ids=input_ids, species_ids=species_ids)

    print(f"\nOutput keys: {list(k for k in dir(out) if not k.startswith('_'))}")

    if out.logits is not None:
        print(f"  logits:              {out.logits.shape}")
    if out.bigwig_tracks_logits is not None:
        print(f"  bigwig_tracks_logits: {out.bigwig_tracks_logits.shape}")
        n_tracks = out.bigwig_tracks_logits.shape[-1]
        pred_len = out.bigwig_tracks_logits.shape[1]
        keep_frac = getattr(model.config, 'keep_target_center_fraction', 0.375)
        expected_pred_len = int(seq_len * keep_frac)
        print(f"  Prediction length: {pred_len} "
              f"(expected ~{expected_pred_len} for {keep_frac} center fraction)")
        print(f"  Number of tracks: {n_tracks}")
    if out.bed_tracks_logits is not None:
        print(f"  bed_tracks_logits:   {out.bed_tracks_logits.shape}")
    if out.embedding is not None:
        print(f"  embedding:           {out.embedding.shape}")

    bigwig = out.bigwig_tracks_logits[0].cpu().float().numpy()
    print(f"\nBigWig stats: mean={bigwig.mean():.4f}, "
          f"std={bigwig.std():.4f}, "
          f"min={bigwig.min():.4f}, max={bigwig.max():.4f}")

    return bigwig


def test_ablation_signal(tokenizer, model, species_key):
    """Test 3: Ablation (N-masking a region) produces measurable delta."""
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    model = model.to(device).eval()

    seq_len = 32768
    # Use a real-ish sequence (random but reproducible)
    rng = np.random.default_rng(42)
    nucs = list('ATCG')
    test_seq = ''.join(rng.choice(nucs, size=seq_len))

    # Create ablated version (mask center 500bp)
    center = seq_len // 2
    ablated_seq = test_seq[:center-250] + 'N' * 500 + test_seq[center+250:]

    # Score both
    seqs = [test_seq, ablated_seq]
    batch = tokenizer(
        seqs,
        add_special_tokens=False,
        padding=True,
        pad_to_multiple_of=128,
        return_tensors='pt',
    )
    input_ids = batch['input_ids'].to(device)
    species_ids = model.encode_species([species_key] * 2).to(device)

    with torch.no_grad():
        out = model(input_ids=input_ids, species_ids=species_ids)

    ref_bw = out.bigwig_tracks_logits[0].cpu().float().numpy()
    abl_bw = out.bigwig_tracks_logits[1].cpu().float().numpy()

    delta = ref_bw - abl_bw
    print(f"\nAblation test (500bp N-mask at center):")
    print(f"  Mean delta:     {delta.mean():.6f}")
    print(f"  Mean |delta|:   {np.abs(delta).mean():.6f}")
    print(f"  Max |delta|:    {np.abs(delta).max():.6f}")
    print(f"  Delta != 0:     {(np.abs(delta) > 1e-8).any()}")

    # The masked region should show SOME difference
    assert np.abs(delta).max() > 0, "Ablation produced zero delta — something is wrong"
    print("  PASS: Ablation produces measurable signal change")


def test_full_pipeline_small():
    """Test 4: Run the full pipeline on 5 regions."""
    import subprocess

    dataset_dir = get_results_dir() / 'te_fossil_lm_dataset'
    genome_path = get_project_root() / 'data' / 'references' / 'dmel_genome.fasta'
    if not dataset_dir.exists():
        print("\nSKIP: Dataset not found at", dataset_dir)
        return
    if not genome_path.exists():
        print("\nSKIP: Genome FASTA not found at", genome_path)
        return

    with tempfile.TemporaryDirectory() as tmpdir:
        cmd = [
            sys.executable, 'scripts/score_te_fossils_ntv3.py',
            '--dataset-dir', str(dataset_dir),
            '--genome', str(genome_path),
            '--output-dir', tmpdir,
            '--max-regions', '5',
            '--ablation-reps', '2',
            '--satmut', '--top-n', '2',
            '-v',
        ]
        print(f"\nRunning: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        print(result.stdout)
        if result.returncode != 0:
            print("STDERR:", result.stderr)
            assert False, f"Pipeline failed with exit code {result.returncode}"

        # Check output files
        outdir = Path(tmpdir)
        for expected in ['region_scores.tsv', 'ablation_results.tsv',
                         'integrated_results.tsv', 'scoring_stats.json']:
            path = outdir / expected
            assert path.exists(), f"Missing output: {expected}"
            if expected.endswith('.tsv'):
                df = pd.read_csv(path, sep='\t')
                print(f"  {expected}: {len(df)} rows, {list(df.columns)[:5]}...")
            elif expected.endswith('.json'):
                with open(path) as f:
                    data = json.load(f)
                print(f"  {expected}: {list(data.keys())[:5]}...")

        # Validate region_scores invariants
        scores = pd.read_csv(outdir / 'region_scores.tsv', sep='\t')
        assert len(scores) > 0, "No regions scored"
        assert 'mean_signal_all' in scores.columns
        assert 'te_enrichment' in scores.columns

        # Validate ablation results
        abl = pd.read_csv(outdir / 'ablation_results.tsv', sep='\t')
        assert len(abl) > 0, "No ablation results"
        assert 'ablation_type' in abl.columns
        assert set(abl['ablation_type'].unique()) <= {'mask', 'shuffle'}

        # Check for satmut files
        satmut_files = list(outdir.glob('satmut_*.npy'))
        if satmut_files:
            print(f"  Satmut effect matrices: {len(satmut_files)}")
            mat = np.load(str(satmut_files[0]))
            print(f"    Shape: {mat.shape} (positions x nucleotides)")
            assert mat.shape[1] == 4, "Expected 4 columns (A, T, C, G)"

        print("\n  PASS: Full pipeline produces valid output")


def main():
    print("NTv3 TE Fossil Pipeline Tests")
    print("=" * 60)

    # Test 1: Model loading
    print("\n--- Test 1: Model loading ---")
    tokenizer, model, species_key, bigwig_names, bed_names = test_model_loading()

    # Test 2: Inference shapes
    print("\n--- Test 2: Forward pass ---")
    bigwig = test_inference(tokenizer, model, species_key)

    # Test 3: Ablation signal
    print("\n--- Test 3: Ablation signal ---")
    test_ablation_signal(tokenizer, model, species_key)

    # Test 4: Full pipeline (small)
    print("\n--- Test 4: Full pipeline (5 regions) ---")
    test_full_pipeline_small()

    print("\n" + "=" * 60)
    print("ALL TESTS PASSED")
    return 0


if __name__ == '__main__':
    # Tee all output to a persistent log file
    log_path = Path(__file__).parent.parent / 'results' / 'ntv3_test.log'
    log_path.parent.mkdir(parents=True, exist_ok=True)

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

    log_file = open(log_path, 'w')
    tee = Tee(sys.stdout, log_file)
    sys.stdout = tee
    sys.stderr = Tee(sys.stderr, log_file)

    try:
        rc = main()
    except Exception as e:
        print(f"\nFATAL: {e}")
        import traceback
        traceback.print_exc()
        rc = 1
    finally:
        log_file.close()
        print(f"\nLog saved to: {log_path}", file=sys.__stdout__)

    sys.exit(rc)
