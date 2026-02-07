# Running the NTv3 Pipeline on Lightning.ai

## Automated launch script

```bash
pip install lightning-sdk
lightning login

# Test mode (50 regions)
bash scripts/launch_ntv3_lightning.sh --test

# Full run
bash scripts/launch_ntv3_lightning.sh
```

The script creates a studio, transfers data, installs deps, runs tests, and launches the pipeline. Edit the `STUDIO_NAME`, `TEAMSPACE`, and `MACHINE` variables at the top if needed.

**NOTE:** You may need to adjust the `SSH_TARGET` variable depending on how your Lightning.ai SSH access is configured. Check your studio settings for the SSH hostname.

---

## What needs to get there

| What | Size | Source | How |
|------|------|--------|-----|
| Scripts + docs | ~2 MB | Git-tracked | `git clone` |
| `results/te_fossil_lm_dataset/` | 691 MB | Gitignored (local only) | `rsync` or `scp` |
| `data/references/dmel_genome.fasta` | 139 MB | Gitignored (local only) | `rsync` or `scp` |
| NTv3 model weights | ~400 MB | HuggingFace | Auto-downloaded on first run |
| **Total upload** | **~830 MB** | | |

**Key insight:** Almost all your data is gitignored. `git clone` gives you scripts and docs only. The genome FASTA, regulatory region data, and prepared dataset must be transferred separately.

## Step-by-step

### 1. Provision a Lightning.ai studio

- GPU: **A100 40GB** (minimum) or A100 80GB (comfortable)
- Storage: 50GB is plenty
- Image: Any PyTorch image works; you'll pip install the rest

### 2. Clone the repo (gets scripts only)

```bash
# On the Lightning.ai machine
git clone git@github.com:jacobboysen/repeat_finder.git
cd repeat_finder
```

### 3. Transfer the data

**Option A: rsync over SSH (recommended)**

Lightning.ai studios expose an SSH endpoint. Find it in the studio settings.

```bash
# On your LOCAL machine
# Replace LIGHTNING_SSH with your studio's SSH address

# Transfer the prepared LM dataset (691 MB)
rsync -avz --progress \
  results/te_fossil_lm_dataset/ \
  LIGHTNING_SSH:~/repeat_finder/results/te_fossil_lm_dataset/

# Transfer the genome FASTA (139 MB)
rsync -avz --progress \
  data/references/dmel_genome.fasta \
  LIGHTNING_SSH:~/repeat_finder/data/references/dmel_genome.fasta
```

**Option B: tar + scp (simpler)**

```bash
# On your LOCAL machine
tar czf te_fossil_data.tar.gz \
  results/te_fossil_lm_dataset/ \
  data/references/dmel_genome.fasta

scp te_fossil_data.tar.gz LIGHTNING_SSH:~/repeat_finder/

# On the REMOTE machine
cd ~/repeat_finder
tar xzf te_fossil_data.tar.gz
rm te_fossil_data.tar.gz
```

**Option C: Lightning.ai persistent storage**

If your studio has a `/teamspace/` or persistent volume mounted, you can upload there once and symlink across studio restarts:

```bash
# Upload once to persistent storage
rsync -avz results/te_fossil_lm_dataset/ /teamspace/te_fossil_lm_dataset/
rsync -avz data/references/dmel_genome.fasta /teamspace/dmel_genome.fasta

# Then in each studio session:
cd ~/repeat_finder
mkdir -p results data/references
ln -s /teamspace/te_fossil_lm_dataset results/te_fossil_lm_dataset
ln -s /teamspace/dmel_genome.fasta data/references/dmel_genome.fasta
```

### 4. Install dependencies

```bash
pip install torch transformers pyfaidx numpy pandas
```

Verify:
```bash
python -c "import torch; print(f'CUDA: {torch.cuda.is_available()}, Device: {torch.cuda.get_device_name(0)}')"
python -c "from transformers import AutoModel; print('transformers OK')"
```

### 5. Run tests

```bash
python scripts/test_ntv3_pipeline.py
```

This will:
- Download NTv3 model weights from HuggingFace (~400MB, cached after first run)
- Verify Drosophila species key works
- Run inference on test sequences
- Run the full pipeline on 5 regions

### 6. Run the real pipeline

```bash
# Quick validation (50 regions, ~10 min)
python scripts/score_te_fossils_ntv3.py -v --max-regions 50

# Full run, all TE-bearing regions (~6-12 hours on A100)
python scripts/score_te_fossils_ntv3.py -v

# Full run + saturation mutagenesis on top 100 (~12-24 hours)
python scripts/score_te_fossils_ntv3.py -v --satmut --top-n 100

# Use 650M model for higher accuracy (needs A100 80GB)
python scripts/score_te_fossils_ntv3.py -v --model InstaDeepAI/NTv3_650M_post
```

**The pipeline checkpoints every 50 regions.** If your studio times out or you need to stop, just re-run the same command — it picks up where it left off.

### 7. Get results back

```bash
# On your LOCAL machine
rsync -avz --progress \
  LIGHTNING_SSH:~/repeat_finder/results/ntv3_te_scoring/ \
  results/ntv3_te_scoring/
```

## Quick reference: files on the GPU box

```
repeat_finder/
├── scripts/
│   ├── score_te_fossils_ntv3.py    # Main pipeline
│   ├── test_ntv3_pipeline.py       # Tests
│   ├── prepare_te_fossils_for_lm.py
│   └── utils/                      # Shared utilities
├── data/
│   └── references/
│       └── dmel_genome.fasta       # 139 MB — transferred manually
├── results/
│   ├── te_fossil_lm_dataset/       # 691 MB — transferred manually
│   │   ├── regions.tsv
│   │   ├── hits.tsv
│   │   ├── sequences.fasta
│   │   ├── masks.npz
│   │   └── dataset_stats.json
│   └── ntv3_te_scoring/            # OUTPUT — transfer back when done
│       ├── region_scores.tsv
│       ├── ablation_results.tsv
│       ├── integrated_results.tsv
│       ├── scoring_stats.json
│       └── satmut_*.npy            # (if --satmut)
└── ~/.cache/huggingface/           # Model weights (auto-downloaded)
```

## Troubleshooting

**"CUDA out of memory"**: Use `--window-size 16384` to halve the context window, or use the 100M model instead of 650M.

**Studio timeout mid-run**: Just re-run the same command. Checkpoint resume is automatic.

**"Could not find Drosophila in species map"**: The script prints all available species keys. Check the output and update if InstaDeep changed the key name.

**Slow first run**: The first forward pass downloads and caches model weights (~400MB). Subsequent runs start immediately.
