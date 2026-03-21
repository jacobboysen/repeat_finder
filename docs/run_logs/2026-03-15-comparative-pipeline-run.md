# Comparative TE Pipeline Run Log — 2026-03-15

## Overview

First run of the dmel x dpse comparative TE fossil pipeline. Produced 4 BLAST
result sets (2 self + 2 cross) plus ortholog lookup table.

**Spec:** `docs/superpowers/specs/2026-03-14-comparative-te-pipeline-design.md`

## Environment

- **Machine:** macOS Darwin 24.6.0, 16 cores
- **Conda env:** `fossil-finder` (Python 3.13, Snakemake 9.x, BLAST 2.16, RepeatMasker 4.2.2)
- **Package:** `fossil_finder 0.1.0` (editable install)

## Input Data

### D. melanogaster (dm6)
- **Genome:** `data/references/dmel_genome.fasta` (146 MB, FlyBase r6.66)
- **Annotation:** `data/references/dmel_annotation_coding.gff3` (37 MB, FlyBase r6.66)
- **RepeatMasker .out:** `data/references/dm6.fa.out` (17 MB, pre-existing, 137,555 annotations)
- **Config:** `src/fossil_finder/config/genomes/dmel_r6.66_rmmode.yaml`

### D. pseudoobscura (UCI_Dpse_MV25)
- **Genome:** `data/references/dpse_genome.fasta` (165 MB, NCBI GCF_009870125.1)
  - Downloaded via `datasets download genome accession GCF_009870125.1`
  - 5 chromosomes (NC_046679.1–NC_046683.1) + 65 unplaced scaffolds
- **Annotation:** `data/references/dpse_annotation.gff3` (109 MB, NCBI RefSeq Gnomon Jan 2025)
  - No explicit UTR features; 3'UTRs inferred via `infer_utrs()` from mRNA minus CDS
  - 25,693 inferred 3'UTRs (median 245 bp, mean 590 bp)
- **RepeatMasker .out:** `data/references/dpse_genome.fa.out` (28 MB, generated this session)
  - Run with: `RepeatMasker -e rmblast -pa 8 -lib data/references/dmel_te_consensus.fasta -qq -no_is`
  - Used dmel TE consensus as custom library (Dfam partition not available locally)
  - 216,164 annotations
- **Config:** `src/fossil_finder/config/genomes/dpse_MV25.yaml`

### Ortholog Table
- **Source:** NCBI `gene_orthologs.gz` (FTP bulk download, 116 MB compressed)
- **Filter:** tax_id 7227 (dmel) ↔ 7237 (dpse)
- **Output:** `data/references/orthologs/dmel_dpse_orthologs.tsv`
  - 11,515 ortholog pairs
  - Columns: dmel_fbgn, dmel_entrez, dmel_symbol, dpse_entrez, dpse_symbol, relationship
  - FBgn↔EntrezGene mapping derived from dmel GFF3 Dbxref attributes

## BLAST Parameters

```yaml
word_size: 7
gapopen: 2
gapextend: 1
penalty: -1
reward: 1
dust: false
evalue: 0.01
max_target_seqs: 10000
```

## Pipeline Execution

### TE Library Construction (RepeatMasker → TE instances)

Both genomes: `.out` file → `build_te_fasta_from_repeatmasker()` → TE instance FASTA → `makeblastdb`

| | dm6 | dpse |
|---|---|---|
| RM annotations | 137,555 | 216,164 |
| TE instances extracted | 119,240 | 199,327 |
| TE FASTA size | ~24 MB | ~14 MB |

**Note:** Both TE databases include Simple_repeat and Low_complexity entries (see Known Issues).

### 3'UTR Extraction

| | dm6 | dpse |
|---|---|---|
| UTR source | Explicit GFF3 `three_prime_UTR` | Inferred from mRNA minus CDS |
| UTR count | 28,539 | 25,693 |
| Total UTR bp | ~11 Mb | ~15 Mb |

### Snakemake Execution

```bash
conda run -n fossil-finder snakemake --snakefile workflows/Snakefile --cores 8 \
  --config 'genome_configs=["...dmel_r6.66_rmmode.yaml","...dpse_MV25.yaml"]' \
  base_dir=data/ output_dir=results/comparative/ \
  --omit-from analyze
```

14 jobs total (2 extract + 1 RM + 2 build_te + 2 makeblastdb + 2 self-BLAST + 2 cross-BLAST + 2 analyze).
Analyze steps skipped via `--omit-from` due to performance bugs (now fixed, not yet re-run).

## Output Files

### Self-BLAST Results

| Run | Output | Size | Lines |
|-----|--------|------|-------|
| dmel self | `results/comparative/dm6/blast_results.tsv` | 1.6 GB | 7,262,695 |
| dpse self | `results/comparative/UCI_Dpse_MV25/blast_results.tsv` | 5.5 GB | 25,282,100 |

### Cross-BLAST Results

| Run | Output | Size | Lines |
|-----|--------|------|-------|
| dmel UTRs vs dpse TEs | `results/comparative/dm6_vs_UCI_Dpse_MV25_te/blast_results.tsv` | 2.6 GB | 12,783,919 |
| dpse UTRs vs dmel TEs | `results/comparative/UCI_Dpse_MV25_vs_dm6_te/blast_results.tsv` | 2.3 GB | 9,972,287 |

### Intermediate Files

| File | Description |
|------|-------------|
| `results/comparative/dm6/regions.fa` | dmel 3'UTR FASTA (12 MB) |
| `results/comparative/dm6/regions.tsv` | dmel region metadata (28,539 regions) |
| `results/comparative/dm6/repeatmasker/te_instances.fasta` | dmel RM-derived TEs (119,240 seqs) |
| `results/comparative/dm6/blastdb/te_db.*` | dmel BLAST database |
| `results/comparative/UCI_Dpse_MV25/regions.fa` | dpse 3'UTR FASTA |
| `results/comparative/UCI_Dpse_MV25/regions.tsv` | dpse region metadata (25,693 regions) |
| `results/comparative/UCI_Dpse_MV25/repeatmasker/te_instances.fasta` | dpse RM-derived TEs (199,327 seqs) |
| `results/comparative/UCI_Dpse_MV25/blastdb/te_db.*` | dpse BLAST database |

## Known Issues

### 1. Simple_repeat/Low_complexity inflation (CRITICAL)

The TE instance databases include Simple_repeat and Low_complexity sequences from
RepeatMasker, which massively inflate BLAST hit counts. Post-hoc analysis shows:

| | dm6 | dpse |
|---|---|---|
| Total BLAST hits | 7.3M | 25.3M |
| Simple_repeat hits | 6.1M (83%) | 22.5M (89%) |
| **Real TE hits** | **1.1M** | **1.1M** |

After excluding Simple_repeat/Low_complexity, both species produce ~1.1M real TE hits.

**Recommendation:** Filter `build_te_fasta_from_repeatmasker()` output to exclude
`class=Simple_repeat` and `class=Low_complexity` before building BLAST databases.
This would reduce dpse database from 199k to ~10k sequences and eliminate ~96% of
uninformative hits.

### 2. Analyze step not run

The `analyze` rule was skipped due to O(N²) performance bugs in `find_overlaps()` and
`classify_hits()` (see Code Changes below). These are now fixed but the analyze step
has not been re-run.

### 3. Duplicate dpse UTR coordinates

NCBI Gnomon annotation produces multiple transcript isoforms sharing the same UTR
coordinates. Of 25,693 dpse UTR regions, only ~16,029 have unique genomic coordinates
(~9,664 duplicates). Consider deduplication in future runs.

## Code Changes Made This Session

| File | Change |
|------|--------|
| `src/fossil_finder/io/gff.py` | Added `infer_utrs()` + auto-inference in `parse_gff3()` |
| `src/fossil_finder/repeatmasker/library.py` | Fixed chr prefix mismatch (chr2L vs 2L) |
| `src/fossil_finder/analysis/repeatmasker.py` | Vectorized `find_overlaps()` and `classify_hits()` |
| `workflows/Snakefile` | Added `cross_blast` rule, `CROSS_PAIRS`, removed `-parse_seqids` |
| `src/fossil_finder/config/genomes/dpse_MV25.yaml` | New dpse config |
| `src/fossil_finder/config/genomes/dmel_r6.66_rmmode.yaml` | New dmel RM-mode config |

## Re-run: Filtered TE Databases (2026-03-15 morning)

Added `exclude_simple=True` to `build_te_fasta_from_repeatmasker()`, filtering out
Simple_repeat, Low_complexity, and Satellite entries before BLAST database construction.

### Filtered TE database sizes

| | Before | After | Removed |
|---|---|---|---|
| dmel | 119,240 seqs | 28,258 seqs | 90,982 (76%) |
| dpse | 199,327 seqs | 9,949 seqs | 189,378 (95%) |

### Filtered BLAST results

| Run | Before | After | Reduction |
|-----|--------|-------|-----------|
| dmel self | 1.6 GB / 7.3M | 325 MB / 1.1M | 6.6x |
| dpse self | 5.5 GB / 25.3M | 265 MB / 959K | 26x |
| dmel→dpse cross | 2.6 GB / 12.8M | 214 MB / 913K | 14x |
| dpse→dmel cross | 2.3 GB / 10.0M | 366 MB / 1.1M | 9x |

Saturated queries (>=10,000 hits): 4 dmel + 1 dpse (vs 177 + 719 before).

Total pipeline time with filtering: ~80 minutes (vs ~6 hours before).

## Wall Clock Times (initial run, with simple repeats)

| Step | Duration |
|------|----------|
| dpse genome download | ~2 min |
| dpse RepeatMasker | ~15 min (rush mode, 8 threads) |
| dmel TE extraction | ~30 sec |
| dpse TE extraction | ~2 min |
| dmel self-BLAST | ~85 min |
| dpse self-BLAST | ~120 min |
| dmel→dpse cross-BLAST | ~60 min |
| dpse→dmel cross-BLAST | ~90 min |
| **Total pipeline** | **~6 hours** |
