# Conservation Analysis Workflow

This directory contains scripts for validating TE hits using phyloP conservation scores from the UCSC 27-way Drosophila alignment.

## Overview

The conservation analysis compares phyloP scores between TE hit regions and UTR background to assess whether hits represent functional sequences under selective constraint.

## Prerequisites

### Data Files

1. **phyloP 27-way conservation scores** (396 MB):
   ```bash
   cd data/references
   curl -L -O https://hgdownload.soe.ucsc.edu/goldenPath/dm6/phyloP27way/dm6.phyloP27way.bw
   ```

2. **UCSC bigWigAverageOverBed tool**:
   ```bash
   mkdir -p data/references/tools
   cd data/references/tools
   curl -L -O http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/bigWigAverageOverBed
   chmod +x bigWigAverageOverBed
   ```

### Input Files Required

- `results/repeatmasker_analysis/blast_hits_novel.tsv` - BLAST hits not in RepeatMasker
- `results/repeatmasker_analysis/blast_hits_known_repeatmasker.tsv` - BLAST hits overlapping RepeatMasker
- `data/references/dmel-all-r6.66.gff` - FlyBase genome annotations

## Workflow Steps

### Step 1: Convert TE Hits to Genomic BED Format

```bash
python scripts/convert_all_hits_to_bed.py
```

This converts UTR-relative BLAST coordinates to genomic coordinates:
- Plus strand: `genomic_pos = UTR_start + query_pos - 1`
- Minus strand: `genomic_pos = UTR_end - query_pos + 1`
- Adds "chr" prefix for UCSC compatibility

**Output**: `results/repeatmasker_analysis/te_hits_all_genomic.bed`

### Step 2: Extract Conservation Scores

```bash
./data/references/tools/bigWigAverageOverBed \
    data/references/dm6.phyloP27way.bw \
    results/repeatmasker_analysis/te_hits_all_genomic.bed \
    results/repeatmasker_analysis/te_hits_all_conservation.tab
```

**Output format** (6 columns):
```
name    size    covered    sum    mean0    mean
```

- `mean0`: Mean including bases with no data (as 0)
- `mean`: Mean over bases with data only (use this one)

### Step 3: Analyze Conservation Distribution

```bash
python scripts/analyze_all_hits_conservation.py
```

Produces summary statistics by:
- Conservation category (highly conserved, conserved, weak, fast-evolving)
- Quality threshold (identity, length)
- Novel vs Known (RepeatMasker)
- TE family

### Step 4: Generate Visualizations

```bash
conda activate bioinformatics-program
python scripts/plot_all_vs_hq_conservation.py
```

**Outputs**:
- `figures/repeatmasker_comparison/15_all_vs_hq_conservation.png`
- `figures/repeatmasker_comparison/16_conservation_by_quality_tier.png`

## Interpreting phyloP Scores

| Score | Interpretation |
|-------|---------------|
| > 2 | Highly conserved (strong purifying selection) |
| 1-2 | Conserved (under constraint) |
| 0-1 | Weakly conserved (weak constraint or neutral) |
| < 0 | Fast-evolving (positive selection or relaxed) |

## Key Finding

**Counter-intuitively, stricter quality filters yield LESS conserved hits:**

| Filter | N hits | % Conserved |
|--------|--------|-------------|
| All hits | 1.8M | 61% |
| ≥80% id | 468K | 54% |
| ≥80%/≥50bp | 16K | 25% |
| ≥80%/≥100bp | 6K | 8% |

**Interpretation**: Long, high-identity hits are more likely to be recent TE insertions evolving neutrally. Short hits in conserved regions may represent ancient domesticated TE sequences or convergent evolution.

## Files Generated

| File | Description |
|------|-------------|
| `te_hits_all_genomic.bed` | All hits in genomic coordinates |
| `te_hits_all_conservation.tab` | Conservation scores for all hits |
| `te_hits_conservation_hq.tab` | Conservation scores for HQ hits only |

## References

- [UCSC phyloP 27-way](https://hgdownload.soe.ucsc.edu/goldenPath/dm6/phyloP27way/)
- phyloP method: Pollard et al. (2010) Genome Research
