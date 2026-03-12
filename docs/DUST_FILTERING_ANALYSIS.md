# DUST Filtering Analysis: Impact on TE Detection

**Date**: 2026-01-22
**Analysis**: Parameter sweep comparing DUST=yes vs DUST=no

## Executive Summary

DUST filtering, intended to mask low-complexity regions, removes **genuine TE-derived sequences** that contain microsatellite-like repeats intrinsic to certain TE families. Analysis of e-value distributions and sequence content demonstrates that DUST=no captures more high-quality hits with strong enrichment over shuffled controls.

**Key Finding**: DUST=no with stringent e-value thresholds (e < 0.001) provides better TE detection than DUST=yes.

## Methods

### Parameter Sweep
- **Sample**: 3,032 UTRs (10% of genome-wide 30,324 UTRs)
- **Controls**: Paired dinucleotide-shuffled sequences for each UTR
- **Parameters tested**: word_size (7, 11), gapopen (2, 5, 10), gapextend (1, 2), penalty (-1, -2), dust (yes, no)
- **Total runs**: 64 BLAST searches

### Key Comparisons
- DUST=yes vs DUST=no at multiple e-value thresholds
- Real vs shuffled enrichment ratios
- Positional distribution along UTRs by length

## Results

### 1. E-value Distribution Analysis

#### Hit Counts at Stringent E-value Thresholds (ws=7, optimal gap params)

| E-value | DUST=yes Real | DUST=no Real | Extra hits | DUST=yes Shuf | DUST=no Shuf |
|---------|---------------|--------------|------------|---------------|--------------|
| < 10    | 463,094       | 1,290,513    | +178%      | 197,241       | 239,763      |
| < 1     | 176,353       | 642,395      | +264%      | 47,395        | 62,920       |
| < 0.01  | 44,175        | 183,543      | +316%      | 3,523         | 5,551        |
| < 0.001 | 29,345        | 109,941      | +275%      | 887           | 1,482        |
| < 1e-5  | 17,642        | 50,186       | +184%      | 116           | 230          |
| < 1e-10 | 10,074        | 15,299       | +52%       | 0             | 0            |

**Critical observation**: At e < 1e-10 (very high quality), DUST=no captures **52% more hits** than DUST=yes, with **zero hits in shuffled controls** for both settings.

#### Real/Shuffled Enrichment Ratios

| E-value | DUST=yes Ratio | DUST=no Ratio |
|---------|----------------|---------------|
| < 10    | 2.35x          | 5.38x         |
| < 1     | 3.72x          | 10.21x        |
| < 0.01  | 12.54x         | 33.06x        |
| < 0.001 | 33.08x         | 74.18x        |
| < 1e-5  | 152.09x        | 218.20x       |
| < 1e-10 | inf            | inf           |

**DUST=no shows higher enrichment at ALL e-value thresholds**, contradicting the assumption that extra hits are noise.

### 2. Sequence Content of DUST-Filtered Hits

Analysis of 5,299 high-quality hits (e < 1e-10) unique to DUST=no:

#### Composition
| Hit Category | A% | T% | G% | C% | AT% |
|--------------|-----|-----|-----|-----|------|
| Unique to DUST=no | 40.5 | 14.2 | 17.6 | 27.7 | **54.7** |
| Shared (both settings) | 31.9 | 31.9 | 16.6 | 19.6 | **63.8** |

**Surprising finding**: Hits removed by DUST are **less AT-rich** (54.7%) than shared hits (63.8%). They contain structured repeats like CAG/GCA rather than simple AT runs.

#### Example Sequences Filtered by DUST

1. **CAG repeats** matching `roo` family TEs:
   ```
   Query:  CAGCAGCATCTGCAACATTAGCAACAGCAGCAGCAGCAACAACAGCAGCAACAACAACAGCAGCAGCAACAACAGCAAAT
   TE:     CAGCAACAACAGCAGCAGCAGAAACAGCAACAGCAGTAGCAACAGCAGCAACAACAACAGCAGCAGCAACAACGACGACG
   E-value: 5.27e-16, Length: 157bp, Identity: 70.7%
   ```

2. **ATAT patterns** matching `17.6` family TEs:
   ```
   Query:  ATATTAATAATTAATTA-TA--CGATGCATATATACATATATATATATAAATATATATGGATATGTATGTGTAAAACGGA
   TE:     ATATAAATAATAAATTAATAAGCGA-AAATTAAAACGTAT-TAAAAGTAAGAATA-ATAAATAAATAAGTGAAAATTCTA
   E-value: 7.42e-12, Length: 212bp, Identity: 66.0%
   ```

#### Top TE Families in DUST-Filtered Hits
| TE ID | Hits | Family |
|-------|------|--------|
| FBti0061482 | 227 | - |
| FBti0215166 | 175 | - |
| FBti0059769 | 89 | roo |
| FBti0059726 | 76 | roo |
| FBti0020154 | 70 | roo |

The hits are from **specific TE families**, not random low-complexity sequences.

### 3. Positional Distribution by UTR Length

#### All Hits (e < 10), Coarse Bins

| UTR Length | DUST | 5' Enrich | 3' Enrich |
|------------|------|-----------|-----------|
| <200bp     | yes  | 0.52x     | 1.52x     |
| <200bp     | no   | 0.50x     | 1.57x     |
| 200-500bp  | yes  | 0.62x     | 1.68x     |
| 200-500bp  | no   | 0.59x     | 1.81x     |
| 500-1kb    | yes  | 0.86x     | 1.31x     |
| 500-1kb    | no   | 1.14x     | 1.38x     |
| 1-2kb      | yes  | 1.06x     | 1.14x     |
| 1-2kb      | no   | 0.92x     | 1.01x     |
| >2kb       | yes  | 1.15x     | 1.05x     |
| >2kb       | no   | 1.13x     | 1.27x     |

#### Fine-Grained Positional Analysis (100bp UTR Length Bins, DUST=no)

Each UTR normalized to relative position 0-100%. Enrichment = Real density / Shuffled density.

| UTR Length | 5' (0-30%) | Mid (30-70%) | 3' (70-100%) | Pattern |
|------------|-----------|-------------|-------------|---------|
| 0-100bp    | 0.49x     | 1.06x       | 1.26x       | 3' bias |
| 100-200bp  | 0.50x     | 1.00x       | 1.53x       | 3' bias |
| 200-300bp  | 0.56x     | 0.85x       | 1.66x       | 3' bias |
| 300-400bp  | 0.79x     | 0.70x       | 1.60x       | 3' bias |
| 400-500bp  | 0.57x     | 0.84x       | 1.67x       | 3' bias |
| 500-600bp  | 0.85x     | 1.03x       | 1.14x       | 3' bias |
| 600-700bp  | **1.46x** | 0.63x       | **1.28x**   | **U-shape** |
| 700-800bp  | 0.86x     | 0.97x       | 1.17x       | 3' bias |
| 800-900bp  | **1.44x** | 0.76x       | 1.08x       | 5' bias |
| 900-1000bp | **1.43x** | 0.66x       | **1.65x**   | **U-shape** |
| 1000-1100bp| 0.49x     | 1.03x       | **2.10x**   | Strong 3' bias |
| 1100-1200bp| 0.53x     | 1.01x       | **1.49x**   | 3' bias |
| 1200-1300bp| **1.13x** | 1.54x       | 0.62x       | 5' bias |
| 1300-1400bp| **2.07x** | 0.94x       | 0.60x       | 5' bias |
| 1400-1500bp| 0.97x     | 0.79x       | **1.41x**   | 3' bias |

**Key patterns by UTR length:**
- **Short UTRs (0-500bp)**: Consistent 3' bias (1.3-1.7x). TE hits concentrate near the poly-A site. 5' end depleted (0.5-0.8x).
- **Medium UTRs (600-1000bp)**: U-shape emerges at 600-700bp and 900-1000bp — both ends enriched, middle depleted. This suggests functional constraints against TE retention in the UTR interior (regulatory elements, miRNA sites).
- **Longer UTRs (>1000bp)**: Pattern becomes variable. Some bins show strong 3' bias (1000-1200bp), others show 5' bias (1200-1400bp). May reflect different TE insertion/retention dynamics at greater distances from both the stop codon and poly-A site.

**Shuffled controls show approximately uniform distribution across all UTR length bins**, confirming these positional patterns are biological signal, not compositional artifacts.

#### High-Quality Hits (e < 0.01)

| UTR Length | DUST | 5' Enrich | 3' Enrich |
|------------|------|-----------|-----------|
| <200bp     | yes  | 3.79x     | 7.46x     |
| <200bp     | no   | 0.94x     | 6.72x     |
| 500-1kb    | yes  | 0.85x     | 2.87x     |
| 500-1kb    | no   | 1.21x     | 4.40x     |
| >2kb       | yes  | 1.09x     | 1.07x     |
| >2kb       | no   | 1.31x     | 1.32x     |

**At stringent e-values, DUST=no shows stronger 3' enrichment** in most UTR length categories, and reveals a U-shaped distribution in medium/long UTRs that DUST=yes obscures.

## Conclusions

1. **DUST filtering removes genuine TE hits**: The sequences filtered by DUST are microsatellite-like repeats (CAG, ATAT) that are intrinsic to specific TE families (especially `roo`).

2. **Higher signal-to-noise with DUST=no**: Real/shuffled enrichment ratios are higher with DUST=no at all e-value thresholds, indicating the extra hits are biological signal, not noise.

3. **Composition evidence**: DUST-filtered hits are less AT-rich (54.7%) than shared hits (63.8%), containing GC-rich structured repeats.

4. **Positional bias preserved**: The 3' enrichment pattern is present (and often stronger) with DUST=no.

5. **Recommended approach**: Use **DUST=no with e-value < 0.001** for TE fossil detection. This captures more genuine TE-derived sequences while maintaining high specificity through stringent e-value filtering.

## Files Generated

- `results/param_sweep_full/` - All 64 BLAST result files (3,032 UTR sample)
- `figures/param_sweep_analysis/evalue_analysis.png` - E-value distribution comparison
- `figures/param_sweep_analysis/position_by_utr_length_dust_comparison.png` - Positional analysis (coarse bins)
- `figures/param_sweep_analysis/position_by_100bp_bins_clean.png` - Positional analysis (100bp UTR bins, decile position)
- `figures/param_sweep_analysis/position_normalized_by_utr_length.png` - Full 4-panel normalized position analysis
- `figures/param_sweep_analysis/param_effects_summary.png` - Parameter effects overview
- `figures/param_sweep_analysis/wordsize_comparison.png` - Word size 7 vs 11
- `figures/param_sweep_analysis/te_family_analysis.png` - TE family enrichment
- `figures/param_sweep_analysis/analysis_summary.tsv` - Numerical summary

## Parameter Recommendations

Based on this analysis:

```
word_size=7
gapopen=2
gapextend=1
penalty=-1
reward=1
dust=no          # Changed from previous recommendation
evalue < 0.001   # Apply stringent post-filtering
```
