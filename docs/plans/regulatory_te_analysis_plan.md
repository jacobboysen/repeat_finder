# Regulatory Region TE Analysis Plan

## Overview
Extend the TE fossil mining pipeline to analyze transposable element content in regulatory regions: promoters, enhancers, silencers, and transcription factor binding sites.

## Data Sources (from FlyBase GFF r6.66)

### 1. Promoters - Extended TSS Regions
- **Source**: TSS features (RAMPAGE + modENCODE)
- **Count**: 43,461 TSS (34,763 with gene associations)
- **Region definition**: TSS-2000bp to TSS+500bp (2.5kb extended promoter)
- **Experimental basis**: RAMPAGE-seq across 36 developmental stages including 10 embryonic timepoints (0-10h hourly)

### 2. Enhancers - CRMs from REDfly
- **Source**: regulatory_region features
- **Count**: 58,481 total (12,520 STARR-seq validated)
- **Gene assignment**: Nearest gene with distance annotation
- **Size**: median 499bp, mean 960bp

### 3. Silencers - Repressed Regions
- **Source**: regulatory_region features with "Repressed" in name
- **Count**: 331 regions (from STARR-seq negative controls in Kc167 cells)
- **Caveat**: Limited coverage. H3K27me3 ChIP-seq data would provide better genome-wide silencer/repressed chromatin annotation but requires additional data download.

### 4. TF Binding Sites - ChIP-seq/ChIP-chip
- **Source**: TF_binding_site features
- **Count**: 240,700 total binding sites
- **Key datasets**:
  - BDTNP early embryo: dl (25k), twi (11k), sna (4k), hb (3k)
  - Timed embryo ChIP: E2-4h, E4-6h, E6-8h, E8-10h, E10-12h
- **Analysis**: Genome-wide first, gene set enrichment later

## Implementation

### Phase 1: Sequence Extraction

#### `scripts/extract_promoters.py`
```
Input:  GFF (TSS features), genome FASTA
Output: promoters_sense.fasta, promoter_metadata.tsv
Region: TSS-2000 to TSS+500 (strand-aware)
Fields: region_id, fbgn, symbol, chrom, start, end, strand, tss_position
```

#### `scripts/extract_enhancers.py`
```
Input:  GFF (regulatory_region features), genome FASTA, gene annotations
Output: enhancers_sense.fasta, enhancer_metadata.tsv
Fields: region_id, chrom, start, end, source, nearest_gene, distance_to_gene,
        associated_gene (if annotated), is_starr_seq
```

#### `scripts/extract_silencers.py`
```
Input:  GFF (regulatory_region with "Repressed" name)
Output: silencers_sense.fasta, silencer_metadata.tsv
Fields: region_id, chrom, start, end, nearest_gene, distance_to_gene
Note:   Flag that H3K27me3 data would improve coverage
```

#### `scripts/extract_tfbs.py`
```
Input:  GFF (TF_binding_site features)
Output: tfbs_sense.fasta, tfbs_metadata.tsv
Fields: region_id, chrom, start, end, tf_name, tf_fbgn, library,
        nearest_gene, distance_to_gene
```

### Phase 2: BLAST Search
Use existing `blast_runner.py` with optimized parameters:
- word_size=7, gapopen=2, gapextend=1, dust=yes

### Phase 3: Analysis

#### `scripts/analyze_regulatory_te.py`
Unified analysis script handling all four region types:

1. **Per-region metrics**:
   - TE hit count, coverage (bp), density
   - TE family breakdown
   - Strand bias

2. **RepeatMasker overlap**:
   - Flag hits overlapping known RepeatMasker annotations
   - Calculate % novel vs known

3. **Gene-level aggregation**:
   - Sum hits across all regions per gene
   - Compare promoter vs enhancer vs silencer enrichment

4. **Output files**:
   - `{region_type}_te_hits.tsv` - Raw BLAST hits with RepeatMasker flag
   - `{region_type}_te_summary.tsv` - Per-region summary
   - `gene_{region_type}_te_summary.tsv` - Per-gene aggregation
   - `{region_type}_repeatmasker_overlap.tsv` - Known vs novel breakdown

## Output Directory Structure
```
results/regulatory_analysis/
├── promoter/
│   ├── promoter_te_hits.tsv
│   ├── promoter_te_summary.tsv
│   ├── gene_promoter_te_summary.tsv
│   └── promoter_repeatmasker_overlap.tsv
├── enhancer/
│   └── (same structure)
├── silencer/
│   └── (same structure)
├── tfbs/
│   └── (same structure)
└── integrated/
    ├── all_regions_gene_summary.tsv
    ├── region_comparison_stats.tsv
    └── analysis_summary.md
```

## Deduplication Strategy
- Promoters: By genomic coordinates (many TSS per gene)
- Enhancers: By genomic coordinates (overlapping CRMs common)
- TFBS: By genomic coordinates (same region bound by multiple TFs)
- Key: (chrom, start, end, te_id, te_start, te_end)

## Key Caveats to Document
1. **Silencer coverage**: Only 331 regions from STARR-seq repressed controls. H3K27me3 ChIP would provide better genome-wide coverage of silenced/repressive chromatin.
2. **Enhancer gene assignment**: Nearest-gene assignment is a simplification; enhancers can regulate distant genes.
3. **TFBS overlap**: Many TFBS overlap each other and with enhancers; coordinate-based deduplication essential.

## Execution Order
1. Extract all four region types (can run in parallel)
2. Run BLAST on each (can run in parallel)
3. Run unified analysis
4. Generate integrated summary
