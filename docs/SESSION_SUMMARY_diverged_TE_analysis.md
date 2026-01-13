# Session Summary: Diverged TE Detection with DUST Filtering

## Problem Statement

The original BLAST parameter sweep was **coalescing on simple sequence repeats** (AT-repeats, poly-A) rather than finding real diverged transposable element (TE) sequences. The user wanted to:

1. Detect ancient TE "fossils" - sequences that have accumulated mutations and indels over evolutionary time
2. Filter out simple/low-complexity repeats
3. Validate findings against proper controls

## Key Discovery

**Root cause**: With `dust=no`, BLAST was matching low-complexity regions that happen to be similar between 3'UTRs and TE sequences. These are not real TE fossils.

**Solution**: Enable DUST filtering (`dust=yes`) to mask simple repeats before alignment.

## Analysis Performed

### 1. Diagnosed the Simple Repeat Problem

Analyzed the original parameter sweep results:
- 90.5% of hits were simple repeats (AT-repeat, poly-A/T)
- Only 9.5% were "complex" sequences
- Almost no hits had gaps (0.3%) - indicating exact matches to simple repeats

Created `scripts/extract_top_hits.py` to visualize actual alignments and confirm the problem.

### 2. Re-ran Parameter Sweep with DUST Filtering

Tested 32 parameter combinations with `dust=yes`:

**Best parameters identified:**
```
word_size=7, gapopen=2, gapextend=1, penalty=-1, reward=1, dust=yes
```

**Results improvement:**
| Metric | Before (dust=no) | After (dust=yes) |
|--------|------------------|------------------|
| Complex sequences | 9.5% | 99.9% |
| Hits with gaps | 0.3% | 60.4% |
| Complex + gapped | 0 | 352 |

### 3. Ran Optimized BLAST on All Control Groups

Created control gene lists and ran BLAST with optimized parameters:

- **germ_plasm_sense/antisense**: 12 canonical germ plasm genes
- **housekeeping_sense**: 8 housekeeping genes (negative control)
- **somatic_sense**: 8 somatic-localized genes
- **cleared_sense**: 10 posteriorly-cleared genes
- **adult_sense**: 9 adult-expressed genes
- **shuffled**: Dinucleotide-shuffled germ plasm sequences (statistical control)

### 4. Validated Against Controls

**Key validation**: Long alignments (≥150bp, ≥200bp) are enriched in real sequences vs shuffled:

| Dataset | ≥150bp Enrichment | ≥200bp Enrichment |
|---------|-------------------|-------------------|
| germ_plasm | 11.3x | 6.3x |
| somatic | 98.3x | 47.7x |
| cleared | 46.5x | 16.3x |
| housekeeping | 0.7x | 0.7x |

**Housekeeping genes show NO enrichment** - confirms they are proper negative controls without TE-derived sequences.

## Key Findings

### Top Diverged TE Fossil Candidates

**Germ plasm genes:**
1. **tudor (tud)**: 256bp Stalker2 match (58.6% identity, 11 gaps)
2. **nanos (nos)**: 223bp Stalker2 match (59.2% identity, 9 gaps)
3. **piwi**: 200-208bp 17.6 element matches (60-62% identity)

**Somatic genes (strongest signal):**
1. **Syt1**: 510bp rover element match (58.6% identity, 29 gaps)

**Cleared genes:**
1. **Krüppel (Kr)**: 290-295bp Stalker2/17.6 matches (59-62% identity)

### TE Family Pattern

All top candidates are **LTR retrotransposons** (Stalker2, 17.6, rover), consistent with known Drosophila TE biology.

### Divergence Characteristics

- 58-62% identity (38-42% diverged from consensus)
- Multiple gaps (8-29 per alignment) indicate accumulated indels
- Long alignments (200-510bp) rule out chance similarity

These patterns are consistent with ancient TE insertions that occurred millions of years ago.

## Files Created/Modified

### New Scripts
- `scripts/extract_top_hits.py` - Extract and display top BLAST alignments
- `scripts/build_control_genelists.py` - Build control gene lists from FlyBase groups

### New Data Files
- `data/gene_lists/somatic_genes_consolidated.tsv`
- `data/gene_lists/cleared_genes_consolidated.tsv`
- `data/gene_lists/adult_genes_consolidated.tsv`
- `data/gene_lists/*_fbgn_ids.txt`

### Results (gitignored, local only)
- `results/dust_sweep/` - Parameter sweep with DUST filtering
- `results/diverged_controls/` - BLAST results for all control groups
- `results/DIVERGED_TE_ANALYSIS_SUMMARY.md`
- `results/DIVERGED_TE_CONTROL_COMPARISON.md`

## Next Steps (if desired)

1. **Conservation analysis**: Check if diverged TE regions are conserved across Drosophila species
2. **TE region mapping**: Determine which parts of TEs are being matched (LTR vs coding)
3. **Functional testing**: Investigate if these regions affect mRNA localization/translation
4. **Publication-ready figures**: Generate signal density plots for top candidates
