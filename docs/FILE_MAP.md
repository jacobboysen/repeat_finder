# Repository File Map

> **Last updated**: 2026-01-14
> **Purpose**: Track what each file contains and how the repo has evolved
> **Maintenance**: Update this file when adding new scripts, data, or results

## Recent Changes (2026-01-14)

### Ancient TE Candidates: Annotate, Don't Filter

**Current recommended output**: `results/utr_te_loci/ancient_te_loci_by_isoform.tsv`

Earlier versions (v1, v2) aggressively filtered candidates based on CDS overlap and multi-family
hits. This was **too aggressive** and made unwarranted biological assumptions:

1. **CDS overlap ≠ invalid UTR hit** - Drosophila has a compact genome with frequent overlapping
   features. A region can be CDS in one isoform and functional UTR in another (dual-function).
   CDS overlap is now **annotated** (`cds_overlap: yes/no`) rather than filtered.

2. **Multi-family ≠ spurious** - Multiple TE families hitting the same region could indicate:
   - Related families sharing sequences (PBS, regulatory motifs)
   - Nested insertions (one TE into another)
   - The original concern was low-complexity sequence, not multi-family per se
   Multi-family is now **annotated** (`n_te_families`, `te_families`) rather than filtered.

**Data files**:
- `results/utr_te_loci/ancient_te_loci_by_isoform.tsv` - **CURRENT** (25,240 loci, annotated not filtered)
- `results/repeatmasker_analysis_v2/` - v2 filtered data (use with caution)
- `results/archive/ancient_te_candidates_v1_FLAWED/` - v1 archived data

**Naming convention**: `{transcript}|{chrom}|{start}-{end}` (one UTR model per isoform)

See `results/utr_te_loci/METHODOLOGY_NOTES.md` for detailed rationale.

---

- **Added** External annotation pipeline for TE enrichment analysis
  - Downloads FlyBase expression, GO terms, and FlyFISH localization
  - Creates functional gene sets based on expression, GO, and localization
  - Performs Fisher's exact test and Mann-Whitney U enrichment analysis
  - See `scripts/download_external_annotations.py` and related scripts
- **Added** `scripts/utils/annotation_loaders.py` - Parsers for external annotations
- **Added** `data/annotations/` - Downloaded and processed annotation files
- **Added** `data/gene_lists/functional/` - Expression/GO/FlyFISH gene sets
- **Added** `figures/enrichment/` - Enrichment analysis visualizations

## Previous Changes (2026-01-14)

- **Added** Phase 2 synteny analysis using MAF alignments
  - Key finding: **98.5% of all hits syntenic, but only 57% of HQ hits**
  - Higher identity = LESS syntenic (86% at 80-85% → 10% at 100%)
  - See `results/repeatmasker_analysis/SYNTENY_ANALYSIS_RESULTS.md`
- **Added** `scripts/conservation_analysis/` - phyloP conservation validation
  - Key finding: **Stricter quality filters yield LESS conserved hits**
  - All hits: 61% conserved → HQ (≥80%/≥50bp): 25% conserved
  - See `results/repeatmasker_analysis/CONSERVATION_ANALYSIS_RESULTS.md`
- **Added** `data/references/dm6.phyloP27way.bw` - UCSC 27-way conservation scores
- **Added** `data/references/maf/` - UCSC 27-way MAF alignments (~1GB)
- **Added** `results/repeatmasker_analysis/` - RepeatMasker vs BLAST comparison
  - Key finding: **83.2% of BLAST hits are NOT in RepeatMasker**
  - Even 100% identity matches (up to 897bp) missed by RepeatMasker
  - See `REPEATMASKER_COMPARISON_SUMMARY.md` for full analysis
- **Added** `scripts/analyze_repeatmasker_overlap.py` - RepeatMasker comparison script
- **Added** `data/references/dm6.fa.out` - UCSC RepeatMasker annotations
- **Added** `scripts/utils/` - Shared utility modules (paths.py, blast_io.py, data_loaders.py)
- **Renamed** compare scripts for clarity:
  - `compare_gene_sets.py` → `rank_gene_sets.py`
  - `compare_genes.py` → `visualize_gene_comparison.py`
  - `compare_5utr_3utr.py` → `analyze_utr_variants.py`
  - `compare_all_gene_sets.py` → `summarize_all_gene_sets.py`
- **Archived** old results to `results/archive/` (BAD_*, INTERMEDIATE_*)
- **Removed** dead code in `src/bioinfo/` (was template code)

---

## Quick Reference

| What you want | Where to find it |
|---------------|------------------|
| Latest germ plasm TE hits | `results/diverged_controls/germ_plasm_sense.tsv` |
| **Synteny analysis** | `results/repeatmasker_analysis/SYNTENY_ANALYSIS_RESULTS.md` |
| **Conservation analysis** | `results/repeatmasker_analysis/CONSERVATION_ANALYSIS_RESULTS.md` |
| **RepeatMasker comparison** | `results/repeatmasker_analysis/REPEATMASKER_COMPARISON_SUMMARY.md` |
| **UTR visualizations (HTML)** | `results/te_annotations/index.html` |
| **TE visualizations (HTML)** | `results/te_annotations/te_index.html` |
| Strand analysis summary | `results/STRAND_ANALYSIS_SUMMARY.md` |
| Analysis summary | `docs/SESSION_SUMMARY_diverged_TE_analysis.md` |
| BLAST parameters | `results/DIVERGED_TE_ANALYSIS_SUMMARY.md` |
| Gene lists | `data/gene_lists/` |
| 3'UTR sequences | `data/queries/{group}/3UTR_sense.fasta` |
| Run BLAST | `scripts/blast_runner.py` |

---

## Directory Structure

```
repeat_finder/
├── data/                    # Input data (tracked, some downloaded)
│   ├── gene_lists/          # Gene sets for different tissue groups
│   ├── queries/             # Extracted 3'UTR sequences for BLAST
│   ├── references/          # Reference genomes, TEs, annotations
│   ├── blastdb/             # BLAST databases (gitignored)
│   └── cache/               # API response caches
├── scripts/                 # Analysis scripts (tracked)
│   └── utils/               # Shared utility modules (paths, blast_io, data_loaders)
├── results/                 # BLAST outputs, analyses (gitignored)
│   └── archive/             # Superseded results (pre_dust/, intermediate/)
├── reports/                 # HTML/MD reports (gitignored)
├── figures/                 # Generated plots (gitignored)
├── docs/                    # Documentation (tracked)
└── workflows/               # Snakemake workflows
```

---

## Scripts

### Shared Utilities (`scripts/utils/`)

| Module | Purpose | Key Functions |
|--------|---------|---------------|
| `paths.py` | Centralized path resolution | `get_project_root()`, `get_results_dir()`, `get_gene_lists_dir()` |
| `blast_io.py` | BLAST result parsing | `BLAST_COLUMNS`, `load_blast_results()`, `classify_strand()` |
| `data_loaders.py` | Data loading utilities | `load_gene_list()`, `parse_fasta()`, `load_te_database()` |
| `annotation_loaders.py` | External annotation parsing | `load_rnaseq_expression()`, `load_go_annotations()`, `load_flyfish_localization()` |

### Core Pipeline Scripts

| Script | Purpose | Input | Output |
|--------|---------|-------|--------|
| `download_references.py` | Download FlyBase/Dfam reference data | URLs | `data/references/` |
| `build_databases.py` | Build BLAST databases | FASTAs | `data/blastdb/` |
| `blast_runner.py` | Run BLAST searches | Query FASTA, DB | TSV results |

### External Annotation Pipeline Scripts

| Script | Purpose | Input | Output |
|--------|---------|-------|--------|
| `download_external_annotations.py` | Download FlyBase expression, GO, FlyFISH | URLs | `data/annotations/raw/` |
| `build_annotation_table.py` | Build unified gene annotation table | Raw annotations | `data/annotations/gene_annotations.tsv` |
| `build_functional_gene_sets.py` | Create gene sets by expression/GO/localization | Annotation table | `data/gene_lists/functional/` |
| `analyze_functional_enrichment.py` | TE enrichment in functional gene sets | Gene sets + TE data | `results/functional_te_enrichment.tsv` |
| `visualize_functional_enrichment.py` | Generate enrichment visualizations | Enrichment results | `figures/enrichment/` |

**Usage:**
```bash
python scripts/download_external_annotations.py
python scripts/build_annotation_table.py
python scripts/build_functional_gene_sets.py
python scripts/analyze_functional_enrichment.py
python scripts/visualize_functional_enrichment.py
```

### Gene List Scripts

| Script | Purpose | Input | Output |
|--------|---------|-------|--------|
| `build_germ_plasm_genelist.py` | Build canonical germ plasm gene list | Hardcoded + FlyBase | `data/gene_lists/germ_plasm_*` |
| `build_housekeeping_genelist.py` | Build housekeeping control list | FlyBase API | `data/gene_lists/housekeeping_*` |
| `build_control_genelists.py` | Build somatic/cleared/adult lists | FlyBase API | `data/gene_lists/{group}_*` |
| `list_gene_groups.py` | Query available FlyBase gene groups | FlyBase API | stdout |

### 3'UTR Extraction Scripts

| Script | Purpose | Input | Output |
|--------|---------|-------|--------|
| `extract_germ_plasm_3utrs.py` | Extract 3'UTRs for gene lists | Gene list + reference FASTA | `data/queries/{group}/` |
| `shuffle_sequences.py` | Dinucleotide shuffle for controls | FASTA | `*_shuffled.fasta` |
| `filter_queries.py` | Filter queries by criteria | FASTA | Filtered FASTA |

### Analysis Scripts

| Script | Purpose | Input | Output |
|--------|---------|-------|--------|
| `analyze_genome_wide_te.py` | Genome-wide TE analysis | BLAST TSV | Top/bottom genes, strand bias |
| `analyze_full_transcripts.py` | Full transcript TE analysis | Reference FASTAs | HTML visualizations |
| `analyze_te_families.py` | TE family distribution | BLAST TSV | Family stats |
| `analyze_te_regions.py` | Map hits to TE regions | BLAST TSV + GFF | Region enrichment |
| `analyze_repeatmasker_overlap.py` | Compare BLAST hits vs RepeatMasker | BLAST TSV + dm6.fa.out | Known/novel classification |
| `rank_gene_sets.py` | Rank gene sets by TE metrics | Results TSVs | Console tables |
| `analyze_utr_variants.py` | Compare 5'UTR vs 3'UTR | Results TSVs | Console tables |
| `summarize_all_gene_sets.py` | Summary across all gene sets | Results TSVs | Console tables |
| `te_parameter_sweep.py` | Test BLAST parameters | Queries + DB | `results/parameter_sweep/` |
| `calculate_te_signal_density.py` | Position-wise TE signal | BLAST TSV | Density arrays |
| `detect_te_clusters.py` | Find TE signal hotspots | Density data | Cluster coords |

### Conservation Analysis Scripts (`scripts/conservation_analysis/`)

| Script | Purpose | Input | Output |
|--------|---------|-------|--------|
| `convert_all_hits_to_bed.py` | Convert TE hits to genomic BED format | BLAST TSVs + GFF | `te_hits_all_genomic.bed` |
| `analyze_all_hits_conservation.py` | Analyze phyloP conservation scores | Conservation tab | Summary stats |
| `plot_all_vs_hq_conservation.py` | Visualize conservation by quality | Conservation tab | PNG figures |
| `fast_synteny_analysis.py` | Analyze MAF synteny for all hits | MAF + BED | Synteny TSV |
| `fast_synteny_hq.py` | Analyze MAF synteny for HQ hits | MAF + BED | Summary stats |
| `summarize_synteny.py` | Summarize synteny results | Synteny TSV | Tables |
| `plot_synteny_results.py` | Visualize synteny findings | Synteny TSV | PNG figures |
| `README.md` | Workflow documentation | - | - |

**Key Findings**:
- Conservation: Stricter quality filters yield LESS conserved hits (61% → 8%)
- Synteny: 98.5% of all hits syntenic, but only 57% of HQ hits

### Visualization Scripts

| Script | Purpose | Output |
|--------|---------|--------|
| `plot_te_signal.py` | Per-gene signal density plots | `figures/te_signal/` |
| `plot_parameter_sweep.py` | Parameter sweep heatmaps | `figures/parameter_sweep/` |
| `visualize_gene_comparison.py` | Cross-gene comparison plots | `figures/` |
| `generate_te_fossil_report.py` | Comprehensive HTML/MD report | `reports/` |

### Utility Scripts

| Script | Purpose | Output |
|--------|---------|--------|
| `extract_top_hits.py` | Display top BLAST alignments | Formatted text |
| `summarize_results.py` | Quick BLAST result summary | Stats |
| `list_databases.py` | List available BLAST databases | stdout |

---

## Data Files

### External Annotations (`data/annotations/`)

| Path | Contents | Source |
|------|----------|--------|
| `raw/gene_rpkm_report.tsv` | RNA-Seq RPKM by tissue/stage | FlyBase |
| `raw/gene_rpkm_matrix.tsv` | RNA-Seq RPKM matrix | FlyBase |
| `raw/gene_association.gaf` | GO annotations (GAF format) | FlyBase |
| `raw/gene_group_data.tsv` | Gene group memberships | FlyBase |
| `raw/flyfish_localization.csv` | RNA localization patterns | FlyFISH |
| `gene_annotations.tsv` | Unified per-gene annotation table | Built |

### Functional Gene Sets (`data/gene_lists/functional/`)

| Pattern | Type | Description |
|---------|------|-------------|
| `expr_{tissue}_high_fbgn_ids.txt` | Expression | High RPKM (>10) in tissue |
| `expr_{tissue}_enriched_fbgn_ids.txt` | Expression | 2x enriched vs mean |
| `expr_ubiquitous_fbgn_ids.txt` | Expression | Expressed in 5+ tissues |
| `expr_maternal_fbgn_ids.txt` | Expression | High early embryo, low later |
| `expr_germline_specific_fbgn_ids.txt` | Expression | High ovary, low elsewhere |
| `go_{category}_fbgn_ids.txt` | GO | Genes with specific GO terms |
| `flyfish_{pattern}_fbgn_ids.txt` | FlyFISH | Genes with localization pattern |
| `gene_sets_summary.tsv` | Summary | All sets with counts |

### Gene Lists (`data/gene_lists/`)

| File Pattern | Contents | Genes |
|--------------|----------|-------|
| `germ_plasm_genes_consolidated.tsv` | Canonical germ plasm genes | nos, osk, pgc, gcl, vas, aub, tud, piwi, AGO3, dhd, CycB, Hsp83 |
| `germ_plasm_genes_tier1.txt` | Priority subset | nos, osk, pgc, gcl, vas, aub, CycB |
| `housekeeping_genes_consolidated.tsv` | Negative control genes | RpL32, Act5C, Gapdh1, etc. |
| `somatic_genes_consolidated.tsv` | Somatically-localized genes | bcd, nos, osk, stau, etc. |
| `cleared_genes_consolidated.tsv` | Posteriorly-cleared genes | hb, Kr, kni, etc. |
| `adult_genes_consolidated.tsv` | Adult-expressed genes | Various |
| `*_fbgn_ids.txt` | FlyBase gene IDs only | FBgn IDs |
| `*_status.json` | Build metadata | Timestamps, counts |

### Query Sequences (`data/queries/{group}/`)

Each group (germ_plasm, housekeeping, somatic, cleared, adult) contains:

| File | Contents |
|------|----------|
| `3UTR_sense.fasta` | All 3'UTRs, sense strand |
| `3UTR_antisense.fasta` | Reverse complement (control) |
| `3UTR_sense_tier1.fasta` | Priority genes only |
| `3UTR_antisense_tier1.fasta` | Priority genes, antisense |
| `3UTR_shuffled.fasta` | Dinucleotide-shuffled (germ_plasm only) |
| `isoform_map.json` | Gene → transcript mapping |

### Reference Data (`data/references/`)

| File | Source | Contents |
|------|--------|----------|
| `dmel_3utr.fasta` | FlyBase | 30,324 3'UTR sequences |
| `dmel_5utr.fasta` | FlyBase | 5'UTR sequences |
| `dmel_cds.fasta` | FlyBase | Coding sequences |
| `dmel_te_flybase.fasta` | FlyBase | 5,734 TE instances |
| `dmel_genome.fasta` | FlyBase | Genome assembly |
| `dmel_annotation.gff` | FlyBase | Gene annotations |
| `dmel-all-r6.66.gff` | FlyBase | Full GFF (extracted from tar) |
| `dm6.fa.out` | UCSC | RepeatMasker annotations (137,555 repeats) |
| `dm6.phyloP27way.bw` | UCSC | phyloP conservation scores, 27-way alignment (396 MB) |
| `tools/bigWigAverageOverBed` | UCSC | Tool for extracting conservation scores |
| `gene_info.tsv` | FlyBase | FBgn/FBtr/FBpp mapping |
| `gene_groups.tsv` | FlyBase | Gene group definitions |

---

## Results (gitignored, local only)

### Current Best Results

| File | What it contains | Parameters |
|------|------------------|------------|
| `results/diverged_controls/germ_plasm_sense.tsv` | **CURRENT** - Germ plasm TE hits with DUST | word_size=7, gapopen=2, gapextend=1, dust=yes |
| `results/diverged_controls/germ_plasm_antisense.tsv` | Antisense control | Same |
| `results/diverged_controls/housekeeping_sense.tsv` | Negative control | Same |
| `results/diverged_controls/somatic_sense.tsv` | Somatic comparison | Same |
| `results/diverged_controls/cleared_sense.tsv` | Cleared comparison | Same |
| `results/diverged_controls/adult_sense.tsv` | Adult comparison | Same |
| `results/diverged_controls/shuffled.tsv` | Statistical control | Same |

### Archived Results (`results/archive/`)

Old and intermediate results have been moved to `results/archive/` for cleanliness:

| Directory | Contents | Why archived |
|-----------|----------|--------------|
| `archive/pre_dust/` | BAD_* files | Pre-DUST filtering (90% simple repeats) |
| `archive/intermediate/` | INTERMEDIATE_* files | Development iterations |

These are preserved for reference but should not be used for analysis.

### Analysis Outputs

| File | Contents |
|------|----------|
| `results/DIVERGED_TE_ANALYSIS_SUMMARY.md` | Parameter comparison, top candidates |
| `results/DIVERGED_TE_CONTROL_COMPARISON.md` | Cross-group enrichment analysis |
| `results/STRAND_ANALYSIS_SUMMARY.md` | Sense vs antisense orientation analysis |
| `results/germ_plasm_strand_summary.tsv` | Per-gene strand statistics |
| `results/te_fossil_candidates.txt` | Curated candidate list |
| `results/te_fossils_diverged.txt` | Top 50 diverged hits with alignments |
| `results/top_hits_*.txt` | Per-gene hit summaries |

### RepeatMasker Comparison (`results/repeatmasker_analysis/`)

| File | Contents |
|------|----------|
| `REPEATMASKER_COMPARISON_SUMMARY.md` | Full analysis writeup with case study |
| `blast_hits_known_repeatmasker.tsv` | BLAST hits overlapping RepeatMasker (432,351) |
| `blast_hits_novel.tsv` | BLAST hits NOT in RepeatMasker (2,141,575) |
| `three_prime_UTR_repeatmasker_overlaps.tsv` | RepeatMasker annotations in UTR regions |

### UTR TE Loci (`results/utr_te_loci/`) - CURRENT

**Recommended output - annotated, not filtered.**

| File | Contents |
|------|----------|
| `ancient_te_loci_by_isoform.tsv` | **25,240 loci** across 14,423 UTR isoforms |
| `METHODOLOGY_NOTES.md` | Rationale for annotate-not-filter approach |

Key features:
- One UTR model per transcript isoform (preserves all isoforms)
- CDS overlap **annotated** (`cds_overlap`, `cds_overlap_transcripts`) not filtered
- Multi-family hits **annotated** (`n_te_families`, `te_families`) not filtered
- Naming: `{transcript}|{chrom}|{start}-{end}`

### Ancient TE Candidates v2 (`results/repeatmasker_analysis_v2/`) - SUPERSEDED

**Overly aggressive filtering - use `utr_te_loci/` instead.**

| File | Contents |
|------|----------|
| `ancient_te_candidates_clean.tsv` | 11,100 candidates (CDS-filtered, single-family only) |
| `ancient_te_candidates_suspicious.tsv` | 13,112 multi-family loci (wrongly excluded) |
| `FILTERING_REPORT_v2.md` | Filtering methodology |

This version filtered out CDS overlaps and multi-family hits, which was too aggressive.
See `utr_te_loci/METHODOLOGY_NOTES.md` for why this approach was revised.

### HTML Visualizations (`results/te_annotations/`)

Interactive HTML visualizations with colored sequences and nucleotide alignments.

**UTR-centric (where do TEs match on each UTR?):**

| File | Contents |
|------|----------|
| `index.html` | Index linking to all UTR visualizations |
| `{gene}_3UTR_TE_annotation.html` | Per-gene: sequence with TE regions colored by strand |

Genes covered: nos, osk, piwi, tud, vas, gcl, aub, pgc, CycB, Kr

**TE-centric (where do UTRs match on each TE?):**

| File | Contents |
|------|----------|
| `te_index.html` | Index + database coverage stats |
| `tes/{TE_ID}_{family}_UTR_matches.html` | Per-TE: sequence with UTR match regions highlighted |

Top TEs visualized: mdg1, roo, antonia, Stalker2, gypsy12, Max, F, HMS-Beagle

**Visualization features:**
- Color-coded by strand: blue=sense, yellow=antisense, green=both
- TE class labels (LTR, LINE, DNA, Helitron)
- Nucleotide-level alignments with match/mismatch coloring
- Per-region breakdown with supporting evidence

---

## Reports (`reports/`)

Multiple timestamped versions exist. Latest is most current:

| Pattern | Contents |
|---------|----------|
| `te_fossil_analysis_{timestamp}.html` | Full HTML report |
| `te_fossil_analysis_{timestamp}.md` | Markdown version |

---

## Evolution / Version History

### Version 1: Simple Repeat Problem (superseded)
- **Parameters**: dust=no
- **Problem**: 90% of hits were AT-repeats, poly-A
- **Files**: `results/parameter_sweep/`, `results/controls/`

### Version 2: DUST-Filtered (current)
- **Parameters**: dust=yes, word_size=7, gapopen=2, gapextend=1
- **Improvement**: 99.9% complex sequences, 60% with gaps
- **Files**: `results/diverged_controls/`

### Version 3: Strand-Aware Analysis
- Added sense vs antisense orientation tracking
- Created HTML visualizations for UTRs and TEs
- **Key finding**: Gene-specific strand biases (piwi 80% antisense, tud 67% sense)
- **Files**: `results/te_annotations/`, `results/STRAND_ANALYSIS_SUMMARY.md`

### Version 4: TE Structural Region Analysis
- Downloaded Bergman Lab TE consensus sequences with structural annotations
- Built separate BLAST database against consensus TEs
- **Key finding**: UTRs preferentially hit LTR regions (92.7%) over CDS (20.2%)
- **Key finding**: piRNA pathway genes (piwi, AGO3, aub) show strongest LTR enrichment
- **New resources**:
  - `data/references/dmel_te_consensus.fasta` - 127 consensus sequences
  - `data/references/te_annotations.gff` - LTR, CDS positions per family
  - `data/blastdb/dmel_te_consensus.*` - Consensus BLAST database
  - `scripts/analyze_te_regions.py` - Structural region analysis
  - `results/TE_REGION_ENRICHMENT_ANALYSIS.md` - Detailed findings
  - `results/{group}_te_regions.tsv` - Per-group statistics

### Version 5: RepeatMasker Comparison (current)
- Downloaded UCSC dm6 RepeatMasker annotations
- Compared BLAST TE hits against RepeatMasker to classify "known" vs "novel"
- **Key finding**: **83.2% of BLAST hits are NOT in RepeatMasker**
- **Key finding**: Even 100% identity matches (up to 897bp) missed by RepeatMasker
- **Why**: RepeatMasker only annotates ~27% of some TEs (ends only, misses internal regions)
- **New resources**:
  - `data/references/dm6.fa.out` - RepeatMasker annotations (137,555 repeats)
  - `scripts/analyze_repeatmasker_overlap.py` - Comparison script
  - `results/repeatmasker_analysis/` - Analysis outputs
  - `results/repeatmasker_analysis/REPEATMASKER_COMPARISON_SUMMARY.md` - Full writeup

---

## Key Statistics

### TE Database Coverage (germ plasm UTRs)
| Threshold | TEs matched | % of 5,734 total |
|-----------|-------------|------------------|
| ≥50bp | 653 | 11.4% |
| ≥100bp | 138 | 2.4% |
| ≥150bp | 55 | 1.0% |
| ≥200bp | 18 | 0.3% |

### Strand Orientation by Dataset
| Dataset | % Sense | % Antisense |
|---------|---------|-------------|
| germ_plasm | 39% | 61% |
| somatic | 69% | 31% |
| cleared | 55% | 45% |
| shuffled | 40% | 60% |

### TE Structural Region Enrichment
| Dataset | LTR % | CDS % | LTR/CDS Ratio |
|---------|-------|-------|---------------|
| germ_plasm | 92.7% | 20.2% | 4.6x |
| housekeeping | 81.5% | 16.3% | 5.0x |
| somatic | 102.1% | 25.3% | 4.0x |
| shuffled | 76.6% | 36.2% | 2.1x |

*LTRs are regulatory regions; CDS are coding regions*
*Shuffled controls show random baseline*

### Top TE Families in Germ Plasm UTRs
1. roo (LTR) - 190 hits
2. 1360 (DNA) - 187 hits
3. mdg1 (LTR) - 105 hits
4. S (DNA) - 83 hits
5. opus (LTR) - 73 hits

### RepeatMasker Comparison (genome-wide 3'UTRs)
| Category | Hits | Percentage |
|----------|------|------------|
| Known (in RepeatMasker) | 432,351 | 16.8% |
| Novel (NOT in RepeatMasker) | 2,141,575 | **83.2%** |

Novel hits by identity threshold:
| Threshold | Novel Hits | Transcripts |
|-----------|------------|-------------|
| ≥95% | 31,355 | 6,020 |
| ≥90% | 109,670 | 14,010 |
| ≥80% | 561,685 | 27,010 |
| ≥70% | 1,506,940 | 29,281 |

---

## BLAST Output Format

All TSV files use 17-column format:

```
qseqid, sseqid, pident, length, mismatch, gapopen,
qstart, qend, sstart, send, evalue, bitscore,
qlen, slen, qseq, sseq, strand
```

The `strand` column indicates TE orientation:
- `plus`: sstart < send (TE in sense orientation)
- `minus`: sstart > send (TE in antisense orientation)

---

## Maintenance Notes

When adding new files:
1. Update the relevant section in this document
2. If superseding old results, move old files to "Superseded" section
3. Document the parameters/methods used
4. Update the "Quick Reference" table if adding key outputs
