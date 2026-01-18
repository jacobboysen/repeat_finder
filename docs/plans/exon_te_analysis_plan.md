# Exon TE Analysis Implementation Plan

## Overview
Extend the TE fossil mining pipeline to analyze individual exons for TE similarity, with annotations for UTR overlap, TE ORF type, and intrinsic disorder.

## Approach Decision
- **Exon source**: GFF parsing (not CDS FASTA) - enables individual exon tracking and UTR overlap annotation
- **Antisense control**: NO - not useful for coding sequences due to codon constraint; TE orientation detectable from sstart/send

## Scope
- **Input**: Genome-wide exon sequences extracted from GFF + genome
- **Output**: BLAST results with annotations for UTR position, TE domain type, and disorder scores
- **Gene sets**: Genome-wide first, then compare against existing gene sets

---

## Phase 1: Exon Extraction Pipeline

### 1.1 Create exon extraction script
**File**: `scripts/extract_exons.py`

**Approach**: Parse GFF to get exon coordinates, then extract sequences from genome FASTA (unlike UTR scripts which use pre-extracted FlyBase FASTAs).

**Tasks**:
1. Parse `data/references/dmel-all-r6.66.gff`:
   - Collect exon records (type=exon, Parent=FBtr)
   - Build transcript→gene mapping from mRNA records
   - Collect UTR coordinates for overlap detection
2. Load genome sequence from `data/references/dmel_genome.fasta` using BioPython
3. Extract exon sequences with headers containing:
   - Transcript ID (FBtr), parent gene (FBgn)
   - Exon number within transcript (1-indexed, sorted by genomic position)
   - Genomic coordinates (chr:start..end:strand)
   - Position flag: `first_exon`, `last_exon`, `internal_exon`
   - UTR overlap annotation
4. Generate sense and antisense versions (following existing pattern)
5. Deduplicate identical sequences (using SequenceCache pattern from extract_germ_plasm_3utrs.py)

**Header format**:
```
>FBtr0346766_exon1 gene=DIP-lambda FBgn=FBgn0267428 position=first_exon utr_overlap=5utr loc=2R:326768..326939:+ length=172
```

**Output**: `data/queries/genome_wide/exons_sense.fasta`

### 1.2 UTR overlap detection logic
**Within extraction script**:
- Parse 5'UTR and 3'UTR coordinates from GFF for each transcript
- For first exon on + strand (or last on - strand): check 5'UTR overlap
- For last exon on + strand (or first on - strand): check 3'UTR overlap
- Annotation values: `utr_overlap=5utr|3utr|both|none`

**Edge case handling**:
- Single-exon transcripts: annotate as `position=single_exon`
- UTR-only exons (no CDS overlap): include but flag as `cds_overlap=false`

---

## Phase 2: BLAST Execution

### 2.1 Run BLAST with established parameters
**Command**:
```bash
blastn -query data/queries/genome_wide/exons_sense.fasta \
  -db data/blastdb/dmel_te_flybase \
  -word_size 7 -gapopen 2 -gapextend 1 -dust yes -evalue 10 \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qseq sseq'
```

**Output**: `results/exon_analysis/genome_wide_exons.tsv`

**Note**: No antisense control needed - TE orientation (sense/antisense relative to gene) is determined from sstart vs send in BLAST output.

---

## Phase 3: TE Domain Annotation

### 3.1 Create TE domain classifier
**File**: `scripts/utils/te_domain_classifier.py`

**Approach**: Position-based heuristics using relative position within TE consensus

**Domain map** (for LTR retrotransposons):
```python
LTR_DOMAINS = {
    '5_ltr': (0.0, 0.10),    # First 10% of element
    'gag': (0.10, 0.35),     # ~10-35%
    'pol': (0.35, 0.90),     # ~35-90% (includes RT, RNaseH, integrase)
    '3_ltr': (0.90, 1.0),    # Last 10%
}

LINE_DOMAINS = {
    'orf1': (0.0, 0.30),     # First 30%
    'orf2_rt': (0.30, 0.70), # RT domain
    'orf2_en': (0.70, 1.0),  # Endonuclease domain
}

DNA_DOMAINS = {
    'transposase': (0.15, 0.85),  # Central coding region
    'tir': [(0.0, 0.15), (0.85, 1.0)],  # Terminal inverted repeats
}
```

**Function**:
```python
def classify_te_domain(te_id: str, sstart: int, send: int, te_length: int) -> str:
    """Return domain type: 'gag', 'pol', 'ltr', 'transposase', etc."""
```

### 3.2 Integrate with existing TE classification
- Use `parse_te_class()` from `data_loaders.py` to get superfamily (LTR/LINE/DNA)
- Apply appropriate domain map based on superfamily
- Add columns to output: `te_domain`, `te_class`

---

## Phase 4: Disorder Prediction Integration

### 4.1 Download protein sequences
**File**: `scripts/download_protein_sequences.py`

- Download FlyBase canonical proteins (or use existing if available)
- Map FBgn → protein sequence
- Output: `data/references/dmel_proteins.fasta`

### 4.2 Predict disorder via IUPred3 REST API
**File**: `scripts/predict_disorder.py`

**Tasks**:
1. Batch upload proteins to IUPred3 API (https://iupred3.elte.hu/)
2. Parse per-residue disorder scores
3. Calculate summary metrics per gene:
   - `disorder_fraction`: % residues with score > 0.5
   - `max_disordered_region`: longest contiguous disordered stretch
4. Output: `data/annotations/gene_disorder_scores.tsv`

### 4.3 Create disorder loader utility
**File**: `scripts/utils/disorder_loaders.py`

```python
def load_disorder_scores(path) -> Dict[str, float]:
    """Return {fbgn: disorder_fraction}"""

def get_disorder_category(fraction: float) -> str:
    """Return 'low' (<0.2), 'medium' (0.2-0.5), 'high' (>0.5)"""
```

---

## Phase 5: Analysis Script

### 5.1 Create exon analysis script
**File**: `scripts/analyze_exon_te.py`

**Outputs**:
1. `results/exon_analysis/exon_te_summary.tsv` - Per-exon statistics
2. `results/exon_analysis/gene_exon_te_summary.tsv` - Per-gene aggregation
3. `results/exon_analysis/te_domain_distribution.tsv` - Hits by TE domain type
4. `results/exon_analysis/utr_overlap_comparison.tsv` - First/last vs internal exons

**Key analyses**:
- Compare TE density: internal exons vs UTR-overlapping exons
- TE domain distribution: are exon hits enriched in gag/pol/env regions?
- Disorder correlation: do high-disorder genes have more exon TE hits?
- Strand bias: sense vs antisense in exons (compare to UTR patterns)

### 5.2 Generate comparison with UTR results
- Load existing 3'UTR and 5'UTR results
- Compare TE density across region types
- Identify genes with hits in multiple region types

---

## File Structure

```
scripts/
├── extract_exons.py                    # NEW: Exon extraction from GFF
├── predict_disorder.py                 # NEW: IUPred3 integration
├── analyze_exon_te.py                  # NEW: Main analysis script
├── download_protein_sequences.py       # NEW: Get FlyBase proteins
└── utils/
    ├── te_domain_classifier.py         # NEW: TE domain annotation
    └── disorder_loaders.py             # NEW: Disorder data loading

data/
├── queries/genome_wide/
│   └── exons_sense.fasta               # NEW: Extracted exons with position annotations
├── references/
│   └── dmel_proteins.fasta             # NEW: Protein sequences (for disorder)
└── annotations/
    └── gene_disorder_scores.tsv        # NEW: Disorder predictions

results/exon_analysis/
├── genome_wide_exons.tsv               # BLAST output
├── exon_te_summary.tsv                 # Per-exon analysis
├── gene_exon_te_summary.tsv            # Per-gene aggregation
├── te_domain_distribution.tsv          # Domain breakdown
└── utr_overlap_comparison.tsv          # Position-based comparison (terminal vs internal)
```

---

## Implementation Order

1. **extract_exons.py** - Get sequences with position annotations
2. **Run BLAST** - Generate raw results
3. **te_domain_classifier.py** - Add TE domain annotation utility
4. **analyze_exon_te.py** (basic) - Initial analysis without disorder
5. **download_protein_sequences.py** - Get proteins
6. **predict_disorder.py** - Run IUPred3
7. **disorder_loaders.py** - Loader utility
8. **analyze_exon_te.py** (full) - Add disorder integration

---

## Verification

1. **Exon extraction**: Check header format, verify first/last exon detection against known genes
2. **BLAST**: Compare hit counts to UTR analyses (expect similar or lower density in coding regions)
3. **TE domains**: Spot-check domain assignments against known TE structures
4. **Disorder**: Validate predictions against known disordered Drosophila proteins (e.g., Bicoid has disordered regions)
5. **End-to-end**: Run full pipeline on germ_plasm gene subset first as sanity check

---

## Expected Findings (Hypotheses to Test)

1. **Lower TE density in exons vs UTRs** - Coding regions under stronger selection
2. **Internal exons have fewer hits than terminal exons** - UTR overlap effect
3. **TE pol domain enrichment** - Pol has conserved RT/integrase motifs similar to host proteins
4. **Disordered regions correlate with TE similarity** - Low-complexity/flexible regions may tolerate TE-like sequences
5. **Gene expression proteins** - RNA-binding proteins may share motifs with TE ORFs

---

## Estimated Scale

- **Input**: ~200,000+ exons genome-wide (based on GFF counts)
- **BLAST output**: Expect 2-5M hits (similar to UTR analyses)
- **Processing time**: ~30 min for exon extraction, ~1-2 hours for BLAST (sense only), ~30 min for analysis
- **Disorder prediction**: ~2-4 hours for IUPred3 API (rate-limited)

---

## Summary

This plan adds exon analysis to complement the existing 5'UTR and 3'UTR pipelines:

1. **New capability**: Individual exon-level TE detection with position tracking
2. **Key annotations**: UTR overlap, TE domain type (gag/pol/env), disorder scores
3. **Comparison enabled**: Internal vs terminal exons, exons vs UTRs, disorder correlation
4. **Follows existing patterns**: Same BLAST parameters, utility module structure, output format
