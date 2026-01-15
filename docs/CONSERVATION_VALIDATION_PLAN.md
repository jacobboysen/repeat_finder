# Conservation & Synteny Validation Plan for TE Hits

**Goal**: Validate BLAST TE hits in 3'UTRs using cross-species conservation to distinguish real ancient TE fossils from noise.

**Rationale**: If a TE-derived sequence is conserved across Drosophila species that diverged 5-40 million years ago, it's likely:
1. A real insertion (not alignment artifact)
2. Under functional constraint (regulatory role)
3. Present before species divergence (true "fossil")

---

## Available Resources

### 1. UCSC Conservation Scores (easiest)

| Resource | URL | Description |
|----------|-----|-------------|
| phyloP 27-way | `hgdownload.soe.ucsc.edu/goldenPath/dm6/phyloP27way/` | Per-base conservation scores (bigWig, 396MB) |
| phastCons 27-way | `hgdownload.soe.ucsc.edu/goldenPath/dm6/phastCons27way/` | HMM-based conserved elements |
| MAF alignments | `hgdownload.soe.ucsc.edu/goldenPath/dm6/multiz27way/maf/` | Actual sequence alignments |

**Species included** (27-way): 23 Drosophila + house fly, mosquito, honeybee, beetle

### 2. FlyBase Ortholog Data

- Orthology calls for D. simulans, D. yakuba, D. erecta, etc.
- Can extract orthologous 3'UTR sequences
- Download from FlyBase FTP or use DIOPT

### 3. Raw Genome Sequences

- D. simulans, D. yakuba genomes available from NCBI/FlyBase
- Can build custom BLAST databases

---

## Validation Approaches

### Approach 1: Conservation Score Overlay (Recommended First)

**Idea**: Compare phyloP conservation scores between TE-hit regions and non-hit regions within the same UTRs.

**Steps**:
1. Download `dm6.phyloP27way.bw` (396MB)
2. For each UTR with TE hits:
   - Extract phyloP scores for the entire UTR
   - Calculate mean score for TE-hit regions
   - Calculate mean score for non-hit regions
   - Compare distributions
3. If TE-hit regions are MORE conserved → supports functional TE fossil
4. If TE-hit regions are LESS conserved → suggests noise or neutral

**Tools needed**:
```bash
# Download bigWig
wget https://hgdownload.soe.ucsc.edu/goldenPath/dm6/phyloP27way/dm6.phyloP27way.bw

# Extract scores for regions (UCSC tools)
bigWigAverageOverBed dm6.phyloP27way.bw regions.bed output.tab
```

**Output**: Conservation score per TE hit, can use as filter/ranking

---

### Approach 2: Synteny Check via MAF Alignments

**Idea**: Look at actual sequence alignments to see if TE sequence is present in orthologous positions.

**Steps**:
1. Download MAF files for relevant chromosomes
2. For each high-confidence TE hit:
   - Extract the MAF block covering that UTR region
   - Check if D. simulans/yakuba have alignable sequence
   - If yes, check if that sequence also matches the TE
3. Classify hits as:
   - **Syntenic**: TE sequence present in multiple species → ancient insertion
   - **Mel-specific**: Only in D. melanogaster → recent insertion or artifact
   - **Diverged**: Aligned but sequence differs → under selection

**Tools needed**:
```bash
# Download MAF files
rsync -avz rsync://hgdownload.cse.ucsc.edu/goldenPath/dm6/multiz27way/maf/ ./maf/

# Parse with mafTools or BioPython
```

**Advantage**: Can see actual sequence, not just scores
**Disadvantage**: Larger data, more complex parsing

---

### Approach 3: Ortholog UTR BLAST

**Idea**: Extract orthologous 3'UTRs from related species and BLAST against same TE database.

**Steps**:
1. Get ortholog mapping: D. mel gene → D. sim/yak gene
2. Extract 3'UTR sequences from D. simulans/yakuba
3. BLAST orthologous UTRs against same TE database
4. Compare: Does the ortholog also hit the same TE family?

**Interpretation**:
- Same TE family hit in ortholog → ancient conserved insertion
- Different TE family → independent insertions (less interesting)
- No hit in ortholog → mel-specific or noise

**Data needed**:
- FlyBase ortholog table
- D. simulans/yakuba 3'UTR sequences (or extract from genome + GFF)

---

### Approach 4: TE Presence in Other Species Genomes

**Idea**: Check if the TE family itself exists in other species.

**Steps**:
1. For each TE family with UTR hits (roo, INE-1, 1360, etc.)
2. Check if this TE family is present in D. simulans/yakuba genomes
3. If TE family is mel-specific → hits could be recent insertions
4. If TE family is ancestral → supports ancient fossil hypothesis

**Note**: Most major TE families (roo, gypsy, copia) are shared across Drosophila

---

## Recommended Implementation Order

### Phase 1: Quick Conservation Score Analysis (1-2 hours)

```python
# Pseudocode
1. Download dm6.phyloP27way.bw
2. Convert UTR TE hit coordinates to BED format
3. Use bigWigAverageOverBed to get mean conservation per hit
4. Compare to background (shuffled regions or non-hit UTR regions)
5. Plot: conservation score vs hit quality (length, identity)
```

**Expected outcome**:
- If high-quality TE hits show elevated conservation → strong validation
- If no difference → need deeper investigation

### Phase 2: MAF-Based Synteny (2-4 hours)

```python
# Pseudocode
1. Download MAF files for chromosomes with top hits
2. For top 100 hits (≥80% id, ≥50bp):
   - Extract MAF block
   - Check D. simulans alignment
   - Check D. yakuba alignment
3. Classify each hit: syntenic vs mel-specific
4. Report: "X% of high-quality hits are syntenically conserved"
```

### Phase 3: Ortholog UTR Comparison (4-8 hours)

```python
# Pseudocode
1. Download FlyBase ortholog table
2. Download D. simulans 3'UTR sequences
3. For genes with top TE hits:
   - Find D. simulans ortholog
   - Extract ortholog's 3'UTR
   - BLAST against TE database
4. Compare hit patterns between orthologs
```

---

## Expected Results & Interpretation

| Scenario | Conservation | Synteny | Interpretation |
|----------|--------------|---------|----------------|
| High conservation, syntenic | phyloP > 0.5 | Present in sim/yak | **Strong evidence**: Ancient functional TE fossil |
| High conservation, not syntenic | phyloP > 0.5 | Mel-specific | Possible recent but constrained insertion |
| Low conservation, syntenic | phyloP < 0 | Present in sim/yak | Ancient but neutrally evolving |
| Low conservation, not syntenic | phyloP < 0 | Mel-specific | **Likely noise or recent neutral insertion** |

---

## Quick Start Commands

```bash
# 1. Download conservation data
cd data/references
wget https://hgdownload.soe.ucsc.edu/goldenPath/dm6/phyloP27way/dm6.phyloP27way.bw

# 2. Get UCSC tools
wget http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/bigWigAverageOverBed
chmod +x bigWigAverageOverBed

# 3. Create BED file of TE hit regions (need to convert UTR coords to genomic)
# (script needed)

# 4. Extract conservation scores
./bigWigAverageOverBed dm6.phyloP27way.bw te_hits.bed te_hits_conservation.tab
```

---

## References

- [UCSC phyloP 27-way](https://hgdownload.soe.ucsc.edu/goldenPath/dm6/phyloP27way/)
- [UCSC multiz 27-way MAF](https://hgdownload.soe.ucsc.edu/goldenPath/dm6/multiz27way/)
- [FlyBase Orthologs](https://flybase.org/commentaries/2012_04/improved_orthologies.html)
- [Evolution of genes and genomes on the Drosophila phylogeny](https://www.nature.com/articles/nature06341)
