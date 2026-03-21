# Comparative Analysis Scripts

Standalone one-off scripts for dmel x dpse comparative TE fossil analysis.
**Not part of the reusable fossil_finder pipeline** — project-specific exploration.

## Scripts

1. `conserved_te_hits.py` — Find TE hits shared between orthologous dmel/dpse UTRs
2. `dpse_only_hits.py` — Find TE hits in dpse UTRs absent from orthologous dmel UTRs
3. `alignment_degradation.py` — Compare alignment quality for shared TE hits across species

## Data Dependencies

- BLAST results in `results/comparative/`
- Ortholog table at `data/references/orthologs/dmel_dpse_orthologs.tsv`
- Region metadata in `results/comparative/{genome}/regions.tsv`
- GFF annotations for gene-to-transcript mapping
