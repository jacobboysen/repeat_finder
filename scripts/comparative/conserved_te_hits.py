#!/usr/bin/env python3
"""Find TE hits conserved between orthologous dmel/dpse 3'UTRs.

A "conserved hit" means the same dmel TE instance is hit by both the dmel
UTR and its orthologous dpse UTR. Two sub-categories:

  - Same TE region: both species' UTRs hit the same part of the TE
    (overlapping sstart-send ranges), suggesting the TE insertion predates
    the species split (~25 Mya).

  - Different TE region: both UTRs hit the same TE instance but at different
    coordinates, suggesting independent or rearranged insertions from the
    same TE family.

Uses the dmel self-BLAST and dpse->dmel cross-BLAST (both query against the
dmel RM-derived TE database, so sseqid is directly comparable).

Output: results/comparative/analysis/conserved_te_hits.tsv
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from _shared import (
    RESULTS_DIR,
    build_ortholog_map,
    build_region_to_gene_dmel,
    build_region_to_gene_dpse,
    load_blast,
    load_orthologs,
    DATA_DIR,
)


def main():
    out_dir = RESULTS_DIR / "analysis"
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading BLAST results...")
    # dmel UTRs vs dmel TEs (self)
    dmel_self = load_blast(RESULTS_DIR / "dm6" / "blast_results.tsv")
    # dpse UTRs vs dmel TEs (cross)
    dpse_cross = load_blast(RESULTS_DIR / "UCI_Dpse_MV25_vs_dm6_te" / "blast_results.tsv")

    print(f"  dmel self: {len(dmel_self):,} hits")
    print(f"  dpse->dmel cross: {len(dpse_cross):,} hits")

    # Build region -> gene mappings
    print("Building gene mappings...")
    dmel_r2g = build_region_to_gene_dmel(
        RESULTS_DIR / "dm6" / "regions.tsv",
        DATA_DIR / "references" / "dmel_annotation_coding.gff3",
    )
    dpse_r2g = build_region_to_gene_dpse(
        RESULTS_DIR / "UCI_Dpse_MV25" / "regions.tsv",
        DATA_DIR / "references" / "dpse_annotation.gff3",
    )

    # Add gene columns
    dmel_self["gene_id"] = dmel_self["qseqid"].map(dmel_r2g)
    dpse_cross["gene_id"] = dpse_cross["qseqid"].map(dpse_r2g)

    # Build ortholog map: dmel FBgn -> dpse gene-LOC*
    orthologs = load_orthologs()
    dmel_to_dpse = build_ortholog_map(orthologs)
    dpse_to_dmel = {v: k for k, v in dmel_to_dpse.items()}

    # Map dpse genes back to dmel orthologs for joining
    dpse_cross["dmel_ortholog"] = dpse_cross["gene_id"].map(dpse_to_dmel)

    # Filter to hits where we have ortholog pairs
    dmel_with_ortho = dmel_self[dmel_self["gene_id"].isin(dmel_to_dpse)].copy()
    dpse_with_ortho = dpse_cross[dpse_cross["dmel_ortholog"].notna()].copy()

    print(f"  dmel hits with orthologs: {len(dmel_with_ortho):,}")
    print(f"  dpse hits with orthologs: {len(dpse_with_ortho):,}")

    # Find shared TE hits: same sseqid hit by orthologous gene pairs
    # Key: (dmel_gene, te_instance)
    print("Finding conserved TE hits...")

    dmel_te_hits = dmel_with_ortho.groupby(["gene_id", "sseqid"]).agg(
        dmel_sstart=("sstart", "min"),
        dmel_send=("send", "max"),
        dmel_pident=("pident", "max"),
        dmel_evalue=("evalue", "min"),
        dmel_qseq=("qseq", "first"),
        dmel_sseq=("sseq", "first"),
        dmel_qseqid=("qseqid", "first"),
        dmel_nhits=("sseqid", "size"),
    ).reset_index()

    dpse_te_hits = dpse_with_ortho.groupby(["dmel_ortholog", "sseqid"]).agg(
        dpse_sstart=("sstart", "min"),
        dpse_send=("send", "max"),
        dpse_pident=("pident", "max"),
        dpse_evalue=("evalue", "min"),
        dpse_qseq=("qseq", "first"),
        dpse_sseq=("sseq", "first"),
        dpse_qseqid=("qseqid", "first"),
        dpse_nhits=("sseqid", "size"),
    ).reset_index()

    # Merge on (dmel_gene, TE instance)
    conserved = dmel_te_hits.merge(
        dpse_te_hits,
        left_on=["gene_id", "sseqid"],
        right_on=["dmel_ortholog", "sseqid"],
        how="inner",
    )

    print(f"  Conserved (gene, TE) pairs: {len(conserved):,}")

    if conserved.empty:
        print("No conserved hits found.")
        return

    # Classify: same TE region vs different TE region
    # Overlapping if dmel [sstart, send] overlaps dpse [sstart, send]
    conserved["te_region_overlap"] = (
        (conserved["dmel_sstart"] <= conserved["dpse_send"]) &
        (conserved["dmel_send"] >= conserved["dpse_sstart"])
    )
    conserved["overlap_bp"] = (
        conserved[["dmel_send", "dpse_send"]].min(axis=1) -
        conserved[["dmel_sstart", "dpse_sstart"]].max(axis=1) + 1
    ).clip(lower=0)

    same_region = conserved["te_region_overlap"].sum()
    diff_region = (~conserved["te_region_overlap"]).sum()

    print(f"  Same TE region (overlapping): {same_region:,}")
    print(f"  Different TE region: {diff_region:,}")

    # Add gene symbols from ortholog table
    fbgn_to_symbol = dict(zip(orthologs["dmel_fbgn"], orthologs["dmel_symbol"]))
    conserved["dmel_symbol"] = conserved["gene_id"].map(fbgn_to_symbol)

    # Parse TE metadata from sseqid header if available
    # (RM_N format — we'd need the FASTA headers for repeat_name)

    # Sort by overlap size
    conserved = conserved.sort_values("overlap_bp", ascending=False)

    # Select output columns
    out_cols = [
        "gene_id", "dmel_symbol", "sseqid",
        "te_region_overlap", "overlap_bp",
        "dmel_sstart", "dmel_send", "dmel_pident", "dmel_evalue",
        "dpse_sstart", "dpse_send", "dpse_pident", "dpse_evalue",
        "dmel_qseqid", "dpse_qseqid",
        "dmel_nhits", "dpse_nhits",
    ]
    conserved[out_cols].to_csv(
        out_dir / "conserved_te_hits.tsv", sep="\t", index=False,
    )

    # Summary stats
    n_genes = conserved["gene_id"].nunique()
    n_tes = conserved["sseqid"].nunique()
    print(f"\nSummary:")
    print(f"  Unique genes with conserved TE hits: {n_genes}")
    print(f"  Unique TE instances involved: {n_tes}")
    print(f"  Mean dmel pident: {conserved['dmel_pident'].mean():.1f}%")
    print(f"  Mean dpse pident: {conserved['dpse_pident'].mean():.1f}%")

    # Also save same-region subset
    same = conserved[conserved["te_region_overlap"]]
    if not same.empty:
        same[out_cols].to_csv(
            out_dir / "conserved_te_hits_same_region.tsv", sep="\t", index=False,
        )
        print(f"\n  Same-region hits saved ({len(same):,} rows)")
        print(f"  Mean overlap: {same['overlap_bp'].mean():.0f} bp")

    print(f"\nOutput: {out_dir / 'conserved_te_hits.tsv'}")


if __name__ == "__main__":
    main()
