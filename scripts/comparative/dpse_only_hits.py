#!/usr/bin/env python3
"""Find dmel TE hits in dpse UTRs that are ABSENT from orthologous dmel UTRs.

These are TE sequences that dpse 3'UTRs match (in the dmel TE database) but
the orthologous dmel 3'UTR does not. Possible explanations:
  - Lineage-specific TE insertion in dpse (post-divergence)
  - TE loss/deletion in the dmel lineage
  - UTR annotation differences (dmel UTR boundary doesn't cover the region)

Uses dmel self-BLAST and dpse->dmel cross-BLAST.

Output: results/comparative/analysis/dpse_only_te_hits.tsv
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
    dmel_self = load_blast(RESULTS_DIR / "dm6" / "blast_results.tsv")
    dpse_cross = load_blast(RESULTS_DIR / "UCI_Dpse_MV25_vs_dm6_te" / "blast_results.tsv")

    print(f"  dmel self: {len(dmel_self):,} hits")
    print(f"  dpse->dmel cross: {len(dpse_cross):,} hits")

    # Gene mappings
    print("Building gene mappings...")
    dmel_r2g = build_region_to_gene_dmel(
        RESULTS_DIR / "dm6" / "regions.tsv",
        DATA_DIR / "references" / "dmel_annotation_coding.gff3",
    )
    dpse_r2g = build_region_to_gene_dpse(
        RESULTS_DIR / "UCI_Dpse_MV25" / "regions.tsv",
        DATA_DIR / "references" / "dpse_annotation.gff3",
    )

    dmel_self["gene_id"] = dmel_self["qseqid"].map(dmel_r2g)
    dpse_cross["gene_id"] = dpse_cross["qseqid"].map(dpse_r2g)

    # Ortholog mapping
    orthologs = load_orthologs()
    dmel_to_dpse = build_ortholog_map(orthologs)
    dpse_to_dmel = {v: k for k, v in dmel_to_dpse.items()}

    dpse_cross["dmel_ortholog"] = dpse_cross["gene_id"].map(dpse_to_dmel)
    dpse_with_ortho = dpse_cross[dpse_cross["dmel_ortholog"].notna()].copy()

    # Build set of (dmel_gene, TE_instance) pairs from dmel self-BLAST
    dmel_gene_te = set(
        zip(dmel_self["gene_id"].values, dmel_self["sseqid"].values)
    )

    print(f"  dmel (gene, TE) pairs: {len(dmel_gene_te):,}")

    # Find dpse hits where the orthologous dmel gene does NOT hit the same TE
    print("Finding dpse-only hits...")
    dpse_with_ortho["in_dmel"] = [
        (ortho, te) in dmel_gene_te
        for ortho, te in zip(
            dpse_with_ortho["dmel_ortholog"].values,
            dpse_with_ortho["sseqid"].values,
        )
    ]

    dpse_only = dpse_with_ortho[~dpse_with_ortho["in_dmel"]].copy()
    dpse_shared = dpse_with_ortho[dpse_with_ortho["in_dmel"]].copy()

    print(f"  dpse hits in ortholog UTRs: {len(dpse_with_ortho):,}")
    print(f"  Also in dmel (shared): {len(dpse_shared):,}")
    print(f"  dpse-only: {len(dpse_only):,}")

    if dpse_only.empty:
        print("No dpse-only hits found.")
        return

    # Add gene symbol
    fbgn_to_symbol = dict(zip(orthologs["dmel_fbgn"], orthologs["dmel_symbol"]))
    dpse_only["dmel_symbol"] = dpse_only["dmel_ortholog"].map(fbgn_to_symbol)

    # Summarize
    n_genes = dpse_only["dmel_ortholog"].nunique()
    n_tes = dpse_only["sseqid"].nunique()

    # Group by gene for summary
    gene_summary = dpse_only.groupby("dmel_ortholog").agg(
        dmel_symbol=("dmel_symbol", "first"),
        n_te_hits=("sseqid", "nunique"),
        n_blast_hits=("sseqid", "size"),
        mean_pident=("pident", "mean"),
        best_evalue=("evalue", "min"),
    ).sort_values("n_te_hits", ascending=False)

    gene_summary.to_csv(
        out_dir / "dpse_only_gene_summary.tsv", sep="\t",
    )

    # Save full hit table (top hits per gene-TE pair for manageability)
    dpse_only_best = dpse_only.sort_values("evalue").drop_duplicates(
        subset=["dmel_ortholog", "sseqid"], keep="first",
    )

    out_cols = [
        "dmel_ortholog", "dmel_symbol", "gene_id", "qseqid", "sseqid",
        "pident", "length", "evalue", "bitscore",
        "qstart", "qend", "sstart", "send",
        "qseq", "sseq",
    ]
    dpse_only_best[[c for c in out_cols if c in dpse_only_best.columns]].to_csv(
        out_dir / "dpse_only_te_hits.tsv", sep="\t", index=False,
    )

    print(f"\nSummary:")
    print(f"  Genes with dpse-only TE hits: {n_genes}")
    print(f"  Unique TEs involved: {n_tes}")
    print(f"  Mean pident: {dpse_only['pident'].mean():.1f}%")
    print(f"  Median alignment length: {dpse_only['length'].median():.0f} bp")

    print(f"\n  Top genes by dpse-only TE count:")
    for _, row in gene_summary.head(10).iterrows():
        print(f"    {row['dmel_symbol'] or row.name}: "
              f"{row['n_te_hits']} TEs, "
              f"mean pident {row['mean_pident']:.1f}%")

    print(f"\nOutput: {out_dir / 'dpse_only_te_hits.tsv'}")
    print(f"        {out_dir / 'dpse_only_gene_summary.tsv'}")


if __name__ == "__main__":
    main()
