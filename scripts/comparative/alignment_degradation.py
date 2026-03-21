#!/usr/bin/env python3
"""Compare alignment quality for shared TE hits across dmel and dpse.

For TE instances hit by BOTH the dmel UTR and its orthologous dpse UTR
(at the same TE region), compare the aligned sequences to detect degradation:

  - Does the dpse alignment show lower percent identity than dmel?
  - Are different positions mutated in the two lineages?
  - What's the pattern: uniform decay, or specific regions diverged?

This reveals whether conserved TE fossils are under purifying selection
(both lineages maintain the sequence) or drifting independently.

Uses dmel self-BLAST and dpse->dmel cross-BLAST (same TE db).

Output: results/comparative/analysis/alignment_degradation.tsv
"""

import sys
from pathlib import Path

import numpy as np
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


def alignment_diff(seq1: str, seq2: str) -> dict:
    """Compare two aligned sequences against the same subject.

    Both seq1 and seq2 are query sequences aligned to the same TE region.
    We compare the mutation patterns.
    """
    s1 = str(seq1).upper().replace("-", "")
    s2 = str(seq2).upper().replace("-", "")
    min_len = min(len(s1), len(s2))
    if min_len == 0:
        return {"shared_len": 0, "matches": 0, "mismatches": 0, "match_rate": 0.0}

    matches = sum(1 for i in range(min_len) if s1[i] == s2[i])
    return {
        "shared_len": min_len,
        "matches": matches,
        "mismatches": min_len - matches,
        "match_rate": matches / min_len,
    }


def compare_subject_alignments(dmel_sseq: str, dpse_sseq: str) -> dict:
    """Compare subject-side alignments to detect TE region divergence.

    If both species hit the same TE region, the subject sequence should
    be identical (it's the same TE instance). Differences in the subject
    alignment reflect gaps/indels in different query species.
    """
    s1 = str(dmel_sseq).upper()
    s2 = str(dpse_sseq).upper()
    min_len = min(len(s1), len(s2))
    if min_len == 0:
        return {"sseq_shared_len": 0, "sseq_match_rate": 0.0}

    matches = sum(1 for i in range(min_len) if s1[i] == s2[i])
    return {
        "sseq_shared_len": min_len,
        "sseq_match_rate": matches / min_len,
    }


def main():
    out_dir = RESULTS_DIR / "analysis"
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading BLAST results...")
    dmel_self = load_blast(RESULTS_DIR / "dm6" / "blast_results.tsv")
    dpse_cross = load_blast(RESULTS_DIR / "UCI_Dpse_MV25_vs_dm6_te" / "blast_results.tsv")

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

    # For each (dmel_gene, TE), take the best dmel hit and best dpse hit
    # "Best" = highest bitscore
    print("Selecting best hits per (gene, TE) pair...")

    dmel_best = (
        dmel_self[dmel_self["gene_id"].isin(dmel_to_dpse)]
        .sort_values("bitscore", ascending=False)
        .drop_duplicates(subset=["gene_id", "sseqid"], keep="first")
    )

    dpse_best = (
        dpse_with_ortho
        .sort_values("bitscore", ascending=False)
        .drop_duplicates(subset=["dmel_ortholog", "sseqid"], keep="first")
    )

    # Join: shared (gene, TE) pairs where TE regions overlap
    print("Finding overlapping TE region hits...")
    merged = dmel_best.merge(
        dpse_best,
        left_on=["gene_id", "sseqid"],
        right_on=["dmel_ortholog", "sseqid"],
        suffixes=("_dmel", "_dpse"),
        how="inner",
    )

    if merged.empty:
        print("No shared hits found.")
        return

    # Filter to overlapping TE regions
    overlap = (
        (merged["sstart_dmel"] <= merged["send_dpse"]) &
        (merged["send_dmel"] >= merged["sstart_dpse"])
    )
    shared = merged[overlap].copy()

    print(f"  Total shared (gene, TE) pairs: {len(merged):,}")
    print(f"  With overlapping TE region: {len(shared):,}")

    if shared.empty:
        print("No overlapping TE region hits.")
        return

    # Compare alignments
    print("Comparing alignment quality...")

    # Direct metrics from BLAST
    shared["pident_diff"] = shared["pident_dmel"] - shared["pident_dpse"]
    shared["length_diff"] = shared["length_dmel"] - shared["length_dpse"]
    shared["evalue_ratio"] = shared["evalue_dpse"] / shared["evalue_dmel"].clip(lower=1e-180)

    # Compare query sequences (dmel qseq vs dpse qseq aligned to same TE)
    qseq_diffs = shared.apply(
        lambda row: alignment_diff(row["qseq_dmel"], row["qseq_dpse"]),
        axis=1,
        result_type="expand",
    )
    shared = pd.concat([shared, qseq_diffs], axis=1)

    # Compare subject-side alignments
    sseq_diffs = shared.apply(
        lambda row: compare_subject_alignments(row["sseq_dmel"], row["sseq_dpse"]),
        axis=1,
        result_type="expand",
    )
    shared = pd.concat([shared, sseq_diffs], axis=1)

    # Add gene symbol
    fbgn_to_symbol = dict(zip(orthologs["dmel_fbgn"], orthologs["dmel_symbol"]))
    shared["dmel_symbol"] = shared["gene_id_dmel"].map(fbgn_to_symbol)

    # Classify degradation pattern
    shared["degradation_class"] = "neutral"
    # dpse significantly worse
    shared.loc[shared["pident_diff"] > 5, "degradation_class"] = "dpse_degraded"
    # dmel significantly worse (unexpected — dpse closer to TE?)
    shared.loc[shared["pident_diff"] < -5, "degradation_class"] = "dmel_degraded"
    # Both high identity
    shared.loc[
        (shared["pident_dmel"] >= 80) & (shared["pident_dpse"] >= 80),
        "degradation_class",
    ] = "both_conserved"

    # Output
    out_cols = [
        "gene_id_dmel", "dmel_symbol", "sseqid",
        "pident_dmel", "pident_dpse", "pident_diff",
        "length_dmel", "length_dpse",
        "evalue_dmel", "evalue_dpse",
        "sstart_dmel", "send_dmel", "sstart_dpse", "send_dpse",
        "shared_len", "matches", "mismatches", "match_rate",
        "sseq_shared_len", "sseq_match_rate",
        "degradation_class",
        "qseq_dmel", "qseq_dpse", "sseq_dmel", "sseq_dpse",
    ]
    out = shared[[c for c in out_cols if c in shared.columns]]
    out.to_csv(out_dir / "alignment_degradation.tsv", sep="\t", index=False)

    # Summary
    print(f"\nAlignment Degradation Summary ({len(shared):,} shared hits):")
    print(f"  Mean dmel pident: {shared['pident_dmel'].mean():.1f}%")
    print(f"  Mean dpse pident: {shared['pident_dpse'].mean():.1f}%")
    print(f"  Mean pident difference (dmel - dpse): {shared['pident_diff'].mean():.1f}%")
    print(f"  Mean query sequence match rate: {shared['match_rate'].mean():.3f}")
    print()

    class_counts = shared["degradation_class"].value_counts()
    print("  Degradation classes:")
    for cls, count in class_counts.items():
        pct = count / len(shared) * 100
        print(f"    {cls}: {count:,} ({pct:.1f}%)")

    # Pident distribution comparison
    print(f"\n  Percent identity distribution:")
    for label, col in [("dmel", "pident_dmel"), ("dpse", "pident_dpse")]:
        vals = shared[col]
        print(f"    {label}: median={vals.median():.1f}%, "
              f"mean={vals.mean():.1f}%, "
              f"std={vals.std():.1f}%")

    # Are the most conserved hits (both_conserved) in specific genes?
    both_conserved = shared[shared["degradation_class"] == "both_conserved"]
    if not both_conserved.empty:
        top_genes = both_conserved["dmel_symbol"].value_counts().head(10)
        print(f"\n  Top genes with both-species-conserved TE hits:")
        for gene, count in top_genes.items():
            print(f"    {gene}: {count} conserved hits")

    print(f"\nOutput: {out_dir / 'alignment_degradation.tsv'}")


if __name__ == "__main__":
    main()
