#!/usr/bin/env python3
"""Run all 7 analysis workstreams on dmel BLAST results.

Grain of salt: evalue=10, word_size=7, dust=no — very permissive params.
98% of genes show TE hits. Quality tiers are critical for interpretation.

Execution order per ANALYSIS_ARCHITECTURE.md:
  WS1 → (WS3, WS5, WS7 parallel) → WS4 → Report
  (WS6 conservation and WS2 motifs deferred — need external data)

Usage:
    python scripts/run_full_analysis.py
"""

import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# ── Setup ────────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from fossil_finder.io.blast import BLAST_COLUMNS, load_blast_results, classify_strand

RESULTS_DIRS = {
    "3utr": PROJECT_ROOT / "results" / "dmel_3utr_e10",
    "5utr": PROJECT_ROOT / "results" / "dmel_5utr_e10",
}

GENE_LISTS_DIR = PROJECT_ROOT / "data" / "gene_lists"
GENE_SETS = {
    "germ_plasm": GENE_LISTS_DIR / "germ_plasm_fbgn_ids.txt",
    "housekeeping": GENE_LISTS_DIR / "housekeeping_fbgn_ids.txt",
    "somatic": GENE_LISTS_DIR / "somatic_fbgn_ids.txt",
    "cleared": GENE_LISTS_DIR / "cleared_fbgn_ids.txt",
    "adult": GENE_LISTS_DIR / "adult_fbgn_ids.txt",
}
DATA_DIR = PROJECT_ROOT / "data"
GFF_PATH = DATA_DIR / "references" / "dmel-all-r6.66.gff"
TE_FASTA = DATA_DIR / "references" / "dmel_te_flybase.fasta"
TE_CONSENSUS = DATA_DIR / "references" / "dmel_te_consensus.fasta"

# Quality tier thresholds
STRICT = {"pident": 85, "length": 100}
MODERATE = {"pident": 75, "length": 50}


def load_gene_set(path):
    """Load gene IDs from file, one per line, skip comments."""
    genes = set()
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                genes.add(line.split("\t")[0].split()[0])
    return genes


def build_transcript_to_gene_map():
    """Build FBtr → FBgn mapping from GFF."""
    if not GFF_PATH.exists():
        print(f"  WARNING: GFF not found at {GFF_PATH}")
        return {}
    t2g = {}
    with open(GFF_PATH) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 9 or parts[2] != "mRNA":
                continue
            attrs = parts[8]
            tid = gid = None
            for pair in attrs.split(";"):
                if pair.strip().startswith("ID="):
                    tid = pair.strip().split("=", 1)[1]
                elif pair.strip().startswith("Parent="):
                    gid = pair.strip().split("=", 1)[1]
            if tid and gid:
                t2g[tid] = gid
    return t2g


def build_query_to_gene(regions_path, t2g):
    """Build region_id → FBgn mapping."""
    q2g = {}
    regions = pd.read_csv(regions_path, sep="\t")
    for _, r in regions.iterrows():
        region_id = r["region_id"]
        parents = [p.strip() for p in str(r["parent_id"]).split(",")]
        for parent in parents:
            gene_id = t2g.get(parent, parent)
            if gene_id.startswith("FBgn"):
                q2g[region_id] = gene_id
                break
    return q2g


def build_consensus_class_map():
    """Build family_name → TE class from consensus FASTA headers.

    Consensus headers use RepBase format: >family_name#class/subclass
    e.g. >gypsy#LTR/Gypsy, >jockey#LINE/Jockey, >hobo#DNA/hAT
    Returns {lowercase_name: top_level_class} e.g. {"gypsy": "LTR", "hobo": "DNA"}
    """
    class_map = {}
    if not TE_CONSENSUS.exists():
        print(f"  WARNING: Consensus FASTA not found at {TE_CONSENSUS}")
        return class_map
    with open(TE_CONSENSUS) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            if "#" not in header:
                continue
            name, classification = header.split("#", 1)
            top_class = classification.split("/")[0]
            # Normalize RC (rolling circle) → Helitron
            if top_class == "RC":
                top_class = "Helitron"
            class_map[name.lower()] = top_class
    return class_map


def build_te_name_map():
    """Build FBti → family name mapping from TE FASTA headers."""
    te_map = {}
    if not TE_FASTA.exists():
        print(f"  WARNING: TE FASTA not found at {TE_FASTA}")
        return te_map
    import re
    with open(TE_FASTA) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            te_id = line.split()[0][1:]
            name_match = re.search(r"name=([^;]+)", line)
            if name_match:
                te_map[te_id] = name_match.group(1).strip()
    return te_map


def _strip_instance_suffix(name):
    """Strip FlyBase instance suffix from TE family name.

    FlyBase uses: family{}instance e.g. 'gypsy{}123', '1360{}Eph[1520]'
    We strip everything from '{}' onward to get the base family name.
    """
    idx = name.find("{}")
    return name[:idx] if idx >= 0 else name


def infer_te_class_from_name(te_name, consensus_map=None):
    """Classify TE using consensus FASTA classes, with pattern fallback.

    Strategy:
    1. Strip {}N instance suffix from name
    2. Exact match against consensus map (case-insensitive)
    3. Try name + "-element" suffix (FlyBase uses "Doc2" but consensus has "Doc2-element")
    4. Fall back to keyword pattern matching
    """
    base_name = _strip_instance_suffix(te_name)
    base_lower = base_name.lower()

    if consensus_map:
        # Exact match
        if base_lower in consensus_map:
            return consensus_map[base_lower]
        # Try with -element suffix
        if base_lower + "-element" in consensus_map:
            return consensus_map[base_lower + "-element"]

    # Known Drosophila families not in consensus FASTA
    # (verified against FlyBase/RepBase annotations)
    FLYBASE_OVERRIDES = {
        "h": "DNA",            # HB-related hAT superfamily
        "het-tag": "LINE",     # telomeric, HeT-A-related
        "antonia": "LTR",      # LTR/Gypsy
        "ninja-dsim-like": "LTR",  # LTR/Gypsy (ninja family)
        "y": "LINE",           # telomeric Y-element
    }
    if base_lower in FLYBASE_OVERRIDES:
        return FLYBASE_OVERRIDES[base_lower]

    # Pattern-based fallback
    tn = base_lower
    if any(k in tn for k in ["gypsy", "copia", "bel", "pao", "mdg", "tirant",
                               "blood", "idefix", "zam", "tabor", "springer",
                               "invader", "opus", "roo", "micropia", "17.6",
                               "297", "412", "1731"]):
        return "LTR"
    elif any(k in tn for k in ["jockey", "cr1", "r1_", "r2_", "rte",
                                 "doc", "f_element", "f-element", "g_element",
                                 "g-element", "i_element", "i-element",
                                 "het-a", "tart", "tahre", "baggins",
                                 "juan", "line", "x-element"]):
        return "LINE"
    elif any(k in tn for k in ["mariner", "tc1", "hobo", "pogo", "p_element",
                                 "p-element", "piggyback", "transib", "bari",
                                 "s_element", "s-element", "1360", "tc3",
                                 "looper", "dna"]):
        return "DNA"
    elif any(k in tn for k in ["helitron", "dine", "ine-1"]):
        return "Helitron"
    elif any(k in tn for k in ["sine", "alu"]):
        return "SINE"
    elif any(k in tn for k in ["ltr", "gag", "env"]):
        return "LTR"
    else:
        return "Unknown"


def assign_quality_tier(df):
    """Add 'tier' column: strict / moderate / relaxed."""
    strict_mask = (df["pident"] >= STRICT["pident"]) & (df["length"] >= STRICT["length"])
    moderate_mask = (df["pident"] >= MODERATE["pident"]) & (df["length"] >= MODERATE["length"])
    df["tier"] = "relaxed"
    df.loc[moderate_mask, "tier"] = "moderate"
    df.loc[strict_mask, "tier"] = "strict"
    return df


def compute_edit_stats(df):
    """Compute per-hit edit distance stats from qseq/sseq columns (vectorized)."""
    if "qseq" not in df.columns or "sseq" not in df.columns:
        return df

    # Use vectorized mismatch/gap counting from the BLAST columns directly
    # mismatch and gapopen are already in the BLAST output
    df["mismatch_count"] = df["mismatch"]
    df["gap_count"] = df["gapopen"]
    df["edit_distance"] = df["mismatch_count"] + df["gap_count"]
    df["mismatch_rate"] = df["mismatch_count"] / df["length"].clip(lower=1)
    df["gap_rate"] = df["gap_count"] / df["length"].clip(lower=1)

    # Ti/Tv from a small subsample (character-level analysis is expensive)
    n_titv = min(10_000, len(df))
    subsample = df.sample(n=n_titv, random_state=42)
    purines = set("AGag")
    pyrimidines = set("CTct")
    total_ti = total_tv = 0
    for _, row in subsample[["qseq", "sseq"]].iterrows():
        qs, ss = str(row["qseq"]), str(row["sseq"])
        for q, s in zip(qs, ss):
            if q != s and q != "-" and s != "-":
                if (q in purines and s in purines) or (q in pyrimidines and s in pyrimidines):
                    total_ti += 1
                else:
                    total_tv += 1
    global_ti_tv = total_ti / max(total_tv, 1)
    df["ti_tv_ratio"] = global_ti_tv  # approximate, same for all rows
    return df


def compute_positional_stats(df):
    """Add normalized UTR position and TE position columns."""
    # UTR position: where on the UTR does the hit midpoint fall?
    hit_midpoint = (df["qstart"] + df["qend"]) / 2
    df["normalized_utr_pos"] = hit_midpoint / df["qlen"].clip(lower=1)
    df["normalized_utr_pos"] = df["normalized_utr_pos"].clip(0, 1)

    # TE position: where on the TE element?
    te_mid = (df["sstart"].abs() + df["send"].abs()) / 2  # abs for minus strand
    df["te_normalized_mid"] = te_mid / df["slen"].clip(lower=1)
    df["te_normalized_mid"] = df["te_normalized_mid"].clip(0, 1)
    return df


# ── WS1: Match Distributions & Patterns ─────────────────────────────────────

def run_ws1(df, utr_type, out_dir, q2g, te_names):
    """WS1: Quality tiers, edit stats, positional profiles, TE breadth."""
    print(f"\n{'='*60}")
    print(f"WS1: Match Distributions — {utr_type}")
    print(f"{'='*60}")
    t0 = time.time()

    # 1.1 Quality tiers
    df = assign_quality_tier(df)
    tier_counts = df["tier"].value_counts()
    tier_pcts = (tier_counts / len(df) * 100).round(2)
    print(f"\nQuality tier distribution ({len(df):,} total hits):")
    for tier in ["strict", "moderate", "relaxed"]:
        n = tier_counts.get(tier, 0)
        pct = tier_pcts.get(tier, 0)
        print(f"  {tier:>10}: {n:>10,} ({pct:>5.1f}%)")

    # 1.1b Alignment quality summary per tier
    tier_summary = []
    for tier in ["strict", "moderate", "relaxed"]:
        subset = df[df["tier"] == tier]
        if len(subset) == 0:
            continue
        tier_summary.append({
            "tier": tier,
            "n_hits": len(subset),
            "mean_pident": subset["pident"].mean(),
            "median_pident": subset["pident"].median(),
            "mean_length": subset["length"].mean(),
            "median_length": subset["length"].median(),
            "mean_evalue": subset["evalue"].mean(),
            "median_evalue": subset["evalue"].median(),
            "mean_bitscore": subset["bitscore"].mean(),
            "q25_pident": subset["pident"].quantile(0.25),
            "q75_pident": subset["pident"].quantile(0.75),
            "q25_length": subset["length"].quantile(0.25),
            "q75_length": subset["length"].quantile(0.75),
        })
    tier_df = pd.DataFrame(tier_summary)
    tier_df.to_csv(out_dir / "alignment_stats_summary.tsv", sep="\t", index=False)

    # 1.2 Edit distance (sample for speed — full 7.5M is slow for string ops)
    n_sample = min(500_000, len(df))
    sample = df.sample(n=n_sample, random_state=42).copy()
    sample = compute_edit_stats(sample)

    if "mismatch_rate" in sample.columns:
        edit_summary = []
        for tier in ["strict", "moderate", "relaxed"]:
            s = sample[sample["tier"] == tier]
            if len(s) == 0:
                continue
            edit_summary.append({
                "tier": tier,
                "n_sampled": len(s),
                "mean_mismatch_rate": s["mismatch_rate"].mean(),
                "median_mismatch_rate": s["mismatch_rate"].median(),
                "mean_gap_rate": s["gap_rate"].mean(),
                "mean_edit_distance": s["edit_distance"].mean(),
                "mean_ti_tv_ratio": s["ti_tv_ratio"].dropna().mean(),
                "median_ti_tv_ratio": s["ti_tv_ratio"].dropna().median(),
            })
        edit_df = pd.DataFrame(edit_summary)
        edit_df.to_csv(out_dir / "edit_distance_summary.tsv", sep="\t", index=False)
        print(f"\nEdit distance (sampled {n_sample:,}):")
        for _, row in edit_df.iterrows():
            print(f"  {row['tier']:>10}: mismatch_rate={row['mean_mismatch_rate']:.3f}, "
                  f"gap_rate={row['mean_gap_rate']:.3f}, ti/tv={row['mean_ti_tv_ratio']:.2f}")

    # 1.3 Positional bias on UTR
    df = compute_positional_stats(df)
    df["utr_decile"] = (df["normalized_utr_pos"] * 10).astype(int).clip(0, 9)

    pos_profiles = []
    for tier in ["strict", "moderate", "relaxed"]:
        for decile in range(10):
            mask = (df["tier"] == tier) & (df["utr_decile"] == decile)
            n = mask.sum()
            pos_profiles.append({
                "decile": decile,
                "decile_label": f"{decile*10}-{(decile+1)*10}%",
                "tier": tier,
                "utr_type": utr_type,
                "hit_count": int(n),
            })
    pos_df = pd.DataFrame(pos_profiles)
    # Normalize within each tier
    for tier in pos_df["tier"].unique():
        mask = pos_df["tier"] == tier
        total = pos_df.loc[mask, "hit_count"].sum()
        pos_df.loc[mask, "normalized_density"] = pos_df.loc[mask, "hit_count"] / max(total / 10, 1)
    pos_df.to_csv(out_dir / "positional_profiles_utr.tsv", sep="\t", index=False)

    # Report 3' end enrichment
    for tier in ["strict", "moderate", "relaxed"]:
        t = pos_df[pos_df["tier"] == tier]
        if len(t) == 0:
            continue
        last_decile = t[t["decile"] == 9]["hit_count"].sum()
        first_decile = t[t["decile"] == 0]["hit_count"].sum()
        ratio = last_decile / max(first_decile, 1)
        print(f"  {tier} 3'-end/5'-end ratio: {ratio:.2f}x")

    # 1.4 TE positional pileup
    df["te_decile"] = (df["te_normalized_mid"] * 10).astype(int).clip(0, 9)
    top_tes = df["sseqid"].value_counts().head(100).index
    te_pos_profiles = []
    for te_id in top_tes:
        te_hits = df[df["sseqid"] == te_id]
        for decile in range(10):
            n = (te_hits["te_decile"] == decile).sum()
            te_pos_profiles.append({
                "te_id": te_id,
                "decile": decile,
                "hit_count": int(n),
                "coverage_fraction": n / max(len(te_hits), 1),
            })
    te_pos_df = pd.DataFrame(te_pos_profiles)
    te_pos_df.to_csv(out_dir / "positional_profiles_te.tsv", sep="\t", index=False)

    # 1.5 Per-TE unique hit counts & gene breadth
    df["gene_id"] = df["qseqid"].map(q2g)
    mapped = df["gene_id"].notna().sum()
    print(f"\n  Gene mapping: {mapped:,}/{len(df):,} hits mapped ({mapped/len(df)*100:.1f}%)")

    te_breadth = df.groupby("sseqid").agg(
        n_hits=("qseqid", "count"),
        n_unique_queries=("qseqid", "nunique"),
        n_unique_genes=("gene_id", "nunique"),
        mean_pident=("pident", "mean"),
        mean_length=("length", "mean"),
        total_bp=("length", "sum"),
        mean_evalue=("evalue", "mean"),
    ).sort_values("n_hits", ascending=False)

    # Gene entropy per TE — compute from pre-aggregated counts
    # (avoid slow .apply over 7.5M rows)
    gene_te_counts = df.dropna(subset=["gene_id"]).groupby(["sseqid", "gene_id"]).size()
    te_entropies = {}
    for te_id, group in gene_te_counts.groupby(level=0):
        counts = group.values
        if len(counts) <= 1:
            te_entropies[te_id] = 0.0
        else:
            probs = counts / counts.sum()
            te_entropies[te_id] = float(stats.entropy(probs))
    te_breadth["gene_entropy"] = pd.Series(te_entropies)
    te_breadth.to_csv(out_dir / "te_hit_breadth.tsv", sep="\t", index=True)

    print(f"\nTE breadth stats ({len(te_breadth):,} TEs):")
    print(f"  Top TE: {te_breadth.index[0]} ({te_breadth.iloc[0]['n_hits']:,} hits, "
          f"{te_breadth.iloc[0]['n_unique_genes']:.0f} genes)")
    print(f"  Median hits/TE: {te_breadth['n_hits'].median():.0f}")
    print(f"  TEs hitting >100 genes: {(te_breadth['n_unique_genes'] >= 100).sum()}")

    # 1.6 Hit multiplicity
    hits_per_query = df.groupby("qseqid")["sseqid"].count()
    hits_per_gene = df.groupby("gene_id")["sseqid"].count() if "gene_id" in df.columns else pd.Series(dtype=int)
    queries_per_te = df.groupby("sseqid")["qseqid"].nunique()
    genes_per_te = df.groupby("sseqid")["gene_id"].nunique()

    mult_rows = []
    for name, series in [("hits_per_query", hits_per_query), ("hits_per_gene", hits_per_gene),
                         ("queries_per_te", queries_per_te), ("genes_per_te", genes_per_te)]:
        if len(series) == 0:
            continue
        mult_rows.append({
            "metric": name,
            "mean": series.mean(),
            "median": series.median(),
            "q25": series.quantile(0.25),
            "q75": series.quantile(0.75),
            "max": series.max(),
            "n": len(series),
        })
    mult_df = pd.DataFrame(mult_rows)
    mult_df.to_csv(out_dir / "hit_multiplicity.tsv", sep="\t", index=False)

    print(f"\nHit multiplicity:")
    for _, row in mult_df.iterrows():
        print(f"  {row['metric']}: median={row['median']:.0f}, "
              f"mean={row['mean']:.0f}, max={row['max']:.0f}")

    elapsed = time.time() - t0
    print(f"\nWS1 complete in {elapsed:.1f}s")
    return df


# ── WS3: TE Family Enrichment ────────────────────────────────────────────────

def run_ws3(df, gene_sets_loaded, out_dir, te_names, consensus_map=None):
    """WS3: TE family ranking, gene-set enrichment, class distribution."""
    print(f"\n{'='*60}")
    print(f"WS3: TE Family Enrichment")
    print(f"{'='*60}")
    t0 = time.time()

    # 3.1 Global family ranking per tier
    rankings = []
    # Pre-compute strand indicators for fast aggregation
    df["is_sense"] = (df["strand"] == "plus").astype(int)
    for tier in ["strict", "moderate", "relaxed"]:
        subset = df[df["tier"] == tier]
        if len(subset) == 0:
            continue
        family_stats = subset.groupby("sseqid").agg(
            n_hits=("qseqid", "count"),
            total_bp=("length", "sum"),
            n_unique_genes=("gene_id", "nunique"),
            mean_pident=("pident", "mean"),
            mean_evalue=("evalue", "mean"),
            sense_hits=("is_sense", "sum"),
        ).reset_index()
        family_stats["antisense_hits"] = family_stats["n_hits"] - family_stats["sense_hits"]
        family_stats["tier"] = tier
        rankings.append(family_stats)

    if rankings:
        ranking_df = pd.concat(rankings, ignore_index=True)
        ranking_df = ranking_df.sort_values(["tier", "n_hits"], ascending=[True, False])
        ranking_df.to_csv(out_dir / "family_ranking_global.tsv", sep="\t", index=False)

        # Report top 5 per tier
        for tier in ["strict", "moderate", "relaxed"]:
            t = ranking_df[ranking_df["tier"] == tier].head(5)
            print(f"\n  Top 5 TEs ({tier}):")
            for _, row in t.iterrows():
                print(f"    {row['sseqid']}: {row['n_hits']:,} hits, "
                      f"{row['n_unique_genes']:.0f} genes, "
                      f"mean_pident={row['mean_pident']:.1f}")

    # 3.2 Family enrichment by gene set (vectorized)
    if gene_sets_loaded and "gene_id" in df.columns:
        all_genes = set(df["gene_id"].dropna().unique())
        enrichment_frames = []
        for set_name, set_genes in gene_sets_loaded.items():
            set_genes_in_data = set_genes & all_genes
            bg_genes = all_genes - set_genes_in_data
            if not set_genes_in_data:
                continue

            set_mask = df["gene_id"].isin(set_genes_in_data)
            bg_mask = df["gene_id"].isin(bg_genes)

            set_counts = df.loc[set_mask, "sseqid"].value_counts()
            bg_counts = df.loc[bg_mask, "sseqid"].value_counts()

            # Align indices
            all_tes_idx = set_counts.index.union(bg_counts.index)
            set_counts = set_counts.reindex(all_tes_idx, fill_value=0)
            bg_counts = bg_counts.reindex(all_tes_idx, fill_value=0)

            set_freq = set_counts / max(len(set_genes_in_data), 1)
            bg_freq = bg_counts / max(len(bg_genes), 1)
            fold = set_freq / bg_freq.replace(0, 1e-10)

            enrich = pd.DataFrame({
                "te_id": all_tes_idx,
                "gene_set": set_name,
                "count_in_set": set_counts.values,
                "count_in_bg": bg_counts.values,
                "freq_in_set": set_freq.values,
                "freq_in_bg": bg_freq.values,
                "fold_enrichment": fold.values,
                "log2_fold": np.log2(fold.replace(0, 1e-10)).values,
            })
            enrichment_frames.append(enrich)

        if enrichment_frames:
            enrich_df = pd.concat(enrichment_frames, ignore_index=True)
            enrich_df.to_csv(out_dir / "family_enrichment_by_set.tsv", sep="\t", index=False)

            for set_name in gene_sets_loaded:
                top = enrich_df[(enrich_df["gene_set"] == set_name) &
                                (enrich_df["count_in_set"] >= 5)].nlargest(3, "fold_enrichment")
                if len(top) > 0:
                    print(f"\n  Top enriched TEs in {set_name}:")
                    for _, row in top.iterrows():
                        print(f"    {row['te_id']}: {row['fold_enrichment']:.1f}x "
                              f"(n={row['count_in_set']:.0f})")

    # 3.3 TE class distribution (use consensus FASTA for authoritative classification)
    def get_te_class(te_id):
        name = te_names.get(te_id, te_id)
        return infer_te_class_from_name(name, consensus_map)

    df["te_class"] = df["sseqid"].map(get_te_class)
    df["te_name"] = df["sseqid"].map(lambda x: te_names.get(x, x))
    class_dist = []
    for tier in ["strict", "moderate", "relaxed"]:
        subset = df[df["tier"] == tier]
        for cls in ["LTR", "LINE", "DNA", "Helitron", "SINE", "Unknown"]:
            n = (subset["te_class"] == cls).sum()
            class_dist.append({
                "te_class": cls,
                "tier": tier,
                "n_hits": int(n),
                "pct_of_tier": n / max(len(subset), 1) * 100,
            })
    class_df = pd.DataFrame(class_dist)
    class_df.to_csv(out_dir / "class_distribution_by_tier.tsv", sep="\t", index=False)

    print(f"\nTE class distribution (all tiers):")
    for cls in ["LTR", "LINE", "DNA", "Helitron", "SINE", "Unknown"]:
        n = (df["te_class"] == cls).sum()
        print(f"  {cls:>10}: {n:>10,} ({n/len(df)*100:.1f}%)")

    # 3.4 Family x strand bias
    strand_by_te = df.groupby("sseqid").agg(
        sense=("is_sense", "sum"),
        antisense=("is_sense", lambda x: (x == 0).sum()),
    )
    strand_by_te["total"] = strand_by_te["sense"] + strand_by_te["antisense"]
    strand_by_te["sense_pct"] = strand_by_te["sense"] / strand_by_te["total"].clip(lower=1) * 100

    # Binomial test for strand bias (min 20 hits) — vectorized via binom_test
    significant = strand_by_te[strand_by_te["total"] >= 20].copy()
    if len(significant) > 0:
        # Use scipy.stats.binom.sf for vectorized p-values (two-sided)
        from scipy.stats import binom
        k = significant["sense"].values.astype(int)
        n = significant["total"].values.astype(int)
        # Two-sided: 2 * min(P(X >= k), P(X <= k))
        p_upper = binom.sf(k - 1, n, 0.5)
        p_lower = binom.cdf(k, n, 0.5)
        significant["p_value"] = 2 * np.minimum(p_upper, p_lower)
        significant["p_value"] = significant["p_value"].clip(upper=1.0)

        from statsmodels.stats.multitest import multipletests
        _, fdr, _, _ = multipletests(significant["p_value"], method="fdr_bh")
        significant["fdr"] = fdr
        sig_biased = significant[significant["fdr"] < 0.05]
        print(f"\n  TEs with significant strand bias (FDR<0.05): {len(sig_biased)}/{len(significant)}")
        significant.to_csv(out_dir / "family_strand_bias.tsv", sep="\t", index=True)

    elapsed = time.time() - t0
    print(f"\nWS3 complete in {elapsed:.1f}s")
    return df


# ── WS5: RepeatMasker Overlap ────────────────────────────────────────────────

def run_ws5(df, out_dir):
    """WS5: Cross-reference BLAST hits with RepeatMasker annotations."""
    print(f"\n{'='*60}")
    print(f"WS5: RepeatMasker Overlap")
    print(f"{'='*60}")

    rm_path = DATA_DIR / "references" / "dm6.fa.out"
    if not rm_path.exists():
        print("  RepeatMasker file not found — skipping WS5")
        print(f"  Expected: {rm_path}")
        return df

    t0 = time.time()

    # Load RM (skip header lines)
    rm_lines = []
    with open(rm_path) as f:
        for i, line in enumerate(f):
            if i < 3:  # skip header
                continue
            parts = line.split()
            if len(parts) < 15:
                continue
            try:
                # Strip 'chr' prefix to match our coordinate system (2L not chr2L)
                chrom = parts[4]
                if chrom.startswith("chr"):
                    chrom = chrom[3:]
                rm_lines.append({
                    "rm_chrom": chrom,
                    "rm_start": int(parts[5]),
                    "rm_end": int(parts[6]),
                    "rm_name": parts[9],
                    "rm_class": parts[10],
                })
            except (ValueError, IndexError):
                continue

    rm_df = pd.DataFrame(rm_lines)
    print(f"  Loaded {len(rm_df):,} RepeatMasker annotations")

    # Build interval lookup per chrom
    regions_path = out_dir.parent / "regions.tsv"
    if not regions_path.exists():
        print("  No regions.tsv — skipping overlap")
        return df

    regions = pd.read_csv(regions_path, sep="\t")

    # For each BLAST hit, check if it overlaps a RM annotation
    # Map query -> genomic coords
    region_coords = {}
    for _, r in regions.iterrows():
        region_coords[r["region_id"]] = (r["chrom"], int(r["start"]), int(r["end"]))

    # Build RM interval set per chrom using 1kb bins for fast lookup
    BIN_SIZE = 1000
    rm_bins = {}  # chrom -> bin_idx -> [(start, end, name, class)]
    for _, r in rm_df.iterrows():
        chrom = r["rm_chrom"]
        if chrom not in rm_bins:
            rm_bins[chrom] = {}
        start_bin = r["rm_start"] // BIN_SIZE
        end_bin = r["rm_end"] // BIN_SIZE
        for b in range(start_bin, end_bin + 1):
            if b not in rm_bins[chrom]:
                rm_bins[chrom][b] = []
            rm_bins[chrom][b].append((r["rm_start"], r["rm_end"], r["rm_name"], r["rm_class"]))

    # Sample for speed
    n_check = min(50_000, len(df))
    sample = df.sample(n=n_check, random_state=42).copy()

    known_count = 0
    novel_count = 0
    rm_class_counter = Counter()

    for _, hit in sample.iterrows():
        qid = hit["qseqid"]
        if qid not in region_coords:
            novel_count += 1
            continue

        chrom, reg_start, reg_end = region_coords[qid]
        hit_start_genomic = reg_start + int(hit["qstart"])
        hit_end_genomic = reg_start + int(hit["qend"])

        found = False
        if chrom in rm_bins:
            start_bin = hit_start_genomic // BIN_SIZE
            end_bin = hit_end_genomic // BIN_SIZE
            for b in range(start_bin, end_bin + 1):
                for rm_start, rm_end, rm_name, rm_class in rm_bins[chrom].get(b, []):
                    if hit_start_genomic < rm_end and hit_end_genomic > rm_start:
                        found = True
                        rm_class_counter[rm_class] += 1
                        break
                if found:
                    break

        if found:
            known_count += 1
        else:
            novel_count += 1

    total = known_count + novel_count
    print(f"\n  Sampled {n_check:,} hits:")
    print(f"    Known (overlaps RM): {known_count:,} ({known_count/total*100:.1f}%)")
    print(f"    Novel (no RM match): {novel_count:,} ({novel_count/total*100:.1f}%)")

    if rm_class_counter:
        print(f"\n  RM class breakdown (known hits):")
        for cls, count in rm_class_counter.most_common(10):
            print(f"    {cls}: {count:,}")

    rm_stats = {
        "n_sampled": n_check,
        "known_hits": known_count,
        "novel_hits": novel_count,
        "known_pct": known_count / total * 100,
        "rm_class_counts": dict(rm_class_counter.most_common()),
    }
    with open(out_dir / "ws5_rm_overlap.json", "w") as f:
        json.dump(rm_stats, f, indent=2)

    elapsed = time.time() - t0
    print(f"\nWS5 complete in {elapsed:.1f}s")
    return df


# ── WS7: GO & Functional Enrichment ─────────────────────────────────────────

def run_ws7(df, gene_sets_loaded, out_dir):
    """WS7: Gene set comparisons — density, hit counts, strand bias."""
    print(f"\n{'='*60}")
    print(f"WS7: Functional Enrichment & Gene Set Comparisons")
    print(f"{'='*60}")
    t0 = time.time()

    # Load gene stats
    gene_stats_path = out_dir.parent / "gene_stats.tsv"
    if not gene_stats_path.exists():
        print("  No gene_stats.tsv — skipping WS7")
        return

    gs = pd.read_csv(gene_stats_path, sep="\t", index_col=0)
    all_genes = set(gs.index)

    results = {}
    for set_name, set_genes in gene_sets_loaded.items():
        in_set = set_genes & all_genes
        out_set = all_genes - in_set
        if not in_set:
            continue

        in_stats = gs.loc[gs.index.isin(in_set)]
        out_stats = gs.loc[gs.index.isin(out_set)]

        # Mann-Whitney on density
        if len(in_stats) > 0 and len(out_stats) > 0:
            try:
                u_stat, mw_p = stats.mannwhitneyu(
                    in_stats["density"], out_stats["density"], alternative="two-sided"
                )
            except ValueError:
                u_stat, mw_p = 0, 1.0
        else:
            u_stat, mw_p = 0, 1.0

        # Fisher's exact on TE-positive
        te_pos_in = (in_stats["hit_count"] > 0).sum()
        te_neg_in = len(in_stats) - te_pos_in
        te_pos_out = (out_stats["hit_count"] > 0).sum()
        te_neg_out = len(out_stats) - te_pos_out
        try:
            odds, fisher_p = stats.fisher_exact([[te_pos_in, te_neg_in],
                                                  [te_pos_out, te_neg_out]])
        except ValueError:
            odds, fisher_p = 1.0, 1.0

        # Quality tier breakdown for this gene set
        set_hits = df[df["gene_id"].isin(in_set)] if "gene_id" in df.columns else pd.DataFrame()
        tier_breakdown = {}
        if len(set_hits) > 0 and "tier" in set_hits.columns:
            tc = set_hits["tier"].value_counts()
            for tier in ["strict", "moderate", "relaxed"]:
                tier_breakdown[f"{tier}_hits"] = int(tc.get(tier, 0))

        results[set_name] = {
            "n_genes_in_set": len(in_set),
            "n_genes_in_data": len(in_stats),
            "mean_density_set": float(in_stats["density"].mean()),
            "mean_density_bg": float(out_stats["density"].mean()),
            "median_density_set": float(in_stats["density"].median()),
            "median_density_bg": float(out_stats["density"].median()),
            "density_ratio": float(in_stats["density"].mean() / max(out_stats["density"].mean(), 1)),
            "mannwhitney_p": float(mw_p),
            "te_positive_pct_set": float(te_pos_in / max(len(in_stats), 1) * 100),
            "te_positive_pct_bg": float(te_pos_out / max(len(out_stats), 1) * 100),
            "fisher_odds": float(odds) if not np.isinf(odds) else "Inf",
            "fisher_p": float(fisher_p),
            "mean_hit_count_set": float(in_stats["hit_count"].mean()),
            "mean_hit_count_bg": float(out_stats["hit_count"].mean()),
            **tier_breakdown,
        }

        print(f"\n  {set_name} ({len(in_stats)} genes):")
        print(f"    Mean density: {in_stats['density'].mean():.0f} vs {out_stats['density'].mean():.0f} (bg)")
        print(f"    Density ratio: {results[set_name]['density_ratio']:.2f}x")
        print(f"    MW p-value: {mw_p:.2e}")
        print(f"    TE+ rate: {results[set_name]['te_positive_pct_set']:.1f}% vs {results[set_name]['te_positive_pct_bg']:.1f}%")
        if tier_breakdown:
            print(f"    Strict/Moderate/Relaxed: "
                  f"{tier_breakdown.get('strict_hits', 0):,} / "
                  f"{tier_breakdown.get('moderate_hits', 0):,} / "
                  f"{tier_breakdown.get('relaxed_hits', 0):,}")

    with open(out_dir / "ws7_gene_set_enrichment.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

    elapsed = time.time() - t0
    print(f"\nWS7 complete in {elapsed:.1f}s")


# ── WS4: Unusual Matches ────────────────────────────────────────────────────

def run_ws4(df, out_dir):
    """WS4: Flag unusual/interesting BLAST hits."""
    print(f"\n{'='*60}")
    print(f"WS4: Unusual Matches")
    print(f"{'='*60}")
    t0 = time.time()

    # 4.1 Ultra-high identity hits (near-perfect matches = recent insertions?)
    perfect = df[(df["pident"] >= 95) & (df["length"] >= 200)]
    print(f"\n  Near-perfect hits (pident>=95, len>=200): {len(perfect):,}")
    if len(perfect) > 0:
        print(f"  Unique TEs in perfect hits: {perfect['sseqid'].nunique()}")
        print(f"  Top TEs:")
        for te, n in perfect["sseqid"].value_counts().head(5).items():
            print(f"    {te}: {n} hits")
        perfect_summary = perfect.groupby("sseqid").agg(
            n_hits=("qseqid", "count"),
            n_genes=("gene_id", "nunique"),
            mean_pident=("pident", "mean"),
            mean_length=("length", "mean"),
        ).sort_values("n_hits", ascending=False)
        perfect_summary.to_csv(out_dir / "ws4_near_perfect_hits.tsv", sep="\t", index=True)

    # 4.2 Extremely long alignments (> 500bp = large fossil fragments)
    long_hits = df[df["length"] >= 500]
    print(f"\n  Long alignments (>=500bp): {len(long_hits):,}")
    if len(long_hits) > 0:
        print(f"  Mean pident: {long_hits['pident'].mean():.1f}")
        print(f"  Max length: {long_hits['length'].max()}")
        long_summary = long_hits.groupby("sseqid").agg(
            n_hits=("qseqid", "count"),
            max_length=("length", "max"),
            mean_pident=("pident", "mean"),
        ).sort_values("max_length", ascending=False).head(20)
        long_summary.to_csv(out_dir / "ws4_long_alignments.tsv", sep="\t", index=True)

    # 4.3 Genes with extreme TE density (outliers)
    gene_stats_path = out_dir.parent / "gene_stats.tsv"
    if gene_stats_path.exists():
        gs = pd.read_csv(gene_stats_path, sep="\t", index_col=0)
        q99 = gs["density"].quantile(0.99)
        outliers = gs[gs["density"] > q99].sort_values("density", ascending=False)
        print(f"\n  Density outliers (>99th pct, density>{q99:.0f}): {len(outliers)}")
        if len(outliers) > 0:
            print(f"  Top 5:")
            for gene_id, row in outliers.head(5).iterrows():
                print(f"    {gene_id}: density={row['density']:.0f}, "
                      f"hits={row['hit_count']:.0f}, len={row['query_len']:.0f}")
            outliers.to_csv(out_dir / "ws4_density_outliers.tsv", sep="\t", index=True)

    # 4.4 Genes with NO TE hits (interesting in context of evalue=10)
    if gene_stats_path.exists():
        gs = pd.read_csv(gene_stats_path, sep="\t", index_col=0)
        no_hits = gs[gs["hit_count"] == 0]
        print(f"\n  Genes with ZERO hits (even at evalue=10): {len(no_hits)}")
        if len(no_hits) > 0:
            print(f"  Mean UTR length: {no_hits['query_len'].mean():.0f}bp")
            print(f"  Min UTR length: {no_hits['query_len'].min():.0f}bp")
            no_hits.to_csv(out_dir / "ws4_zero_hit_genes.tsv", sep="\t", index=True)

    # 4.5 Antisense-only genes (all hits antisense)
    if "gene_id" in df.columns:
        gene_strand = df.groupby("gene_id").agg(
            sense=("strand", lambda x: (x == "plus").sum()),
            antisense=("strand", lambda x: (x == "minus").sum()),
        )
        gene_strand["total"] = gene_strand["sense"] + gene_strand["antisense"]
        antisense_only = gene_strand[(gene_strand["sense"] == 0) & (gene_strand["total"] >= 10)]
        print(f"\n  Antisense-only genes (0 sense, >=10 hits): {len(antisense_only)}")
        if len(antisense_only) > 0:
            antisense_only.to_csv(out_dir / "ws4_antisense_only_genes.tsv", sep="\t", index=True)

    elapsed = time.time() - t0
    print(f"\nWS4 complete in {elapsed:.1f}s")


# ── Summary Report ───────────────────────────────────────────────────────────

def write_summary(utr_type, out_dir):
    """Compile all workstream outputs into a summary."""
    print(f"\n{'='*60}")
    print(f"Summary Report — {utr_type}")
    print(f"{'='*60}")

    summary = {"utr_type": utr_type, "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")}

    # Collect key files
    for fname in sorted(out_dir.glob("*.tsv")) + sorted(out_dir.glob("*.json")):
        if fname.name.startswith("blast_results"):
            continue
        if fname.suffix == ".tsv":
            try:
                n_rows = sum(1 for _ in open(fname)) - 1
                summary[f"file_{fname.stem}"] = f"{n_rows:,} rows"
            except Exception:
                pass
        elif fname.suffix == ".json":
            summary[f"file_{fname.stem}"] = "json"

    # Save analysis summary
    with open(out_dir / "analysis_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\n  Output files in {out_dir}:")
    for f in sorted(out_dir.glob("*")):
        if f.name.startswith("blast_results") or f.name == "regions.fa":
            continue
        size = f.stat().st_size
        if size > 1_000_000:
            print(f"    {f.name} ({size/1_000_000:.1f}MB)")
        elif size > 1000:
            print(f"    {f.name} ({size/1000:.0f}KB)")
        else:
            print(f"    {f.name} ({size}B)")


def run_ws6(df, regions_path, out_dir):
    """WS6: Conservation analysis — phyloP scoring of BLAST hits.

    Converts hits to genomic BED coords, scores with bigWigAverageOverBed,
    analyzes conservation by quality tier, known/novel status, and pident.
    """
    print(f"\n{'='*60}")
    print(f"WS6: Conservation (phyloP)")
    print(f"{'='*60}")
    t0 = time.time()

    PHYLOP_BW = DATA_DIR / "references" / "dm6.phyloP27way.bw"
    BIGWIG_TOOL = DATA_DIR / "references" / "tools" / "bigWigAverageOverBed"

    if not PHYLOP_BW.exists():
        print("  SKIP: phyloP bigWig not found")
        return
    if not BIGWIG_TOOL.exists():
        print("  SKIP: bigWigAverageOverBed not found")
        return

    # Build region lookup from regions.tsv
    if not regions_path.exists():
        print("  SKIP: regions.tsv not found")
        return
    regions = pd.read_csv(regions_path, sep="\t")
    region_lookup = {}
    for _, r in regions.iterrows():
        region_lookup[r["region_id"]] = {
            "chrom": r["chrom"],
            "start": int(r["start"]),
            "end": int(r["end"]),
            "strand": r["strand"],
        }

    # Convert BLAST hits to genomic BED coordinates (vectorized)
    # Join region metadata onto hits
    region_df = pd.DataFrame.from_dict(region_lookup, orient="index")
    region_df.index.name = "qseqid"
    region_df = region_df.rename(columns={
        "chrom": "r_chrom", "start": "r_start", "end": "r_end", "strand": "r_strand"
    })

    # Work on strict+moderate first (smaller, higher value), then sample relaxed
    tiers_to_score = df[df["tier"].isin(["strict", "moderate"])].copy()
    relaxed = df[df["tier"] == "relaxed"]
    n_relaxed_sample = min(200_000, len(relaxed))
    if n_relaxed_sample > 0:
        relaxed_sample = relaxed.sample(n=n_relaxed_sample, random_state=42)
        tiers_to_score = pd.concat([tiers_to_score, relaxed_sample], ignore_index=True)

    print(f"  Scoring {len(tiers_to_score):,} hits "
          f"(strict+moderate: {(tiers_to_score['tier'] != 'relaxed').sum():,}, "
          f"relaxed sample: {n_relaxed_sample:,})")

    # Merge region coords
    hits = tiers_to_score.merge(region_df, left_on="qseqid", right_index=True, how="left")
    hits = hits.dropna(subset=["r_chrom"])

    # Compute genomic BED coords
    # + strand: genomic_start = region_start + qstart - 1  (1-based)
    # - strand: genomic_start = region_end - qend + 1  (1-based)
    plus_mask = hits["r_strand"] == "+"
    hits.loc[plus_mask, "g_start"] = hits.loc[plus_mask, "r_start"] + hits.loc[plus_mask, "qstart"] - 2  # 0-based
    hits.loc[plus_mask, "g_end"] = hits.loc[plus_mask, "r_start"] + hits.loc[plus_mask, "qend"] - 1  # exclusive

    hits.loc[~plus_mask, "g_start"] = hits.loc[~plus_mask, "r_end"] - hits.loc[~plus_mask, "qend"]  # 0-based
    hits.loc[~plus_mask, "g_end"] = hits.loc[~plus_mask, "r_end"] - hits.loc[~plus_mask, "qstart"] + 1  # exclusive

    hits["g_start"] = hits["g_start"].astype(int)
    hits["g_end"] = hits["g_end"].astype(int)
    hits["bed_chrom"] = "chr" + hits["r_chrom"]

    # Assign unique hit IDs for merging back
    hits["hit_id"] = [f"hit_{i}" for i in range(len(hits))]

    # Write BED file (sorted by chrom, start for speed)
    bed_path = out_dir / "ws6_hits.bed"
    bed_df = hits[["bed_chrom", "g_start", "g_end", "hit_id"]].copy()
    bed_df = bed_df.sort_values(["bed_chrom", "g_start"])

    # Filter out any invalid coordinates
    bed_df = bed_df[(bed_df["g_start"] >= 0) & (bed_df["g_end"] > bed_df["g_start"])]
    bed_df.to_csv(bed_path, sep="\t", index=False, header=False)
    print(f"  Wrote {len(bed_df):,} BED intervals to {bed_path.name}")

    # Run bigWigAverageOverBed
    import subprocess
    phylop_out = out_dir / "ws6_phylop_raw.tab"
    cmd = [str(BIGWIG_TOOL), str(PHYLOP_BW), str(bed_path), str(phylop_out)]
    print(f"  Running bigWigAverageOverBed...")
    t1 = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ERROR: bigWigAverageOverBed failed: {result.stderr[:500]}")
        return
    print(f"  bigWigAverageOverBed completed in {time.time()-t1:.1f}s")

    # Parse results: columns are name, size, covered, sum, mean0, mean
    # mean0 = mean treating uncovered bases as 0
    # mean = mean of only covered bases
    phylop_scores = pd.read_csv(
        phylop_out, sep="\t", header=None,
        names=["hit_id", "size", "covered", "sum", "mean0", "mean"]
    )

    # Merge back
    hits = hits.merge(phylop_scores[["hit_id", "mean0", "mean", "covered", "size"]],
                      on="hit_id", how="left")
    hits["phylop_mean"] = hits["mean0"]  # mean0 treats uncovered as 0 (conservative)
    hits["phylop_coverage"] = hits["covered"] / hits["size"].clip(lower=1)

    # ── 6.1 PhyloP by quality tier ──
    print("\n  PhyloP by quality tier:")
    tier_stats = []
    for tier in ["strict", "moderate", "relaxed"]:
        subset = hits[hits["tier"] == tier]
        if len(subset) == 0:
            continue
        stats_row = {
            "tier": tier,
            "n_hits": len(subset),
            "mean_phylop": float(subset["phylop_mean"].mean()),
            "median_phylop": float(subset["phylop_mean"].median()),
            "std_phylop": float(subset["phylop_mean"].std()),
            "pct_positive": float((subset["phylop_mean"] > 0).mean() * 100),
            "pct_negative": float((subset["phylop_mean"] < 0).mean() * 100),
            "mean_coverage": float(subset["phylop_coverage"].mean()),
        }
        tier_stats.append(stats_row)
        print(f"    {tier}: mean={stats_row['mean_phylop']:.4f}, "
              f"median={stats_row['median_phylop']:.4f}, "
              f"positive={stats_row['pct_positive']:.1f}%, "
              f"coverage={stats_row['mean_coverage']:.1%}")
    pd.DataFrame(tier_stats).to_csv(out_dir / "phylop_by_quality.tsv", sep="\t", index=False)

    # ── 6.3 Known (RM) vs Novel conservation ──
    if "in_repeatmasker" in hits.columns:
        print("\n  Known vs Novel conservation:")
        known = hits[hits["in_repeatmasker"] == True]["phylop_mean"]
        novel = hits[hits["in_repeatmasker"] == False]["phylop_mean"]
        kn_stats = {}
        for label, data in [("known", known), ("novel", novel)]:
            if len(data) > 0:
                kn_stats[label] = {
                    "n": int(len(data)),
                    "mean_phylop": float(data.mean()),
                    "median_phylop": float(data.median()),
                    "std_phylop": float(data.std()),
                }
                print(f"    {label}: n={len(data):,}, mean={data.mean():.4f}, "
                      f"median={data.median():.4f}")
        if len(known) > 0 and len(novel) > 0:
            mw_stat, mw_p = stats.mannwhitneyu(known, novel, alternative="two-sided")
            kn_stats["mannwhitney_U"] = float(mw_stat)
            kn_stats["mannwhitney_p"] = float(mw_p)
            print(f"    Mann-Whitney U={mw_stat:.0f}, p={mw_p:.2e}")
        import json
        with open(out_dir / "known_vs_novel_conservation.json", "w") as f:
            json.dump(kn_stats, f, indent=2)

    # ── 6.5 Pident vs conservation correlation ──
    print("\n  Pident vs conservation correlation:")
    corr_results = {}
    valid = hits.dropna(subset=["phylop_mean", "pident"])
    if len(valid) >= 10:
        rho, p = stats.spearmanr(valid["pident"], valid["phylop_mean"])
        corr_results["overall"] = {"rho": float(rho), "p": float(p), "n": len(valid)}
        print(f"    Overall: Spearman rho={rho:.4f}, p={p:.2e}, n={len(valid):,}")

    for tier in ["strict", "moderate", "relaxed"]:
        subset = valid[valid["tier"] == tier]
        if len(subset) >= 10:
            rho, p = stats.spearmanr(subset["pident"], subset["phylop_mean"])
            corr_results[tier] = {"rho": float(rho), "p": float(p), "n": len(subset)}
            print(f"    {tier}: rho={rho:.4f}, p={p:.2e}, n={len(subset):,}")

    import json
    with open(out_dir / "quality_paradox_stats.json", "w") as f:
        json.dump(corr_results, f, indent=2)

    # Save per-hit phyloP data (for downstream WS2 motif analysis)
    phylop_cols = ["qseqid", "sseqid", "pident", "length", "tier", "hit_id",
                   "r_chrom", "g_start", "g_end", "phylop_mean", "phylop_coverage"]
    if "in_repeatmasker" in hits.columns:
        phylop_cols.append("in_repeatmasker")
    if "gene_id" in hits.columns:
        phylop_cols.append("gene_id")
    available_cols = [c for c in phylop_cols if c in hits.columns]
    hits[available_cols].to_csv(out_dir / "phylop_by_hit.tsv", sep="\t", index=False)

    # Cleanup temp files
    bed_path.unlink(missing_ok=True)
    phylop_out.unlink(missing_ok=True)

    print(f"\nWS6 complete in {time.time()-t0:.1f}s")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("TE Fossil Analysis — Full Workstream Run")
    print("CAVEAT: evalue=10, word_size=7, dust=no (very permissive)")
    print("=" * 60)

    # Load gene sets
    gene_sets_loaded = {}
    for name, path in GENE_SETS.items():
        if path.exists():
            gene_sets_loaded[name] = load_gene_set(path)
            print(f"  Gene set {name}: {len(gene_sets_loaded[name])} genes")
        else:
            print(f"  Gene set {name}: NOT FOUND at {path}")

    # Build lookup tables
    print("\nBuilding transcript→gene mapping from GFF...")
    t0 = time.time()
    t2g = build_transcript_to_gene_map()
    print(f"  {len(t2g):,} transcript→gene mappings in {time.time()-t0:.1f}s")

    print("Building TE consensus class map...")
    consensus_map = build_consensus_class_map()
    print(f"  {len(consensus_map):,} consensus families with class annotations")

    print("Building TE name mapping from FASTA headers...")
    te_names = build_te_name_map()
    print(f"  {len(te_names):,} TE ID→name mappings")

    # Save TE lookup table
    te_lookup_rows = []
    for te_id, name in te_names.items():
        te_lookup_rows.append({
            "te_id": te_id,
            "te_name": name,
            "te_class": infer_te_class_from_name(name, consensus_map),
        })
    if te_lookup_rows:
        te_lookup_df = pd.DataFrame(te_lookup_rows)
        te_lookup_df.to_csv(PROJECT_ROOT / "results" / "te_id_to_family_class.tsv",
                           sep="\t", index=False)
        class_counts = te_lookup_df["te_class"].value_counts()
        print(f"  TE class breakdown: {dict(class_counts)}")

    for utr_type, results_dir in RESULTS_DIRS.items():
        if not (results_dir / "blast_results.tsv").exists():
            print(f"\nSkipping {utr_type} — no blast_results.tsv")
            continue

        print(f"\n\n{'#'*60}")
        print(f"# Processing {utr_type.upper()}")
        print(f"{'#'*60}")

        out_dir = results_dir / "analysis"
        out_dir.mkdir(exist_ok=True)

        # Build query→gene mapping for this UTR type
        regions_path = results_dir / "regions.tsv"
        q2g = build_query_to_gene(regions_path, t2g) if regions_path.exists() else {}
        print(f"  Query→gene mapping: {len(q2g):,} regions mapped")

        # Load BLAST results
        print(f"\nLoading BLAST results from {results_dir / 'blast_results.tsv'}...")
        t0 = time.time()
        df = load_blast_results(results_dir / "blast_results.tsv")
        print(f"  Loaded {len(df):,} hits in {time.time()-t0:.1f}s")

        # Run workstreams
        df = run_ws1(df, utr_type, out_dir, q2g, te_names)
        df = run_ws3(df, gene_sets_loaded, out_dir, te_names, consensus_map)
        run_ws5(df, out_dir)
        run_ws7(df, gene_sets_loaded, out_dir)
        run_ws6(df, regions_path, out_dir)
        run_ws4(df, out_dir)
        write_summary(utr_type, out_dir)

    print(f"\n\n{'='*60}")
    print("ALL DONE")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
