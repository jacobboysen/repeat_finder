"""K-mer motif enrichment analysis for TE fossil mining.

Provides functions to extract k-mers from sequences, compare frequencies
between real TE hits and shuffled controls, compute positional profiles,
and test motif enrichment in gene sets.
"""

import numpy as np
import pandas as pd
from scipy import stats


def _benjamini_hochberg(p_values):
    """Benjamini-Hochberg FDR correction.

    Args:
        p_values: Array-like of raw p-values.

    Returns:
        numpy array of FDR-corrected q-values.
    """
    n = len(p_values)
    if n == 0:
        return np.array([])
    sorted_idx = np.argsort(p_values)
    sorted_p = np.array(p_values)[sorted_idx]
    fdr = np.minimum(1.0, sorted_p * n / (np.arange(1, n + 1)))
    # Enforce monotonicity (from end)
    for i in range(n - 2, -1, -1):
        fdr[i] = min(fdr[i], fdr[i + 1])
    result = np.empty(n)
    result[sorted_idx] = fdr
    return result


_VALID_BASES = set("ACGT")


def extract_kmers(sequence: str, k: int = 6) -> list[str]:
    """Extract all k-mers from a sequence, removing gaps first.

    Gaps ('-') are removed before extraction. K-mers containing
    non-ACGT characters are skipped.

    Args:
        sequence: Nucleotide sequence (may contain '-' gaps).
        k: K-mer length (default: 6).

    Returns:
        List of valid k-mers from the ungapped sequence.
    """
    clean_seq = sequence.replace("-", "").upper()

    if len(clean_seq) < k:
        return []

    kmers = []
    for i in range(len(clean_seq) - k + 1):
        kmer = clean_seq[i : i + k]
        if all(base in _VALID_BASES for base in kmer):
            kmers.append(kmer)
    return kmers


def count_kmers_in_sequences(
    sequences: list[str], k: int = 6
) -> dict[str, int]:
    """Count k-mer frequencies across a list of sequences.

    Args:
        sequences: List of nucleotide sequences.
        k: K-mer length (default: 6).

    Returns:
        Dict mapping each k-mer to its total count.
    """
    counts: dict[str, int] = {}
    for seq in sequences:
        for kmer in extract_kmers(seq, k):
            counts[kmer] = counts.get(kmer, 0) + 1
    return counts


def count_kmers_from_blast(
    df: pd.DataFrame, k: int = 6, seq_col: str = "qseq"
) -> dict[str, int]:
    """Count k-mers from BLAST results using a specified sequence column.

    Args:
        df: BLAST results DataFrame.
        k: K-mer length (default: 6).
        seq_col: Column containing aligned sequences (default: "qseq").

    Returns:
        Dict mapping each k-mer to its total count.
    """
    if df.empty or seq_col not in df.columns:
        return {}

    counts: dict[str, int] = {}
    for seq in df[seq_col].dropna():
        for kmer in extract_kmers(str(seq), k):
            counts[kmer] = counts.get(kmer, 0) + 1
    return counts


def compute_kmer_enrichment(
    real_counts: dict[str, int],
    shuffled_counts_list: list[dict[str, int]],
    min_count: int = 10,
) -> pd.DataFrame:
    """Compare real k-mer counts against shuffled replicates.

    Computes enrichment ratio, z-score, p-value, and FDR for each k-mer.

    Args:
        real_counts: K-mer counts from real data.
        shuffled_counts_list: List of k-mer count dicts from shuffled
            replicates.
        min_count: Minimum total count (real + all shuffled) to include
            a k-mer (default: 10).

    Returns:
        DataFrame sorted by p_value with columns: motif, real_count,
        shuf_mean, shuf_std, enrichment, log2_enrichment, z_score,
        p_value, fdr.
    """
    if not real_counts and not shuffled_counts_list:
        return pd.DataFrame(
            columns=[
                "motif", "real_count", "shuf_mean", "shuf_std",
                "enrichment", "log2_enrichment", "z_score", "p_value", "fdr",
            ]
        )

    # Gather all k-mers
    all_kmers: set[str] = set(real_counts.keys())
    for sc in shuffled_counts_list:
        all_kmers.update(sc.keys())

    n_reps = len(shuffled_counts_list)
    pseudocount = 1.0
    results = []

    for kmer in all_kmers:
        real_count = real_counts.get(kmer, 0)
        shuf_counts = [sc.get(kmer, 0) for sc in shuffled_counts_list]
        shuf_mean = float(np.mean(shuf_counts)) if shuf_counts else 0.0
        shuf_std = (
            float(np.std(shuf_counts, ddof=1)) if n_reps > 1 else 0.0
        )

        total_count = real_count + sum(shuf_counts)
        if total_count < min_count:
            continue

        enrichment = (real_count + pseudocount) / (shuf_mean + pseudocount)
        log2_enrichment = float(np.log2(enrichment))

        if shuf_std > 0:
            z_score = (real_count - shuf_mean) / shuf_std
        else:
            z_score = (
                float(np.sign(real_count - shuf_mean)) * 10.0
                if real_count != shuf_mean
                else 0.0
            )

        p_value = float(2.0 * (1.0 - stats.norm.cdf(abs(z_score))))

        results.append(
            {
                "motif": kmer,
                "real_count": real_count,
                "shuf_mean": shuf_mean,
                "shuf_std": shuf_std,
                "enrichment": enrichment,
                "log2_enrichment": log2_enrichment,
                "z_score": z_score,
                "p_value": p_value,
            }
        )

    if not results:
        return pd.DataFrame(
            columns=[
                "motif", "real_count", "shuf_mean", "shuf_std",
                "enrichment", "log2_enrichment", "z_score", "p_value", "fdr",
            ]
        )

    df = pd.DataFrame(results)
    df["fdr"] = _benjamini_hochberg(df["p_value"].values)
    df = df.sort_values("p_value").reset_index(drop=True)
    return df


def find_motif_positions(sequence: str, motif: str) -> list[int]:
    """Find all start positions (0-indexed) of a motif in a gap-cleaned sequence.

    Case-insensitive matching.

    Args:
        sequence: Nucleotide sequence (may contain '-' gaps).
        motif: Motif string to search for.

    Returns:
        List of 0-indexed start positions in the ungapped sequence.
    """
    clean_seq = sequence.replace("-", "").upper()
    motif_upper = motif.upper()

    positions = []
    start = 0
    while True:
        pos = clean_seq.find(motif_upper, start)
        if pos == -1:
            break
        positions.append(pos)
        start = pos + 1
    return positions


def compute_motif_positional_profile(
    df: pd.DataFrame, motifs: list[str], n_bins: int = 10
) -> pd.DataFrame:
    """Compute positional distribution of motifs within query sequences.

    For each motif, positions are normalized to [0, 1] relative to query
    length, then binned.

    Args:
        df: BLAST results DataFrame with qseq, qstart, qlen columns.
        motifs: List of motifs to profile.
        n_bins: Number of position bins (default: 10).

    Returns:
        DataFrame with columns: motif, bin_start, bin_end, count, density.
    """
    required_cols = {"qseq", "qstart", "qlen"}
    if df.empty or not required_cols.issubset(df.columns):
        rows = []
        for motif in motifs:
            for i in range(n_bins):
                rows.append(
                    {
                        "motif": motif,
                        "bin_start": i / n_bins,
                        "bin_end": (i + 1) / n_bins,
                        "count": 0,
                        "density": 0.0,
                    }
                )
        return pd.DataFrame(rows)

    bin_edges = np.linspace(0, 1, n_bins + 1)
    all_rows = []

    for motif in motifs:
        rel_positions: list[float] = []

        for _, row in df.iterrows():
            qseq = str(row["qseq"])
            qstart = float(row["qstart"])
            qlen = float(row["qlen"])
            if qlen == 0:
                continue

            positions = find_motif_positions(qseq, motif)
            for pos in positions:
                rel_pos = (qstart + pos) / qlen
                rel_pos = max(0.0, min(1.0, rel_pos))
                rel_positions.append(rel_pos)

        if rel_positions:
            counts, _ = np.histogram(rel_positions, bins=bin_edges)
        else:
            counts = np.zeros(n_bins, dtype=int)

        total = int(counts.sum())

        for i in range(n_bins):
            all_rows.append(
                {
                    "motif": motif,
                    "bin_start": bin_edges[i],
                    "bin_end": bin_edges[i + 1],
                    "count": int(counts[i]),
                    "density": float(counts[i]) / total if total > 0 else 0.0,
                }
            )

    return pd.DataFrame(all_rows)


def compute_motif_gene_set_enrichment(
    df: pd.DataFrame,
    motifs: list[str],
    query_to_gene: dict[str, str],
    gene_sets: dict[str, set[str]],
    seq_col: str = "qseq",
) -> pd.DataFrame:
    """Test whether motifs are enriched in specific gene sets.

    For each motif x gene_set combination, runs Fisher's exact test on the
    2x2 table: (genes_with_motif_in_set, genes_without_motif_in_set,
    genes_with_motif_not_in_set, genes_without_motif_not_in_set).

    Args:
        df: BLAST results DataFrame.
        motifs: List of motifs to test.
        query_to_gene: Mapping from query ID (qseqid) to gene ID.
        gene_sets: Named gene sets, mapping set name to gene IDs.
        seq_col: Column containing aligned sequences (default: "qseq").

    Returns:
        DataFrame with columns: motif, gene_set, odds_ratio, p_value,
        n_in_set, n_with_motif.
    """
    if df.empty or seq_col not in df.columns or "qseqid" not in df.columns:
        rows = []
        for motif in motifs:
            for gs_name in gene_sets:
                rows.append(
                    {
                        "motif": motif,
                        "gene_set": gs_name,
                        "odds_ratio": float("nan"),
                        "p_value": float("nan"),
                        "n_in_set": 0,
                        "n_with_motif": 0,
                    }
                )
        return pd.DataFrame(rows)

    # Map each query to its gene
    all_genes: set[str] = set()
    gene_has_motif: dict[str, set[str]] = {m: set() for m in motifs}

    for _, row in df.iterrows():
        qseqid = str(row["qseqid"])
        gene = query_to_gene.get(qseqid)
        if gene is None:
            continue
        all_genes.add(gene)

        seq = str(row[seq_col])
        for motif in motifs:
            if find_motif_positions(seq, motif):
                gene_has_motif[motif].add(gene)

    results = []
    for motif in motifs:
        genes_with = gene_has_motif[motif]
        genes_without = all_genes - genes_with

        for gs_name, gs_genes in gene_sets.items():
            a = len(genes_with & gs_genes)      # motif + in set
            b = len(genes_without & gs_genes)    # no motif + in set
            c = len(genes_with - gs_genes)       # motif + not in set
            d = len(genes_without - gs_genes)    # no motif + not in set

            if a + b + c + d == 0:
                odds_ratio = float("nan")
                p_value = float("nan")
            else:
                odds_ratio, p_value = stats.fisher_exact(
                    [[a, b], [c, d]], alternative="two-sided"
                )

            results.append(
                {
                    "motif": motif,
                    "gene_set": gs_name,
                    "odds_ratio": odds_ratio,
                    "p_value": p_value,
                    "n_in_set": a + b,
                    "n_with_motif": a + c,
                }
            )

    return pd.DataFrame(results)
