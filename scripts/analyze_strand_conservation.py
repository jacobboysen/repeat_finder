#!/usr/bin/env python3
"""
Analyze correlation between strand bias and conservation/synteny.

Questions:
1. Are sense-biased UTRs more conserved than antisense-biased?
2. Are sense-biased TEs more conserved?
3. Any striking patterns by gene or TE family?
"""

import sys
from collections import defaultdict
import re

def load_strand_bias_utrs(filepath):
    """Load UTR-level strand bias."""
    data = {}
    with open(filepath) as f:
        header = next(f).strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            row = dict(zip(header, parts))
            transcript = row['fbtr']
            data[transcript] = {
                'fbgn': row['fbgn'],
                'total_hits': int(row['total_hits']),
                'sense_pct': float(row['sense_pct']),
                'anti_pct': float(row['anti_pct'])
            }
    return data

def load_strand_bias_tes(filepath):
    """Load TE-level strand bias."""
    data = {}
    with open(filepath) as f:
        header = next(f).strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            row = dict(zip(header, parts))
            te_id = row['te_id']
            data[te_id] = {
                'total_hits': int(row['total_hits']),
                'sense_pct': float(row['sense_pct']),
                'anti_pct': float(row['anti_pct']),
                'total_bp': int(row['total_bp'])
            }
    return data

def load_conservation_by_transcript(filepath):
    """Load conservation scores grouped by transcript."""
    data = defaultdict(list)
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split('\t')
            name = parts[0]
            phyloP = float(parts[5])

            name_parts = name.split('|')
            if len(name_parts) >= 2:
                transcript = name_parts[1]
                data[transcript].append(phyloP)
    return data

def load_conservation_by_te(filepath):
    """Load conservation scores grouped by TE."""
    data = defaultdict(list)
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split('\t')
            name = parts[0]
            phyloP = float(parts[5])

            name_parts = name.split('|')
            if len(name_parts) >= 3:
                te_id = name_parts[2]
                data[te_id].append(phyloP)
    return data

def load_synteny_by_transcript(filepath):
    """Load synteny scores grouped by transcript."""
    data = defaultdict(list)
    with open(filepath) as f:
        header = next(f)
        for line in f:
            parts = line.strip().split('\t')
            name = parts[3]
            any_cov = float(parts[11])

            name_parts = name.split('|')
            if len(name_parts) >= 2:
                transcript = name_parts[1]
                data[transcript].append(any_cov)
    return data

def load_synteny_by_te(filepath):
    """Load synteny scores grouped by TE."""
    data = defaultdict(list)
    with open(filepath) as f:
        header = next(f)
        for line in f:
            parts = line.strip().split('\t')
            name = parts[3]
            any_cov = float(parts[11])

            name_parts = name.split('|')
            if len(name_parts) >= 3:
                te_id = name_parts[2]
                data[te_id].append(any_cov)
    return data

def load_te_info(te_fasta):
    """Load TE family names."""
    te_info = {}
    with open(te_fasta) as f:
        for line in f:
            if line.startswith('>'):
                parts = line.split()
                te_id = parts[0][1:]
                name_match = re.search(r'name=([^;]+)', line)
                if name_match:
                    name = name_match.group(1)
                    family = re.sub(r'\{[^}]*\}.*', '', name)
                    te_info[te_id] = {'name': name, 'family': family}
    return te_info

def main():
    print("Loading data...", file=sys.stderr)

    # Load strand bias
    utr_strand = load_strand_bias_utrs('results/strand_bias_by_utr.tsv')
    te_strand = load_strand_bias_tes('results/strand_bias_by_te.tsv')

    # Load conservation
    utr_cons = load_conservation_by_transcript('results/repeatmasker_analysis/te_hits_all_conservation.tab')
    te_cons = load_conservation_by_te('results/repeatmasker_analysis/te_hits_all_conservation.tab')

    # Load synteny
    utr_syn = load_synteny_by_transcript('results/repeatmasker_analysis/te_hits_all_synteny_sampled.tsv')
    te_syn = load_synteny_by_te('results/repeatmasker_analysis/te_hits_all_synteny_sampled.tsv')

    # Load TE info
    te_info = load_te_info('/Users/jacobboysen/git_repos/repeat_finder/data/references/dmel_te_flybase.fasta')

    print(f"  UTR strand bias: {len(utr_strand)}", file=sys.stderr)
    print(f"  TE strand bias: {len(te_strand)}", file=sys.stderr)
    print(f"  UTR conservation: {len(utr_cons)}", file=sys.stderr)
    print(f"  TE conservation: {len(te_cons)}", file=sys.stderr)

    # =====================================================================
    # UTR-LEVEL ANALYSIS: Strand bias vs Conservation
    # =====================================================================
    print("\n" + "="*70)
    print("UTR-LEVEL: STRAND BIAS vs CONSERVATION")
    print("="*70)

    # Categorize UTRs by strand bias
    sense_biased = []  # >70% sense
    antisense_biased = []  # >70% antisense
    balanced = []  # 30-70% sense

    for transcript, strand_data in utr_strand.items():
        if strand_data['total_hits'] < 10:  # Require minimum hits
            continue

        if transcript in utr_cons and len(utr_cons[transcript]) > 0:
            mean_cons = sum(utr_cons[transcript]) / len(utr_cons[transcript])

            entry = {
                'transcript': transcript,
                'sense_pct': strand_data['sense_pct'],
                'mean_cons': mean_cons,
                'n_hits': strand_data['total_hits']
            }

            if strand_data['sense_pct'] >= 70:
                sense_biased.append(entry)
            elif strand_data['sense_pct'] <= 30:
                antisense_biased.append(entry)
            else:
                balanced.append(entry)

    print(f"\nUTRs with ≥10 hits and conservation data:")
    print(f"  Sense-biased (≥70%): {len(sense_biased)}")
    print(f"  Antisense-biased (≤30%): {len(antisense_biased)}")
    print(f"  Balanced (30-70%): {len(balanced)}")

    if sense_biased and antisense_biased and balanced:
        sense_mean = sum(e['mean_cons'] for e in sense_biased) / len(sense_biased)
        anti_mean = sum(e['mean_cons'] for e in antisense_biased) / len(antisense_biased)
        bal_mean = sum(e['mean_cons'] for e in balanced) / len(balanced)

        print(f"\nMean conservation by strand bias:")
        print(f"  Sense-biased UTRs:     {sense_mean:.3f}")
        print(f"  Antisense-biased UTRs: {anti_mean:.3f}")
        print(f"  Balanced UTRs:         {bal_mean:.3f}")

        # Statistical significance hint
        if sense_mean > anti_mean:
            diff_pct = 100 * (sense_mean - anti_mean) / anti_mean
            print(f"\n  >> Sense-biased UTRs are {diff_pct:.1f}% MORE conserved!")
        else:
            diff_pct = 100 * (anti_mean - sense_mean) / sense_mean
            print(f"\n  >> Antisense-biased UTRs are {diff_pct:.1f}% MORE conserved!")

    # =====================================================================
    # UTR-LEVEL: Strand bias vs Synteny
    # =====================================================================
    print("\n" + "-"*50)
    print("UTR-LEVEL: STRAND BIAS vs SYNTENY")
    print("-"*50)

    sense_syn = []
    anti_syn = []
    bal_syn = []

    for transcript, strand_data in utr_strand.items():
        if strand_data['total_hits'] < 10:
            continue

        if transcript in utr_syn and len(utr_syn[transcript]) > 0:
            mean_syn = sum(utr_syn[transcript]) / len(utr_syn[transcript])
            pct_syntenic = 100 * sum(1 for s in utr_syn[transcript] if s >= 0.5) / len(utr_syn[transcript])

            entry = {
                'transcript': transcript,
                'sense_pct': strand_data['sense_pct'],
                'mean_syn': mean_syn,
                'pct_syntenic': pct_syntenic
            }

            if strand_data['sense_pct'] >= 70:
                sense_syn.append(entry)
            elif strand_data['sense_pct'] <= 30:
                anti_syn.append(entry)
            else:
                bal_syn.append(entry)

    if sense_syn and anti_syn:
        sense_pct_syn = sum(e['pct_syntenic'] for e in sense_syn) / len(sense_syn)
        anti_pct_syn = sum(e['pct_syntenic'] for e in anti_syn) / len(anti_syn)
        bal_pct_syn = sum(e['pct_syntenic'] for e in bal_syn) / len(bal_syn) if bal_syn else 0

        print(f"\nMean % hits syntenic by strand bias:")
        print(f"  Sense-biased UTRs:     {sense_pct_syn:.1f}%")
        print(f"  Antisense-biased UTRs: {anti_pct_syn:.1f}%")
        print(f"  Balanced UTRs:         {bal_pct_syn:.1f}%")

    # =====================================================================
    # TE-LEVEL ANALYSIS: Strand bias vs Conservation
    # =====================================================================
    print("\n" + "="*70)
    print("TE-LEVEL: STRAND BIAS vs CONSERVATION")
    print("="*70)

    te_analysis = []
    for te_id, strand_data in te_strand.items():
        if strand_data['total_hits'] < 100:  # Require more hits for TEs
            continue

        if te_id in te_cons and len(te_cons[te_id]) > 0:
            mean_cons = sum(te_cons[te_id]) / len(te_cons[te_id])
            info = te_info.get(te_id, {'name': te_id, 'family': te_id})

            te_analysis.append({
                'te_id': te_id,
                'name': info['name'],
                'family': info['family'],
                'sense_pct': strand_data['sense_pct'],
                'mean_cons': mean_cons,
                'n_hits': strand_data['total_hits']
            })

    # Sort by sense bias
    te_analysis.sort(key=lambda x: -x['sense_pct'])

    print(f"\nTEs with ≥100 hits:")
    print(f"{'TE Name':<25} {'Sense%':>8} {'Cons':>8} {'Hits':>10}")
    print("-" * 55)

    # Top 10 most sense-biased
    print("\nMost SENSE-biased TEs:")
    for entry in te_analysis[:10]:
        print(f"  {entry['name'][:24]:<25} {entry['sense_pct']:>7.1f}% {entry['mean_cons']:>8.2f} {entry['n_hits']:>10,}")

    # Top 10 most antisense-biased
    print("\nMost ANTISENSE-biased TEs:")
    for entry in te_analysis[-10:]:
        print(f"  {entry['name'][:24]:<25} {entry['sense_pct']:>7.1f}% {entry['mean_cons']:>8.2f} {entry['n_hits']:>10,}")

    # Compare conservation
    sense_tes = [e for e in te_analysis if e['sense_pct'] >= 60]
    anti_tes = [e for e in te_analysis if e['sense_pct'] <= 40]

    if sense_tes and anti_tes:
        sense_te_cons = sum(e['mean_cons'] for e in sense_tes) / len(sense_tes)
        anti_te_cons = sum(e['mean_cons'] for e in anti_tes) / len(anti_tes)

        print(f"\nMean conservation by TE strand bias:")
        print(f"  Sense-biased TEs (≥60%):     {sense_te_cons:.3f} (n={len(sense_tes)})")
        print(f"  Antisense-biased TEs (≤40%): {anti_te_cons:.3f} (n={len(anti_tes)})")

    # =====================================================================
    # TE FAMILY-LEVEL ANALYSIS
    # =====================================================================
    print("\n" + "="*70)
    print("TE FAMILY-LEVEL: STRAND BIAS vs CONSERVATION")
    print("="*70)

    family_data = defaultdict(lambda: {'sense_pcts': [], 'conservations': [], 'n_hits': 0})

    for entry in te_analysis:
        family = entry['family']
        family_data[family]['sense_pcts'].append(entry['sense_pct'])
        family_data[family]['conservations'].append(entry['mean_cons'])
        family_data[family]['n_hits'] += entry['n_hits']

    # Summarize by family
    family_summary = []
    for family, data in family_data.items():
        if len(data['sense_pcts']) >= 3:  # Require multiple TE instances
            family_summary.append({
                'family': family,
                'mean_sense': sum(data['sense_pcts']) / len(data['sense_pcts']),
                'mean_cons': sum(data['conservations']) / len(data['conservations']),
                'n_instances': len(data['sense_pcts']),
                'n_hits': data['n_hits']
            })

    family_summary.sort(key=lambda x: -x['mean_sense'])

    print(f"\n{'TE Family':<20} {'Mean Sense%':>12} {'Mean Cons':>10} {'Instances':>10}")
    print("-" * 55)

    for entry in family_summary[:20]:
        print(f"{entry['family']:<20} {entry['mean_sense']:>11.1f}% {entry['mean_cons']:>10.2f} {entry['n_instances']:>10}")

    # Correlation analysis
    print("\n" + "="*70)
    print("CORRELATION: STRAND BIAS vs CONSERVATION")
    print("="*70)

    # Simple correlation for TEs
    if len(te_analysis) > 10:
        sense_pcts = [e['sense_pct'] for e in te_analysis]
        conservations = [e['mean_cons'] for e in te_analysis]

        # Calculate Pearson correlation manually
        n = len(sense_pcts)
        mean_x = sum(sense_pcts) / n
        mean_y = sum(conservations) / n

        num = sum((sense_pcts[i] - mean_x) * (conservations[i] - mean_y) for i in range(n))
        den_x = sum((x - mean_x)**2 for x in sense_pcts) ** 0.5
        den_y = sum((y - mean_y)**2 for y in conservations) ** 0.5

        if den_x > 0 and den_y > 0:
            correlation = num / (den_x * den_y)
            print(f"\nTE-level correlation (sense% vs conservation): r = {correlation:.3f}")

            if correlation > 0.1:
                print("  >> Positive correlation: Sense-biased TEs tend to be MORE conserved")
            elif correlation < -0.1:
                print("  >> Negative correlation: Antisense-biased TEs tend to be MORE conserved")
            else:
                print("  >> Weak/no correlation between strand bias and conservation")

if __name__ == '__main__':
    main()
