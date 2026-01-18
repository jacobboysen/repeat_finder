#!/usr/bin/env python3
"""
Analyze shuffle convergence for duplicated sequences.

For each duplicated sequence, examine:
- How many shuffles were performed (duplication × 10 reps)
- Hit count spread across shuffles
- Position convergence (do shuffles hit same TE positions?)
- Identity and length distributions
"""

from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Set up matplotlib
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10

print("=" * 70)
print("SHUFFLE CONVERGENCE ANALYSIS FOR DUPLICATED SEQUENCES")
print("=" * 70)

# Load sequences and find duplicates
print("\nLoading sequences...")
seqs = {}
current_id = None
current_seq = []
with open('data/references/dmel_3utr.fasta') as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            if current_id:
                seqs[current_id] = ''.join(current_seq)
            current_id = line[1:].split()[0]
            current_seq = []
        else:
            current_seq.append(line)
    if current_id:
        seqs[current_id] = ''.join(current_seq)

# Group by content
by_content = defaultdict(list)
for seq_id, seq in seqs.items():
    by_content[seq].append(seq_id)

# Filter to duplicated sequences (2+ copies)
duplicates = {seq: ids for seq, ids in by_content.items() if len(ids) >= 2}
print(f"Total sequences: {len(seqs):,}")
print(f"Unique sequences: {len(by_content):,}")
print(f"Duplicated sequences (2+ copies): {len(duplicates):,}")

# Analyze each duplicated sequence
print("\nAnalyzing shuffle hits for each duplicated sequence...")

results = []

for seq_idx, (seq, transcript_ids) in enumerate(duplicates.items()):
    n_copies = len(transcript_ids)
    n_shuffles = n_copies * 10  # 10 replicates
    seq_len = len(seq)

    # Collect hits for all shuffled versions
    hits_per_version = defaultdict(list)

    for rep in range(1, 11):
        blast_file = f'results/shuffled_full/replicate_{rep:02d}_blast.tsv'
        with open(blast_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                qseqid = parts[0]
                base_id = qseqid.rsplit('_shuf', 1)[0]
                if base_id in transcript_ids:
                    te_id = parts[1]
                    pident = float(parts[2])
                    length = int(parts[3])
                    sstart = int(parts[8])
                    send = int(parts[9])
                    hits_per_version[qseqid].append({
                        'te': te_id,
                        'pos': (te_id, min(sstart, send), max(sstart, send)),
                        'pident': pident,
                        'length': length
                    })

    # Skip if no hits
    if not hits_per_version:
        continue

    # Calculate metrics
    hit_counts = [len(v) for v in hits_per_version.values()]
    all_hits = [h for v in hits_per_version.values() for h in v]

    # Position convergence
    all_positions = set(h['pos'] for h in all_hits)
    position_counts = Counter(h['pos'] for h in all_hits)

    # Identity and length
    pidents = [h['pident'] for h in all_hits]
    lengths = [h['length'] for h in all_hits]

    # Unique TEs
    unique_tes = set(h['te'] for h in all_hits)

    results.append({
        'seq_len': seq_len,
        'n_copies': n_copies,
        'n_shuffles': n_shuffles,
        'versions_with_hits': len(hits_per_version),
        'total_hits': len(all_hits),
        'unique_positions': len(all_positions),
        'unique_tes': len(unique_tes),
        'hits_per_position': len(all_hits) / len(all_positions) if all_positions else 0,
        'hit_count_min': min(hit_counts),
        'hit_count_max': max(hit_counts),
        'hit_count_mean': np.mean(hit_counts),
        'hit_count_std': np.std(hit_counts),
        'pident_mean': np.mean(pidents),
        'pident_std': np.std(pidents),
        'pident_min': min(pidents),
        'pident_max': max(pidents),
        'length_mean': np.mean(lengths),
        'length_std': np.std(lengths),
        'length_min': min(lengths),
        'length_max': max(lengths),
        'max_position_hits': max(position_counts.values()),
    })

    if (seq_idx + 1) % 100 == 0:
        print(f"  Processed {seq_idx + 1}/{len(duplicates)} duplicated sequences...")

print(f"\nAnalyzed {len(results)} duplicated sequences with hits")

# Summary statistics
print("\n" + "=" * 70)
print("SUMMARY STATISTICS")
print("=" * 70)

# Convert to arrays for easier analysis
n_copies_arr = np.array([r['n_copies'] for r in results])
n_shuffles_arr = np.array([r['n_shuffles'] for r in results])
total_hits_arr = np.array([r['total_hits'] for r in results])
unique_pos_arr = np.array([r['unique_positions'] for r in results])
hits_per_pos_arr = np.array([r['hits_per_position'] for r in results])
pident_std_arr = np.array([r['pident_std'] for r in results])
length_std_arr = np.array([r['length_std'] for r in results])

print(f"\nPosition Convergence (hits per unique position):")
print(f"  Mean: {np.mean(hits_per_pos_arr):.2f}")
print(f"  Median: {np.median(hits_per_pos_arr):.2f}")
print(f"  Max: {np.max(hits_per_pos_arr):.2f}")
print(f"  % with convergence ratio < 1.1: {100*np.mean(hits_per_pos_arr < 1.1):.1f}%")
print(f"  % with convergence ratio < 1.5: {100*np.mean(hits_per_pos_arr < 1.5):.1f}%")

print(f"\nIdentity spread (std dev across shuffles):")
print(f"  Mean std: {np.mean(pident_std_arr):.2f}%")
print(f"  Max std: {np.max(pident_std_arr):.2f}%")

print(f"\nLength spread (std dev across shuffles):")
print(f"  Mean std: {np.mean(length_std_arr):.1f} bp")
print(f"  Max std: {np.max(length_std_arr):.1f} bp")

# Create figures
print("\nGenerating figures...")
fig_dir = Path('figures/shuffle_convergence')
fig_dir.mkdir(parents=True, exist_ok=True)

fig, axes = plt.subplots(2, 3, figsize=(15, 10))

# Panel 1: Convergence ratio distribution
ax = axes[0, 0]
ax.hist(hits_per_pos_arr, bins=50, color='steelblue', edgecolor='white', alpha=0.8)
ax.axvline(1.0, color='red', linestyle='--', linewidth=2, label='Perfect divergence (1.0)')
ax.axvline(np.median(hits_per_pos_arr), color='orange', linestyle='-', linewidth=2,
           label=f'Median: {np.median(hits_per_pos_arr):.2f}')
ax.set_xlabel('Hits per Unique TE Position')
ax.set_ylabel('Number of Duplicated Sequences')
ax.set_title('A. Position Convergence\n(1.0 = every hit at new position)')
ax.legend(fontsize=9)
ax.set_xlim(0.9, min(5, np.percentile(hits_per_pos_arr, 99)))

# Panel 2: Convergence vs number of shuffles
ax = axes[0, 1]
ax.scatter(n_shuffles_arr, hits_per_pos_arr, alpha=0.5, s=20, c='steelblue')
ax.axhline(1.0, color='red', linestyle='--', linewidth=1)
ax.set_xlabel('Number of Shuffles (copies × 10 reps)')
ax.set_ylabel('Hits per Unique Position')
ax.set_title('B. Convergence vs Shuffle Count')
ax.set_ylim(0.9, min(5, np.percentile(hits_per_pos_arr, 99)))

# Panel 3: Identity spread distribution
ax = axes[0, 2]
ax.hist(pident_std_arr, bins=50, color='coral', edgecolor='white', alpha=0.8)
ax.axvline(np.mean(pident_std_arr), color='red', linestyle='-', linewidth=2,
           label=f'Mean: {np.mean(pident_std_arr):.1f}%')
ax.set_xlabel('Identity Std Dev (%)')
ax.set_ylabel('Number of Duplicated Sequences')
ax.set_title('C. Identity Spread Across Shuffles')
ax.legend(fontsize=9)

# Panel 4: Length spread distribution
ax = axes[1, 0]
ax.hist(length_std_arr, bins=50, color='seagreen', edgecolor='white', alpha=0.8)
ax.axvline(np.mean(length_std_arr), color='red', linestyle='-', linewidth=2,
           label=f'Mean: {np.mean(length_std_arr):.1f}bp')
ax.set_xlabel('Length Std Dev (bp)')
ax.set_ylabel('Number of Duplicated Sequences')
ax.set_title('D. Length Spread Across Shuffles')
ax.legend(fontsize=9)

# Panel 5: Unique positions vs shuffles
ax = axes[1, 1]
ax.scatter(n_shuffles_arr, unique_pos_arr, alpha=0.5, s=20, c='purple')
z = np.polyfit(n_shuffles_arr, unique_pos_arr, 1)
p = np.poly1d(z)
x_line = np.linspace(n_shuffles_arr.min(), n_shuffles_arr.max(), 100)
ax.plot(x_line, p(x_line), 'r--', linewidth=2, label=f'Linear fit')
ax.set_xlabel('Number of Shuffles')
ax.set_ylabel('Unique TE Positions Hit')
ax.set_title('E. Position Discovery vs Shuffle Count')
ax.legend(fontsize=9)

# Panel 6: Summary text
ax = axes[1, 2]
ax.axis('off')
summary = f"""
CONVERGENCE ANALYSIS SUMMARY

Duplicated sequences analyzed: {len(results):,}
Total shuffles across all: {sum(n_shuffles_arr):,}

POSITION CONVERGENCE:
  Hits per unique position:
    Mean: {np.mean(hits_per_pos_arr):.2f}
    Median: {np.median(hits_per_pos_arr):.2f}

  {100*np.mean(hits_per_pos_arr < 1.1):.0f}% of sequences show <1.1 ratio
  (nearly every hit at new position)

IDENTITY VARIATION:
  Mean std across shuffles: {np.mean(pident_std_arr):.1f}%
  Range: {np.min(pident_std_arr):.1f} - {np.max(pident_std_arr):.1f}%

LENGTH VARIATION:
  Mean std across shuffles: {np.mean(length_std_arr):.0f} bp
  Range: {np.min(length_std_arr):.0f} - {np.max(length_std_arr):.0f} bp

CONCLUSION:
Shuffles show minimal convergence - each shuffle
finds mostly NEW TE positions, not the same ones.
More shuffles = more positions discovered.
"""
ax.text(0.05, 0.95, summary, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

plt.tight_layout()
plt.savefig(fig_dir / 'shuffle_convergence_analysis.png', dpi=150, bbox_inches='tight',
            facecolor='white')
print(f"  Saved: {fig_dir}/shuffle_convergence_analysis.png")

# Second figure: detailed per-copy analysis
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Group by copy number
copy_groups = defaultdict(list)
for r in results:
    copy_groups[r['n_copies']].append(r)

copy_nums = sorted(copy_groups.keys())

# Panel 1: Convergence by copy number
ax = axes[0, 0]
convergence_by_copy = [np.mean([r['hits_per_position'] for r in copy_groups[c]]) for c in copy_nums]
convergence_std = [np.std([r['hits_per_position'] for r in copy_groups[c]]) for c in copy_nums]
ax.errorbar(copy_nums, convergence_by_copy, yerr=convergence_std, fmt='o-', capsize=3, color='steelblue')
ax.axhline(1.0, color='red', linestyle='--', alpha=0.7)
ax.set_xlabel('Number of Copies (duplication level)')
ax.set_ylabel('Mean Hits per Unique Position')
ax.set_title('A. Convergence by Duplication Level')
ax.set_xlim(0, min(80, max(copy_nums)+5))

# Panel 2: Total positions discovered by copy number
ax = axes[0, 1]
positions_by_copy = [np.mean([r['unique_positions'] for r in copy_groups[c]]) for c in copy_nums]
ax.scatter(copy_nums, positions_by_copy, s=50, c='purple', alpha=0.7)
ax.set_xlabel('Number of Copies')
ax.set_ylabel('Mean Unique Positions Discovered')
ax.set_title('B. Position Discovery by Duplication Level')
ax.set_xlim(0, min(80, max(copy_nums)+5))

# Panel 3: Identity spread by copy number
ax = axes[1, 0]
pident_by_copy = [np.mean([r['pident_std'] for r in copy_groups[c]]) for c in copy_nums]
ax.scatter(copy_nums, pident_by_copy, s=50, c='coral', alpha=0.7)
ax.set_xlabel('Number of Copies')
ax.set_ylabel('Mean Identity Std Dev (%)')
ax.set_title('C. Identity Variation by Duplication Level')
ax.set_xlim(0, min(80, max(copy_nums)+5))

# Panel 4: Sample size per copy number
ax = axes[1, 1]
counts_by_copy = [len(copy_groups[c]) for c in copy_nums]
ax.bar(copy_nums, counts_by_copy, color='gray', edgecolor='black', alpha=0.7)
ax.set_xlabel('Number of Copies')
ax.set_ylabel('Number of Sequences')
ax.set_title('D. Sample Size by Duplication Level')
ax.set_xlim(0, min(80, max(copy_nums)+5))

plt.tight_layout()
plt.savefig(fig_dir / 'shuffle_convergence_by_copies.png', dpi=150, bbox_inches='tight',
            facecolor='white')
print(f"  Saved: {fig_dir}/shuffle_convergence_by_copies.png")

# Save detailed results
results_file = fig_dir / 'convergence_results.tsv'
with open(results_file, 'w') as f:
    f.write("seq_len\tn_copies\tn_shuffles\ttotal_hits\tunique_positions\thits_per_position\t")
    f.write("unique_tes\tpident_mean\tpident_std\tlength_mean\tlength_std\n")
    for r in results:
        f.write(f"{r['seq_len']}\t{r['n_copies']}\t{r['n_shuffles']}\t{r['total_hits']}\t")
        f.write(f"{r['unique_positions']}\t{r['hits_per_position']:.3f}\t{r['unique_tes']}\t")
        f.write(f"{r['pident_mean']:.2f}\t{r['pident_std']:.2f}\t{r['length_mean']:.1f}\t{r['length_std']:.1f}\n")
print(f"  Saved: {results_file}")

print("\n" + "=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)
