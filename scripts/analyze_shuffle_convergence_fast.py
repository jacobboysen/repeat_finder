#!/usr/bin/env python3
"""
Analyze shuffle convergence for duplicated sequences - OPTIMIZED VERSION.

Loads all BLAST data once, then analyzes.
"""

from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 10

print("=" * 70)
print("SHUFFLE CONVERGENCE ANALYSIS (OPTIMIZED)")
print("=" * 70)

# Step 1: Load sequences and find duplicates
print("\n[1/4] Loading sequences and finding duplicates...")
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

# Group by content and create mapping
by_content = defaultdict(list)
for seq_id, seq in seqs.items():
    by_content[seq].append(seq_id)

# Create transcript -> sequence_group mapping
transcript_to_group = {}
group_info = {}  # group_id -> {seq_len, n_copies, transcript_ids}
for group_id, (seq, ids) in enumerate(by_content.items()):
    if len(ids) >= 2:  # Only duplicates
        group_info[group_id] = {
            'seq_len': len(seq),
            'n_copies': len(ids),
            'transcript_ids': set(ids)
        }
        for tid in ids:
            transcript_to_group[tid] = group_id

print(f"  Total sequences: {len(seqs):,}")
print(f"  Duplicated groups (2+ copies): {len(group_info):,}")
print(f"  Transcripts in duplicated groups: {len(transcript_to_group):,}")

# Step 2: Load ALL BLAST hits and index by group
print("\n[2/4] Loading all BLAST hits (this takes ~2-3 min)...")

# Store hits by group: group_id -> list of hits
hits_by_group = defaultdict(list)

for rep in range(1, 11):
    blast_file = f'results/shuffled_full/replicate_{rep:02d}_blast.tsv'
    print(f"  Loading replicate {rep}...", end=' ', flush=True)
    count = 0
    with open(blast_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            qseqid = parts[0]
            base_id = qseqid.rsplit('_shuf', 1)[0]

            if base_id in transcript_to_group:
                group_id = transcript_to_group[base_id]
                te_id = parts[1]
                pident = float(parts[2])
                length = int(parts[3])
                sstart = int(parts[8])
                send = int(parts[9])

                hits_by_group[group_id].append({
                    'version': qseqid,
                    'te': te_id,
                    'pos': (te_id, min(sstart, send), max(sstart, send)),
                    'pident': pident,
                    'length': length
                })
                count += 1
    print(f"{count:,} hits from duplicated seqs")

print(f"  Total hits from duplicated sequences: {sum(len(h) for h in hits_by_group.values()):,}")

# Step 3: Analyze each group
print("\n[3/4] Analyzing convergence for each duplicated sequence...")

results = []
for group_id, info in group_info.items():
    hits = hits_by_group.get(group_id, [])
    if not hits:
        continue

    n_copies = info['n_copies']
    n_shuffles = n_copies * 10
    seq_len = info['seq_len']

    # Hits per version
    hits_per_version = defaultdict(list)
    for h in hits:
        hits_per_version[h['version']].append(h)

    hit_counts = [len(v) for v in hits_per_version.values()]

    # Position convergence
    all_positions = set(h['pos'] for h in hits)
    position_counts = Counter(h['pos'] for h in hits)

    # Identity and length
    pidents = [h['pident'] for h in hits]
    lengths = [h['length'] for h in hits]

    # Unique TEs
    unique_tes = set(h['te'] for h in hits)

    results.append({
        'seq_len': seq_len,
        'n_copies': n_copies,
        'n_shuffles': n_shuffles,
        'versions_with_hits': len(hits_per_version),
        'total_hits': len(hits),
        'unique_positions': len(all_positions),
        'unique_tes': len(unique_tes),
        'hits_per_position': len(hits) / len(all_positions) if all_positions else 0,
        'hit_count_min': min(hit_counts) if hit_counts else 0,
        'hit_count_max': max(hit_counts) if hit_counts else 0,
        'hit_count_mean': np.mean(hit_counts) if hit_counts else 0,
        'hit_count_std': np.std(hit_counts) if hit_counts else 0,
        'pident_mean': np.mean(pidents),
        'pident_std': np.std(pidents),
        'pident_min': min(pidents),
        'pident_max': max(pidents),
        'length_mean': np.mean(lengths),
        'length_std': np.std(lengths),
        'length_min': min(lengths),
        'length_max': max(lengths),
        'max_position_hits': max(position_counts.values()) if position_counts else 0,
    })

print(f"  Analyzed {len(results)} duplicated sequences with hits")

# Summary statistics
print("\n" + "=" * 70)
print("SUMMARY STATISTICS")
print("=" * 70)

n_copies_arr = np.array([r['n_copies'] for r in results])
n_shuffles_arr = np.array([r['n_shuffles'] for r in results])
total_hits_arr = np.array([r['total_hits'] for r in results])
unique_pos_arr = np.array([r['unique_positions'] for r in results])
hits_per_pos_arr = np.array([r['hits_per_position'] for r in results])
pident_std_arr = np.array([r['pident_std'] for r in results])
length_std_arr = np.array([r['length_std'] for r in results])
hit_count_std_arr = np.array([r['hit_count_std'] for r in results])

print(f"\nPosition Convergence (hits per unique position):")
print(f"  Mean: {np.mean(hits_per_pos_arr):.3f}")
print(f"  Median: {np.median(hits_per_pos_arr):.3f}")
print(f"  Min: {np.min(hits_per_pos_arr):.3f}")
print(f"  Max: {np.max(hits_per_pos_arr):.3f}")
print(f"  % with ratio < 1.05: {100*np.mean(hits_per_pos_arr < 1.05):.1f}%")
print(f"  % with ratio < 1.1: {100*np.mean(hits_per_pos_arr < 1.1):.1f}%")
print(f"  % with ratio < 1.5: {100*np.mean(hits_per_pos_arr < 1.5):.1f}%")

print(f"\nHit count spread (std dev across shuffled versions):")
print(f"  Mean std: {np.mean(hit_count_std_arr):.1f} hits")
print(f"  Max std: {np.max(hit_count_std_arr):.1f} hits")

print(f"\nIdentity spread (std dev within each sequence's shuffles):")
print(f"  Mean std: {np.mean(pident_std_arr):.2f}%")
print(f"  Max std: {np.max(pident_std_arr):.2f}%")

print(f"\nLength spread (std dev within each sequence's shuffles):")
print(f"  Mean std: {np.mean(length_std_arr):.1f} bp")
print(f"  Max std: {np.max(length_std_arr):.1f} bp")

# Step 4: Create figures
print("\n[4/4] Generating figures...")
fig_dir = Path('figures/shuffle_convergence')
fig_dir.mkdir(parents=True, exist_ok=True)

fig, axes = plt.subplots(2, 3, figsize=(15, 10))

# Panel 1: Convergence ratio distribution
ax = axes[0, 0]
ax.hist(hits_per_pos_arr, bins=50, color='steelblue', edgecolor='white', alpha=0.8)
ax.axvline(1.0, color='red', linestyle='--', linewidth=2, label='Perfect divergence (1.0)')
ax.axvline(np.median(hits_per_pos_arr), color='orange', linestyle='-', linewidth=2,
           label=f'Median: {np.median(hits_per_pos_arr):.3f}')
ax.set_xlabel('Hits per Unique TE Position')
ax.set_ylabel('Number of Duplicated Sequences')
ax.set_title('A. Position Convergence\n(1.0 = every hit at new position)')
ax.legend(fontsize=9)
ax.set_xlim(0.95, min(3, np.percentile(hits_per_pos_arr, 99)))

# Panel 2: Convergence vs number of shuffles
ax = axes[0, 1]
ax.scatter(n_shuffles_arr, hits_per_pos_arr, alpha=0.4, s=15, c='steelblue')
ax.axhline(1.0, color='red', linestyle='--', linewidth=1)
ax.set_xlabel('Number of Shuffles (copies × 10 reps)')
ax.set_ylabel('Hits per Unique Position')
ax.set_title('B. Convergence vs Shuffle Count')
ax.set_ylim(0.95, min(3, np.percentile(hits_per_pos_arr, 99)))

# Panel 3: Hit count spread (CV)
ax = axes[0, 2]
hit_count_cv = hit_count_std_arr / np.array([r['hit_count_mean'] for r in results])
hit_count_cv = hit_count_cv[~np.isnan(hit_count_cv) & ~np.isinf(hit_count_cv)]
ax.hist(hit_count_cv, bins=50, color='coral', edgecolor='white', alpha=0.8)
ax.axvline(np.median(hit_count_cv), color='red', linestyle='-', linewidth=2,
           label=f'Median CV: {np.median(hit_count_cv):.2f}')
ax.set_xlabel('Hit Count CV (std/mean)')
ax.set_ylabel('Number of Sequences')
ax.set_title('C. Hit Count Variability Across Shuffles')
ax.legend(fontsize=9)

# Panel 4: Identity spread distribution
ax = axes[1, 0]
ax.hist(pident_std_arr, bins=50, color='seagreen', edgecolor='white', alpha=0.8)
ax.axvline(np.mean(pident_std_arr), color='red', linestyle='-', linewidth=2,
           label=f'Mean: {np.mean(pident_std_arr):.1f}%')
ax.set_xlabel('Identity Std Dev (%)')
ax.set_ylabel('Number of Sequences')
ax.set_title('D. Identity Spread Across Shuffles')
ax.legend(fontsize=9)

# Panel 5: Length spread distribution
ax = axes[1, 1]
ax.hist(length_std_arr, bins=50, color='purple', edgecolor='white', alpha=0.8)
ax.axvline(np.mean(length_std_arr), color='red', linestyle='-', linewidth=2,
           label=f'Mean: {np.mean(length_std_arr):.1f}bp')
ax.set_xlabel('Length Std Dev (bp)')
ax.set_ylabel('Number of Sequences')
ax.set_title('E. Length Spread Across Shuffles')
ax.legend(fontsize=9)

# Panel 6: Summary
ax = axes[1, 2]
ax.axis('off')
summary = f"""
CONVERGENCE ANALYSIS SUMMARY

Duplicated sequences analyzed: {len(results):,}

POSITION CONVERGENCE:
  Median hits/position: {np.median(hits_per_pos_arr):.3f}
  {100*np.mean(hits_per_pos_arr < 1.05):.0f}% show ratio < 1.05
  {100*np.mean(hits_per_pos_arr < 1.1):.0f}% show ratio < 1.1

  → Nearly every shuffle hit lands at a
    NEW TE position, not a repeated one

HIT COUNT VARIABILITY:
  Median CV: {np.median(hit_count_cv):.2f}
  → Large variation in # hits per shuffle

IDENTITY VARIATION:
  Mean std: {np.mean(pident_std_arr):.1f}%
  → Moderate spread in alignment quality

LENGTH VARIATION:
  Mean std: {np.mean(length_std_arr):.0f} bp
  → Substantial spread in hit lengths

CONCLUSION:
Shuffles do NOT converge - each shuffle
explores different TE sequence space.
10 shuffles is not saturating.
"""
ax.text(0.05, 0.95, summary, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

plt.tight_layout()
plt.savefig(fig_dir / 'shuffle_convergence_analysis.png', dpi=150, bbox_inches='tight',
            facecolor='white')
print(f"  Saved: {fig_dir}/shuffle_convergence_analysis.png")
plt.close()

# Save results
results_file = fig_dir / 'convergence_results.tsv'
with open(results_file, 'w') as f:
    f.write("seq_len\tn_copies\tn_shuffles\ttotal_hits\tunique_positions\thits_per_position\t")
    f.write("hit_count_mean\thit_count_std\tpident_mean\tpident_std\tlength_mean\tlength_std\n")
    for r in results:
        f.write(f"{r['seq_len']}\t{r['n_copies']}\t{r['n_shuffles']}\t{r['total_hits']}\t")
        f.write(f"{r['unique_positions']}\t{r['hits_per_position']:.4f}\t")
        f.write(f"{r['hit_count_mean']:.1f}\t{r['hit_count_std']:.1f}\t")
        f.write(f"{r['pident_mean']:.2f}\t{r['pident_std']:.2f}\t")
        f.write(f"{r['length_mean']:.1f}\t{r['length_std']:.1f}\n")
print(f"  Saved: {results_file}")

print("\n" + "=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)
