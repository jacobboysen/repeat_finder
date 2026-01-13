#!/usr/bin/env python3
"""
Generate publication-quality TE signal plots for each 3'UTR.

Creates multi-panel figures showing:
1. TE signal density line plot with cluster highlighting
2. Nucleotide-resolution hit map (colored by TE family)
3. Windowed density heatmap (50bp bins)
4. GC content track
"""

import argparse
import pickle
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from Bio import SeqIO


# Color palette for TE families
TE_COLORS = [
    '#e41a1c',  # red
    '#377eb8',  # blue
    '#4daf4a',  # green
    '#984ea3',  # purple
    '#ff7f00',  # orange
    '#ffff33',  # yellow
    '#a65628',  # brown
    '#f781bf',  # pink
    '#999999',  # gray
]


def load_density_data(density_file):
    """Load density data from pickle file."""
    with open(density_file, 'rb') as f:
        return pickle.load(f)


def calculate_gc_content(sequence, window_size=50):
    """Calculate windowed GC content."""
    seq = str(sequence).upper()
    gc = np.zeros(len(seq))

    for i in range(len(seq)):
        start = max(0, i - window_size // 2)
        end = min(len(seq), i + window_size // 2)
        window = seq[start:end]
        if len(window) > 0:
            gc_count = window.count('G') + window.count('C')
            gc[i] = gc_count / len(window)

    return gc


def assign_te_colors(hits):
    """Assign colors to TE families."""
    # Count TE families
    family_counts = defaultdict(int)
    for hit in hits:
        te_id = hit['sseqid']
        # Extract family (use full ID for now, could parse further)
        family_counts[te_id] += 1

    # Sort by count and assign colors
    sorted_families = sorted(family_counts.keys(), key=lambda x: -family_counts[x])
    color_map = {}
    for i, family in enumerate(sorted_families):
        color_map[family] = TE_COLORS[i % len(TE_COLORS)]

    return color_map


def plot_single_sequence(qseqid, data, sequence=None, output_dir=None, figsize=(14, 10)):
    """
    Create multi-panel TE signal plot for a single sequence.

    Args:
        qseqid: Query sequence ID
        data: Density data dictionary for this sequence
        sequence: Optional sequence string for GC content
        output_dir: Output directory for figure files
        figsize: Figure size tuple
    """
    length = data['length']
    positions = np.arange(length)

    # Create figure with subplots
    fig, axes = plt.subplots(4, 1, figsize=figsize, height_ratios=[2, 2, 0.5, 0.5],
                             sharex=True)

    # Panel 1: Signal density line plot
    ax1 = axes[0]
    ax1.fill_between(positions, data['combined'], alpha=0.3, color='#377eb8')
    ax1.plot(positions, data['combined'], color='#377eb8', linewidth=1)

    # Highlight peaks/clusters
    peaks = data.get('peaks', [])
    for peak in peaks:
        pos = peak['position'] - 1  # Convert to 0-indexed
        ax1.axvline(pos, color='red', alpha=0.5, linestyle='--', linewidth=0.8)

    ax1.set_ylabel('TE Signal\n(Combined)', fontsize=10)
    ax1.set_title(f'{qseqid} - TE Signal Density ({length} bp)', fontsize=12, fontweight='bold')
    ax1.set_xlim(0, length)

    # Use log scale if values span many orders of magnitude
    if data['combined'].max() > 0:
        if data['combined'].max() / (data['combined'][data['combined'] > 0].min() + 1e-10) > 1000:
            ax1.set_yscale('symlog', linthresh=0.01)

    # Panel 2: Nucleotide-resolution hit map
    ax2 = axes[1]
    hits = data.get('hits', [])

    if hits:
        # Assign colors to TE families
        color_map = assign_te_colors(hits)

        # Sort hits by start position for stacking
        sorted_hits = sorted(hits, key=lambda x: x['qstart'])

        # Simple stacking algorithm
        levels = []  # Track occupied y-ranges
        hit_y_positions = []

        for hit in sorted_hits:
            start = hit['qstart'] - 1
            end = hit['qend']

            # Find first non-overlapping level
            y_level = 0
            for i, (level_end, _) in enumerate(levels):
                if start >= level_end:
                    y_level = i
                    break
            else:
                y_level = len(levels)

            # Update or add level
            if y_level < len(levels):
                levels[y_level] = (end, hit)
            else:
                levels.append((end, hit))

            hit_y_positions.append((hit, y_level))

        # Draw hits as rectangles
        max_level = max(y for _, y in hit_y_positions) + 1 if hit_y_positions else 1

        for hit, y_level in hit_y_positions:
            start = hit['qstart'] - 1
            width = hit['qend'] - start
            color = color_map.get(hit['sseqid'], '#999999')

            # Opacity based on e-value (better = darker)
            evalue = hit.get('evalue', 1)
            alpha = min(0.9, max(0.2, -np.log10(max(evalue, 1e-50)) / 50))

            rect = mpatches.Rectangle(
                (start, y_level), width, 0.9,
                facecolor=color, alpha=alpha, edgecolor='none'
            )
            ax2.add_patch(rect)

        ax2.set_ylim(-0.5, max_level + 0.5)
        ax2.set_ylabel('BLAST Hits\n(Stacked)', fontsize=10)

        # Create legend for top TE families
        top_families = sorted(color_map.keys(),
                              key=lambda x: sum(1 for h in hits if h['sseqid'] == x),
                              reverse=True)[:5]
        legend_patches = [mpatches.Patch(color=color_map[f], label=f[:15])
                          for f in top_families]
        ax2.legend(handles=legend_patches, loc='upper right', fontsize=8,
                   ncol=min(3, len(legend_patches)))
    else:
        ax2.text(length/2, 0.5, 'No hits', ha='center', va='center',
                 fontsize=12, color='gray')
        ax2.set_ylim(0, 1)
        ax2.set_ylabel('BLAST Hits', fontsize=10)

    ax2.set_xlim(0, length)

    # Panel 3: Windowed density heatmap (50bp bins)
    ax3 = axes[2]
    window_size = 50
    n_windows = (length + window_size - 1) // window_size
    windowed = np.zeros(n_windows)

    for i in range(n_windows):
        start = i * window_size
        end = min((i + 1) * window_size, length)
        windowed[i] = np.mean(data['combined'][start:end])

    # Create heatmap
    windowed_2d = windowed.reshape(1, -1)
    cmap = LinearSegmentedColormap.from_list('te_signal', ['white', '#377eb8', '#d62728'])
    im = ax3.imshow(windowed_2d, aspect='auto', cmap=cmap,
                    extent=[0, length, 0, 1], interpolation='nearest')
    ax3.set_ylabel(f'Binned\n({window_size}bp)', fontsize=10)
    ax3.set_yticks([])

    # Panel 4: GC content
    ax4 = axes[3]
    if sequence:
        gc = calculate_gc_content(sequence)
        ax4.fill_between(positions, gc, alpha=0.5, color='#4daf4a')
        ax4.plot(positions, gc, color='#4daf4a', linewidth=0.5)
        ax4.axhline(0.5, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
        ax4.set_ylim(0, 1)
    else:
        ax4.text(length/2, 0.5, 'Sequence not available', ha='center', va='center',
                 fontsize=10, color='gray')
        ax4.set_ylim(0, 1)

    ax4.set_ylabel('GC\nContent', fontsize=10)
    ax4.set_xlabel('Position (bp)', fontsize=10)
    ax4.set_xlim(0, length)

    plt.tight_layout()

    # Save figure
    if output_dir:
        output_dir.mkdir(parents=True, exist_ok=True)
        # Clean up query ID for filename
        safe_name = qseqid.replace('/', '_').replace(':', '_')
        png_file = output_dir / f'{safe_name}_te_signal.png'
        pdf_file = output_dir / f'{safe_name}_te_signal.pdf'

        plt.savefig(png_file, dpi=150, bbox_inches='tight')
        plt.savefig(pdf_file, bbox_inches='tight')
        plt.close()

        return png_file, pdf_file
    else:
        plt.show()
        return None, None


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        'density_dir',
        type=Path,
        nargs='?',
        help='Directory containing density_data.pkl'
    )
    parser.add_argument(
        '--query-fasta',
        type=Path,
        default=Path('data/queries/germ_plasm/3UTR_sense_tier1.fasta'),
        help='Query FASTA for GC content calculation'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('figures/te_signal'),
        help='Output directory for figures'
    )
    parser.add_argument(
        '--sequence',
        type=str,
        help='Plot only this sequence ID'
    )
    parser.add_argument(
        '--all',
        action='store_true',
        help='Plot all sequences'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )

    args = parser.parse_args()

    # Find density data
    if args.density_dir is None:
        sweep_dir = Path('results/parameter_sweep')
        if sweep_dir.exists():
            subdirs = sorted([d for d in sweep_dir.iterdir() if d.is_dir()], reverse=True)
            if subdirs:
                combo_dirs = sorted([d for d in subdirs[0].iterdir()
                                    if d.is_dir() and d.name.startswith('combo_')])
                if combo_dirs:
                    density_dir = combo_dirs[0] / 'density'
                    if density_dir.exists():
                        args.density_dir = density_dir

        if args.density_dir is None:
            print("Error: No density directory specified and none found", file=sys.stderr)
            return 1

    density_file = args.density_dir / 'density_data.pkl'
    if not density_file.exists():
        print(f"Error: Density file not found: {density_file}", file=sys.stderr)
        return 1

    print("TE Signal Plotting")
    print("=" * 60)
    print(f"Density data: {args.density_dir}")
    print(f"Output: {args.output_dir}")
    print()

    # Load density data
    print("Loading density data...")
    density_dict = load_density_data(density_file)
    print(f"  Loaded data for {len(density_dict)} sequences")

    # Load sequences for GC content
    sequences = {}
    if args.query_fasta.exists():
        print(f"Loading sequences from: {args.query_fasta}")
        for record in SeqIO.parse(args.query_fasta, 'fasta'):
            sequences[record.id] = str(record.seq)

    # Determine which sequences to plot
    if args.sequence:
        to_plot = [args.sequence] if args.sequence in density_dict else []
        if not to_plot:
            print(f"Warning: Sequence '{args.sequence}' not found in density data")
    elif args.all:
        to_plot = list(density_dict.keys())
    else:
        # Default: plot top sequences by signal
        by_signal = sorted(density_dict.keys(),
                          key=lambda x: np.max(density_dict[x]['combined']),
                          reverse=True)
        to_plot = by_signal[:5]

    print(f"\nPlotting {len(to_plot)} sequences...")

    for i, qseqid in enumerate(to_plot):
        print(f"  [{i+1}/{len(to_plot)}] {qseqid}...", end=' ')
        data = density_dict[qseqid]
        seq = sequences.get(qseqid)

        png_file, pdf_file = plot_single_sequence(
            qseqid, data, sequence=seq, output_dir=args.output_dir
        )

        if png_file:
            print(f"saved to {png_file.name}")
        else:
            print("displayed")

    print("\n" + "=" * 60)
    print("Plotting complete!")

    return 0


if __name__ == '__main__':
    sys.exit(main())
