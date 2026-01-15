#!/usr/bin/env python3
"""
Get detailed information for top ancient TE candidates including:
- Gene names and functions
- Sequence alignments
- TE family information
"""

import sys
import os
import re
from collections import defaultdict

def load_gene_info(gff_file):
    """Load gene/transcript info from GFF."""
    gene_info = {}  # transcript -> gene symbol
    gene_to_fbgn = {}  # symbol -> FBgn

    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            if parts[2] == 'mRNA':
                attrs = dict(a.split('=') for a in parts[8].split(';') if '=' in a)
                if 'ID' in attrs and 'Name' in attrs:
                    transcript_id = attrs['ID']
                    name = attrs['Name']
                    # Extract gene symbol (before -R*)
                    gene_symbol = re.sub(r'-R[A-Z]+$', '', name)
                    gene_info[transcript_id] = gene_symbol

                if 'Parent' in attrs and 'Name' in attrs:
                    parent = attrs['Parent']
                    name = attrs['Name']
                    gene_symbol = re.sub(r'-R[A-Z]+$', '', name)
                    gene_to_fbgn[gene_symbol] = parent

    return gene_info, gene_to_fbgn

def load_te_info(te_fasta):
    """Load TE family names from FASTA headers."""
    te_info = {}
    with open(te_fasta) as f:
        for line in f:
            if line.startswith('>'):
                parts = line.split()
                te_id = parts[0][1:]  # Remove >
                # Extract name from header
                name_match = re.search(r'name=([^;]+)', line)
                if name_match:
                    te_info[te_id] = name_match.group(1)
                else:
                    te_info[te_id] = te_id
    return te_info

def load_blast_hits(novel_file, known_file, target_pairs):
    """Load BLAST alignment data for specific transcript-TE pairs."""
    hits = {}

    for filepath in [novel_file, known_file]:
        with open(filepath) as f:
            header = next(f)
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 7:
                    continue

                transcript = parts[0]
                te_id = parts[3]
                key = (transcript, te_id)

                if key in target_pairs:
                    if key not in hits:
                        hits[key] = {
                            'transcript': transcript,
                            'te_id': te_id,
                            'qstart': int(parts[1]),
                            'qend': int(parts[2]),
                            'pident': float(parts[4]),
                            'length': int(parts[5]),
                            'evalue': parts[6]
                        }
    return hits

def load_utr_sequences(fasta_file):
    """Load UTR sequences from FASTA."""
    sequences = {}
    current_id = None
    current_seq = []

    with open(fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())

        if current_id:
            sequences[current_id] = ''.join(current_seq)

    return sequences

def main():
    # Load reference data
    print("Loading reference data...", file=sys.stderr)

    gff_file = '/Users/jacobboysen/git_repos/repeat_finder/data/references/dmel-all-r6.66.gff'
    te_fasta = '/Users/jacobboysen/git_repos/repeat_finder/data/references/dmel_te_flybase.fasta'
    utr_fasta = '/Users/jacobboysen/git_repos/repeat_finder/data/references/dmel_3utr.fasta'

    gene_info, gene_fbgn = load_gene_info(gff_file)
    te_info = load_te_info(te_fasta)
    utr_seqs = load_utr_sequences(utr_fasta)

    print(f"  Loaded {len(gene_info)} transcript->gene mappings", file=sys.stderr)
    print(f"  Loaded {len(te_info)} TE names", file=sys.stderr)
    print(f"  Loaded {len(utr_seqs)} UTR sequences", file=sys.stderr)

    # Load top candidates
    candidates_file = '/Users/jacobboysen/git_repos/repeat_finder/results/repeatmasker_analysis/ancient_te_candidates_top100.tsv'

    print("\nLoading top candidates...", file=sys.stderr)
    candidates = []
    target_pairs = set()

    with open(candidates_file) as f:
        header = next(f).strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            row = dict(zip(header, parts))
            candidates.append(row)
            target_pairs.add((row['transcript'], row['te_id']))

    print(f"  Loaded {len(candidates)} top candidates", file=sys.stderr)

    # Load BLAST hit details
    print("\nLoading BLAST alignments...", file=sys.stderr)
    novel_file = '/Users/jacobboysen/git_repos/repeat_finder/results/repeatmasker_analysis/blast_hits_novel.tsv'
    known_file = '/Users/jacobboysen/git_repos/repeat_finder/results/repeatmasker_analysis/blast_hits_known_repeatmasker.tsv'

    blast_hits = load_blast_hits(novel_file, known_file, target_pairs)
    print(f"  Found BLAST data for {len(blast_hits)} candidates", file=sys.stderr)

    # Generate detailed report
    print("\n" + "="*80)
    print("TOP 50 ANCIENT FUNCTIONAL TE CANDIDATES - DETAILED VIEW")
    print("="*80)

    for i, cand in enumerate(candidates[:50], 1):
        transcript = cand['transcript']
        te_id = cand['te_id']
        phyloP = float(cand['phyloP'])
        syn_sp = int(cand['syntenic_species'])
        pident = float(cand['pident'])
        length = int(cand['length'])
        chrom = cand['chrom']
        start = cand['start']
        end = cand['end']

        # Get gene name
        gene_symbol = gene_info.get(transcript, 'Unknown')
        fbgn = gene_fbgn.get(gene_symbol, '')

        # Get TE family name
        te_name = te_info.get(te_id, te_id)

        # Get BLAST details
        blast = blast_hits.get((transcript, te_id), {})
        qstart = blast.get('qstart', '?')
        qend = blast.get('qend', '?')
        evalue = blast.get('evalue', '?')

        # Get UTR sequence context
        utr_seq = utr_seqs.get(transcript, '')
        if utr_seq and isinstance(qstart, int) and isinstance(qend, int):
            # Extract hit region with context
            ctx_start = max(0, qstart - 11)
            ctx_end = min(len(utr_seq), qend + 10)
            seq_context = utr_seq[ctx_start:ctx_end]
            hit_seq = utr_seq[qstart-1:qend] if qstart > 0 else ''
        else:
            seq_context = ''
            hit_seq = ''

        print(f"\n{'─'*80}")
        print(f"#{i} | {gene_symbol} ({transcript}) | {te_name}")
        print(f"{'─'*80}")
        print(f"  Gene: {gene_symbol} ({fbgn})")
        print(f"  Transcript: {transcript}")
        print(f"  TE: {te_name} ({te_id})")
        print(f"  Location: {chrom}:{start}-{end}")
        print(f"  UTR position: {qstart}-{qend}")
        print(f"  Conservation: phyloP = {phyloP:.2f} (HIGHLY CONSERVED)")
        print(f"  Synteny: {syn_sp} species (sim={cand['sim_cov']}, yak={cand['yak_cov']}, ere={cand['ere_cov']})")
        print(f"  Alignment: {pident:.1f}% identity, {length}bp, E={evalue}")

        if hit_seq:
            print(f"\n  Sequence ({len(hit_seq)}bp):")
            # Print sequence in chunks
            for j in range(0, len(hit_seq), 60):
                chunk = hit_seq[j:j+60]
                print(f"    {j+1:4d}  {chunk}")

    # Summary by gene
    print("\n" + "="*80)
    print("GENES WITH MOST ANCIENT TE CANDIDATES (Top 100)")
    print("="*80)

    gene_counts = defaultdict(list)
    for cand in candidates:
        gene = gene_info.get(cand['transcript'], cand['transcript'])
        gene_counts[gene].append(cand)

    print(f"\n{'Gene':<20} {'Candidates':>10} {'Best phyloP':>12} {'TE families':<30}")
    print("-" * 75)

    for gene, cands in sorted(gene_counts.items(), key=lambda x: -len(x[1]))[:20]:
        best_phylo = max(float(c['phyloP']) for c in cands)
        te_fams = set(te_info.get(c['te_id'], c['te_id'])[:15] for c in cands)
        te_str = ', '.join(list(te_fams)[:3])
        if len(te_fams) > 3:
            te_str += f' +{len(te_fams)-3}'
        print(f"{gene:<20} {len(cands):>10} {best_phylo:>12.2f} {te_str:<30}")

    # TE family summary
    print("\n" + "="*80)
    print("TE FAMILIES IN TOP 100 ANCIENT CANDIDATES")
    print("="*80)

    te_family_counts = defaultdict(int)
    for cand in candidates:
        te_name = te_info.get(cand['te_id'], cand['te_id'])
        # Extract family (before {})
        family = re.sub(r'\{[^}]*\}.*', '', te_name)
        te_family_counts[family] += 1

    print(f"\n{'TE Family':<20} {'Count':>10}")
    print("-" * 35)
    for fam, count in sorted(te_family_counts.items(), key=lambda x: -x[1])[:15]:
        print(f"{fam:<20} {count:>10}")

if __name__ == '__main__':
    main()
