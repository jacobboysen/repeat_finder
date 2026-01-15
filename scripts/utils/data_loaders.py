"""
Data loading utilities for the TE fossil mining pipeline.

Provides consistent loading of gene lists, FASTA files, and TE databases
across all analysis scripts.
"""

import re
from pathlib import Path
from typing import Dict, List, Optional, Union

try:
    from Bio import SeqIO
except ImportError:
    SeqIO = None


def load_gene_list(path: Union[str, Path]) -> List[str]:
    """
    Load FBgn IDs from a gene list file.

    Handles both simple ID files (*_fbgn_ids.txt) and TSV files.

    Args:
        path: Path to gene list file

    Returns:
        List of FBgn IDs
    """
    path = Path(path)
    genes = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            # Handle TSV files - look for FBgn in any column
            if '\t' in line:
                parts = line.split('\t')
                for part in parts:
                    if part.startswith('FBgn'):
                        genes.append(part)
                        break
            elif line.startswith('FBgn'):
                genes.append(line)

    return genes


def load_gene_list_with_symbols(path: Union[str, Path]) -> Dict[str, str]:
    """
    Load gene list with symbols from consolidated TSV file.

    Expected format: symbol<TAB>FBgn<TAB>...

    Args:
        path: Path to consolidated gene list TSV

    Returns:
        Dictionary mapping FBgn ID to gene symbol
    """
    path = Path(path)
    genes = {}

    with open(path) as f:
        # Skip header if present
        first_line = f.readline()
        if not first_line.startswith('FBgn') and '\t' in first_line:
            # This was a header, continue
            pass
        else:
            # Not a header, process this line
            f.seek(0)

        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split('\t')
            if len(parts) >= 2:
                # Find FBgn and symbol
                fbgn = None
                symbol = None

                for i, part in enumerate(parts):
                    if part.startswith('FBgn'):
                        fbgn = part
                        # Symbol is usually the first non-FBgn column
                        for j, p in enumerate(parts):
                            if j != i and not p.startswith('FBgn') and not p.startswith('FBtr'):
                                symbol = p
                                break
                        break

                if fbgn:
                    genes[fbgn] = symbol or fbgn

    return genes


def parse_fasta(path: Union[str, Path]) -> Dict[str, str]:
    """
    Parse FASTA file into dictionary of ID -> sequence.

    Args:
        path: Path to FASTA file

    Returns:
        Dictionary mapping sequence ID to sequence string
    """
    path = Path(path)
    sequences = {}
    current_id = None
    current_seq = []

    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                # Extract ID (first whitespace-separated token, without >)
                current_id = line.strip().split()[0][1:]
                current_seq = []
            else:
                current_seq.append(line.strip())

        # Don't forget the last sequence
        if current_id:
            sequences[current_id] = ''.join(current_seq)

    return sequences


def parse_fasta_by_parent(path: Union[str, Path]) -> Dict[str, Dict]:
    """
    Parse FASTA file and extract parent gene IDs.

    For UTR FASTAs with headers like:
    >FBtr0070000 type=three_prime_UTR; loc=...; length=1019; parent=FBgn0031081;

    Args:
        path: Path to FASTA file

    Returns:
        Dictionary mapping transcript ID to dict with 'sequence', 'length', 'parent'
    """
    path = Path(path)
    transcripts = {}
    current_id = None
    current_info = {}
    current_seq = []

    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if current_id:
                    current_info['sequence'] = ''.join(current_seq)
                    transcripts[current_id] = current_info

                parts = line.strip().split()
                current_id = parts[0][1:]  # Remove >
                current_info = {'length': None, 'parent': None}
                current_seq = []

                # Parse header attributes
                for part in parts[1:]:
                    if part.startswith('length='):
                        try:
                            current_info['length'] = int(part.split('=')[1].rstrip(';'))
                        except (ValueError, IndexError):
                            pass
                    elif part.startswith('parent='):
                        current_info['parent'] = part.split('=')[1].rstrip(';')
            else:
                current_seq.append(line.strip())

        # Don't forget the last sequence
        if current_id:
            current_info['sequence'] = ''.join(current_seq)
            transcripts[current_id] = current_info

    return transcripts


def parse_te_name(description: str) -> Optional[str]:
    """
    Extract TE family name from FlyBase TE description.

    Args:
        description: FASTA header description

    Returns:
        TE family name or None
    """
    # Look for name= pattern
    match = re.search(r'name=([^;]+)', description)
    if match:
        name = match.group(1)
        # Clean up the name - remove instance identifiers like {}555
        name = re.sub(r'\{[^}]*\}\d*', '', name)
        return name.strip()
    return None


def parse_te_class(name: str) -> str:
    """
    Classify TE by type based on common naming patterns.

    Args:
        name: TE family name

    Returns:
        TE class (LTR, LINE, DNA, SINE, or Unknown)
    """
    name_lower = name.lower()

    # LTR retrotransposons
    ltr_families = [
        'gypsy', 'copia', 'bel', 'pao', 'mdg', 'roo', '412', '297',
        'blood', 'accord', 'tirant', 'springer', 'opus', 'diver',
        'quasimodo', 'idefix', 'invader', 'gtwin', 'tabor', 'stalker'
    ]
    for fam in ltr_families:
        if fam in name_lower:
            return 'LTR'

    # Non-LTR retrotransposons (LINEs)
    line_families = [
        'jockey', 'doc', 'i-element', 'f-element', 'g-element',
        'x-element', 'het-a', 'tart', 'tahre', 'r1', 'r2', 'cr1',
        'rt1', 'baggins', 'juan', 'ivk'
    ]
    for fam in line_families:
        if fam in name_lower:
            return 'LINE'

    # DNA transposons
    dna_families = [
        'p-element', 'hobo', 'pogo', 'bari', 's-element',
        'transib', 'tc1', 'mariner', 'piggybac', 'helitron',
        'mite', 'protop', 'dnarep'
    ]
    for fam in dna_families:
        if fam in name_lower:
            return 'DNA'

    # SINEs
    if 'sine' in name_lower or name_lower.startswith('ine-'):
        return 'SINE'

    # Terminal repeats
    if 'ltr' in name_lower:
        return 'LTR'

    return 'Unknown'


def load_te_database(te_fasta: Union[str, Path]) -> Dict[str, Dict]:
    """
    Load TE database from FASTA file.

    Args:
        te_fasta: Path to TE database FASTA

    Returns:
        Dictionary mapping TE ID to dict with 'name', 'class', 'length'
    """
    te_fasta = Path(te_fasta)
    te_info = {}

    if SeqIO is not None:
        # Use BioPython if available
        for record in SeqIO.parse(te_fasta, 'fasta'):
            te_id = record.id
            name = parse_te_name(record.description)
            if name:
                te_class = parse_te_class(name)
                te_info[te_id] = {
                    'name': name,
                    'class': te_class,
                    'length': len(record.seq)
                }
            else:
                te_info[te_id] = {
                    'name': te_id,
                    'class': 'Unknown',
                    'length': len(record.seq)
                }
    else:
        # Fallback: parse manually
        sequences = parse_fasta(te_fasta)

        with open(te_fasta) as f:
            for line in f:
                if line.startswith('>'):
                    parts = line.strip().split()
                    te_id = parts[0][1:]
                    name = parse_te_name(line) or te_id

                    te_info[te_id] = {
                        'name': name,
                        'class': parse_te_class(name),
                        'length': len(sequences.get(te_id, ''))
                    }

    return te_info


def load_strand_bias_data(path: Union[str, Path]) -> Dict[str, Dict]:
    """
    Load strand bias data, aggregated by gene.

    Args:
        path: Path to strand_bias_by_utr.tsv file

    Returns:
        Dictionary mapping FBgn to dict with hit counts and strand percentages
    """
    path = Path(path)
    gene_data = {}

    with open(path) as f:
        # Skip header
        f.readline()

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue

            fbtr = parts[0]
            fbgn = parts[1]
            utr_len = int(parts[2])
            total_hits = int(parts[3])
            sense_hits = int(parts[4])
            anti_hits = int(parts[5])

            if fbgn not in gene_data:
                gene_data[fbgn] = {
                    'total_hits': 0,
                    'sense_hits': 0,
                    'anti_hits': 0,
                    'total_utr_len': 0
                }

            gene_data[fbgn]['total_hits'] += total_hits
            gene_data[fbgn]['sense_hits'] += sense_hits
            gene_data[fbgn]['anti_hits'] += anti_hits
            gene_data[fbgn]['total_utr_len'] += utr_len

    # Calculate percentages
    for fbgn, data in gene_data.items():
        total = data['total_hits']
        if total > 0:
            data['sense_pct'] = 100 * data['sense_hits'] / total
            data['anti_pct'] = 100 * data['anti_hits'] / total
        else:
            data['sense_pct'] = 50
            data['anti_pct'] = 50

    return gene_data


def load_gene_te_density(path: Union[str, Path]) -> Dict[str, Dict]:
    """
    Load gene-level TE density data from top/bottom 100 files.

    Args:
        path: Path to top_100_te_genes.tsv or similar

    Returns:
        Dictionary mapping FBgn to dict with rank, density, hits, etc.
    """
    path = Path(path)
    gene_data = {}

    with open(path) as f:
        # Skip header
        f.readline()

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue

            rank = int(parts[0])
            fbgn = parts[1]
            density = float(parts[2])
            hits = int(parts[3])
            hit_bp = int(parts[4])
            utr_len = int(parts[5])

            gene_data[fbgn] = {
                'rank': rank,
                'density': density,
                'hits': hits,
                'hit_bp': hit_bp,
                'utr_len': utr_len
            }

    return gene_data
