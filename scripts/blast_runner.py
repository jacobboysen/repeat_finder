#!/usr/bin/env python3
"""
Run BLAST queries with flexible configuration and CLI overrides.

Supports loading parameters from config file with command-line overrides.
Creates organized output directories with metadata and logs.
"""

import argparse
import logging
import subprocess
import sys
import yaml
from datetime import datetime
from pathlib import Path


class BlastRunner:
    """Manages BLAST execution with configuration management."""

    def __init__(self, config_file='config.yaml'):
        """Initialize with configuration file."""
        self.config_file = Path(config_file)
        self.config = self.load_config()

    def load_config(self):
        """Load configuration from YAML file."""
        if not self.config_file.exists():
            print(f"Warning: Config file {self.config_file} not found, using defaults")
            return self.default_config()

        with open(self.config_file) as f:
            return yaml.safe_load(f)

    def default_config(self):
        """Return default configuration."""
        return {
            'blast': {
                'program': 'blastn',
                'evalue': 1e-5,
                'word_size': 11,
                'dust': 'yes',
                'perc_identity': 0,
                'max_target_seqs': 500,
                'max_hsps': 10,
                'outfmt': '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen',
                'num_threads': 4
            },
            'query': {
                'segment_type': 'CDS',
                'gene_group': None
            },
            'subject': {
                'database': 'dmel_genome'
            },
            'paths': {
                'references': 'data/references',
                'blastdb': 'data/blastdb',
                'queries': 'data/queries',
                'results': 'results'
            }
        }

    def override_config(self, args):
        """Override configuration with command-line arguments."""
        # BLAST parameters
        blast_params = [
            'program', 'evalue', 'word_size', 'dust', 'perc_identity',
            'max_target_seqs', 'max_hsps', 'outfmt', 'num_threads'
        ]

        for param in blast_params:
            value = getattr(args, param, None)
            if value is not None:
                self.config['blast'][param] = value

        # Query parameters
        if args.segment_type:
            self.config['query']['segment_type'] = args.segment_type
        if args.gene_group:
            self.config['query']['gene_group'] = args.gene_group

        # Subject parameters
        if args.database:
            self.config['subject']['database'] = args.database

    def get_query_file(self):
        """Determine query file path based on configuration."""
        segment_type = self.config['query']['segment_type']
        gene_group = self.config['query']['gene_group']

        queries_dir = Path(self.config['paths']['queries'])

        # If gene group specified, use filtered queries
        if gene_group:
            group_dir = queries_dir / gene_group.replace(' ', '_')
            query_file = group_dir / f"{segment_type.lower()}.fasta"
        else:
            # Use full segment file from references
            refs_dir = Path(self.config['paths']['references'])
            file_map = {
                'CDS': 'dmel_cds.fasta',
                '5UTR': 'dmel_5utr.fasta',
                '3UTR': 'dmel_3utr.fasta',
                'intron': 'dmel_intron.fasta'
            }
            query_file = refs_dir / file_map.get(segment_type, 'dmel_cds.fasta')

        return query_file

    def get_database_path(self):
        """Get full path to BLAST database."""
        db_name = self.config['subject']['database']
        blastdb_dir = Path(self.config['paths']['blastdb'])
        return blastdb_dir / db_name

    def create_output_dir(self):
        """Create output directory with timestamp and metadata."""
        results_dir = Path(self.config['paths']['results'])

        segment_type = self.config['query']['segment_type']
        database = self.config['subject']['database']
        gene_group = self.config['query']['gene_group']

        # Create path: results/{segment_type}/{database}/{timestamp}_{gene_group}/
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        group_suffix = f"_{gene_group.replace(' ', '_')}" if gene_group else "_all"

        output_dir = results_dir / segment_type / database / f"{timestamp}{group_suffix}"
        output_dir.mkdir(parents=True, exist_ok=True)

        return output_dir

    def setup_logging(self, output_dir):
        """Set up logging to file and console."""
        log_file = output_dir / 'run.log'

        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )

        return logging.getLogger(__name__)

    def save_config(self, output_dir):
        """Save configuration used for this run."""
        config_file = output_dir / 'config_used.yaml'

        with open(config_file, 'w') as f:
            yaml.dump(self.config, f, default_flow_style=False)

        return config_file

    def get_query_stats(self, query_file):
        """Get statistics about query sequences."""
        from Bio import SeqIO

        if not query_file.exists():
            return None

        num_seqs = 0
        total_length = 0
        min_length = float('inf')
        max_length = 0

        for record in SeqIO.parse(query_file, 'fasta'):
            num_seqs += 1
            seq_len = len(record.seq)
            total_length += seq_len
            min_length = min(min_length, seq_len)
            max_length = max(max_length, seq_len)

        return {
            'num_seqs': num_seqs,
            'total_length': total_length,
            'mean_length': total_length / num_seqs if num_seqs > 0 else 0,
            'min_length': min_length if num_seqs > 0 else 0,
            'max_length': max_length
        }

    def save_query_stats(self, query_file, output_dir):
        """Save query statistics to file."""
        stats = self.get_query_stats(query_file)

        if not stats:
            return None

        stats_file = output_dir / 'query_stats.txt'

        with open(stats_file, 'w') as f:
            f.write(f"Query Statistics\n")
            f.write(f"{'=' * 60}\n")
            f.write(f"Query file: {query_file}\n")
            f.write(f"Number of sequences: {stats['num_seqs']}\n")
            f.write(f"Total length: {stats['total_length']:,} bp\n")
            f.write(f"Mean length: {stats['mean_length']:.1f} bp\n")
            f.write(f"Min length: {stats['min_length']} bp\n")
            f.write(f"Max length: {stats['max_length']} bp\n")

        return stats_file

    def build_blast_command(self, query_file, db_path, output_file):
        """Build BLAST command from configuration."""
        blast_config = self.config['blast']

        cmd = [
            blast_config['program'],
            '-query', str(query_file),
            '-db', str(db_path),
            '-out', str(output_file),
            '-evalue', str(blast_config['evalue']),
            '-word_size', str(blast_config['word_size']),
            '-dust', blast_config['dust'],
            '-max_target_seqs', str(blast_config['max_target_seqs']),
            '-max_hsps', str(blast_config['max_hsps']),
            '-num_threads', str(blast_config['num_threads']),
            '-outfmt', blast_config['outfmt']
        ]

        # Add perc_identity if > 0
        if blast_config['perc_identity'] > 0:
            cmd.extend(['-perc_identity', str(blast_config['perc_identity'])])

        return cmd

    def run_blast(self, logger):
        """Execute BLAST search."""
        # Get query file
        query_file = self.get_query_file()

        if not query_file.exists():
            logger.error(f"Query file not found: {query_file}")
            logger.info("If filtering by gene group, run filter_queries.py first")
            return 1

        # Get database path
        db_path = self.get_database_path()

        if not db_path.with_suffix('.nhr').exists():
            logger.error(f"Database not found: {db_path}")
            logger.info("Run build_databases.py to create databases")
            return 1

        # Create output directory
        output_dir = self.create_output_dir()
        logger.info(f"Output directory: {output_dir}")

        # Save configuration
        config_file = self.save_config(output_dir)
        logger.info(f"Configuration saved: {config_file}")

        # Save query statistics
        stats_file = self.save_query_stats(query_file, output_dir)
        if stats_file:
            logger.info(f"Query statistics saved: {stats_file}")

        # Build BLAST command
        output_file = output_dir / 'blast_results.tsv'
        cmd = self.build_blast_command(query_file, db_path, output_file)

        logger.info(f"Running BLAST...")
        logger.info(f"  Query: {query_file}")
        logger.info(f"  Database: {db_path}")
        logger.info(f"  Command: {' '.join(cmd)}")

        # Run BLAST
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )

            if result.stderr:
                logger.warning(f"BLAST stderr: {result.stderr}")

            logger.info(f"✓ BLAST completed successfully")
            logger.info(f"  Results: {output_file}")

            # Count results
            if output_file.exists():
                num_lines = sum(1 for _ in open(output_file))
                logger.info(f"  Hits: {num_lines}")

            return 0

        except subprocess.CalledProcessError as e:
            logger.error(f"BLAST failed with exit code {e.returncode}")
            logger.error(e.stderr)
            return 1
        except FileNotFoundError:
            logger.error("BLAST not found. Is BLAST+ installed?")
            logger.info("Install with: conda install -c bioconda blast")
            return 1


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # Config file
    parser.add_argument(
        '--config',
        type=str,
        default='config.yaml',
        help='Configuration file (default: config.yaml)'
    )

    # BLAST parameters
    parser.add_argument('--program', help='BLAST program (e.g., blastn, blastp)')
    parser.add_argument('--evalue', type=float, help='E-value threshold')
    parser.add_argument('--word-size', type=int, dest='word_size', help='Word size')
    parser.add_argument('--dust', choices=['yes', 'no'], help='Dust filter')
    parser.add_argument('--perc-identity', type=float, dest='perc_identity',
                        help='Percent identity threshold')
    parser.add_argument('--max-target-seqs', type=int, dest='max_target_seqs',
                        help='Maximum target sequences')
    parser.add_argument('--max-hsps', type=int, dest='max_hsps',
                        help='Maximum HSPs per subject')
    parser.add_argument('--num-threads', type=int, dest='num_threads',
                        help='Number of threads')
    parser.add_argument('--outfmt', help='Output format string')

    # Query parameters
    parser.add_argument(
        '--segment-type',
        dest='segment_type',
        choices=['CDS', '5UTR', '3UTR', 'intron'],
        help='Genic segment type'
    )
    parser.add_argument(
        '--gene-group',
        dest='gene_group',
        help='Gene group to query'
    )

    # Subject parameters
    parser.add_argument(
        '--database',
        help='BLAST database name (without path)'
    )

    args = parser.parse_args()

    # Create runner and load config
    runner = BlastRunner(args.config)

    # Override config with CLI arguments
    runner.override_config(args)

    # Create output directory and setup logging
    output_dir = runner.create_output_dir()
    logger = runner.setup_logging(output_dir)

    logger.info("=" * 60)
    logger.info("BLAST Runner")
    logger.info("=" * 60)

    # Run BLAST
    exit_code = runner.run_blast(logger)

    logger.info("=" * 60)

    return exit_code


if __name__ == '__main__':
    sys.exit(main())
