import gzip
import os
import tempfile
from pathlib import Path

from modelseedpy.core.msgenome import MSGenome


def _parse_gff_attributes(attr_str: str) -> dict:
    attrs = {}
    for part in attr_str.strip().split(';'):
        part = part.strip()
        if '=' in part:
            k, v = part.split('=', 1)
            attrs[k.strip()] = v.strip()
    return attrs


class ExtractContainerGenome:

    def __init__(self):
        self.contigs = None   # MSGenome — contig sequences from .fna
        self.proteins = None  # MSGenome — protein sequences from .faa
        self.features = None  # list[dict] — all GFF feature records


class ExtractGenomeNCBI:

    def __init__(self):
        pass

    def extract(self, genome_fna: Path, genome_faa: Path, genome_gff: Path) -> ExtractContainerGenome:
        # Load contig sequences
        contigs = MSGenome.from_fasta(str(genome_fna))

        # Load protein sequences
        proteins = MSGenome.from_fasta(str(genome_faa))

        # Parse all features from GFF — no filtering, no transformation
        features = []
        open_fn = gzip.open if str(genome_gff).endswith('.gz') else open
        with open_fn(genome_gff, 'rt') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 9:
                    continue
                attrs = _parse_gff_attributes(parts[8])
                features.append({
                    'contig_id': parts[0],
                    'source': parts[1],
                    'feature_type': parts[2],
                    'start': int(parts[3]),
                    'end': int(parts[4]),
                    'score': parts[5],
                    'strand': parts[6],
                    'phase': parts[7],
                    'attrs': attrs,
                })

        container = ExtractContainerGenome()
        container.contigs = contigs
        container.proteins = proteins
        container.features = features
        return container
