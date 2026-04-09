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


def example_old_stuff():
    from modelseedpy_ext.biodb.ncbi import NCBIAssembly
    from modelseedpy_ext.biodb.ncbi_eutils import NcbiEutils
    from modelseedpy_ext.biodb.ncbi import NCBIAssembly

    client_eutils = NcbiEutils("da5f61d75e802e4effe119ac9868992b5a08")

    query = 'GCF_000146045.2'
    refseq_id = list(client_eutils.esearch('assembly', query)['ids'])[0]
    print(f'{query} has id {refseq_id}')

    """
    TranslationSet None
    TranslationStack    
    GCF_000146045.2 has id 285498
    """

    assembly_data = client_eutils.esummary('assembly', refseq_id)[0][0]

    print(assembly_data)

    """
    output:
    {'RsUid': 285498,
 'GbUid': 285798,
 'AssemblyAccession': 'GCF_000146045.2',
 'LastMajorReleaseAccession': 'GCF_000146045.2',
 'LatestAccession': 'None',
 'ChainId': 146045,
 'AssemblyName': 'R64',
 'UCSCName': 'sacCer3',
 'EnsemblName': 'None',
 'Taxid': 559292,
 'Organism': "Saccharomyces cerevisiae S288C (brewer's yeast)",
 'SpeciesTaxid': 4932,
 'SpeciesName': 'Saccharomyces cerevisiae',
 'AssemblyType': 'haploid',
 'AssemblyStatus': 'Complete Genome',
 'AssemblyStatusSort': '1',
 'WGS': 'None',
 'BioSampleAccn': 'None',
 'BioSampleId': 'None',
 'Primary': '285488',
 'ReleaseLevel': 'Major',
 'ReleaseType': 'Major',
 'AsmReleaseDate_GenBank': '2011/05/27 00:00',
 'AsmReleaseDate_RefSeq': '2011/05/27 00:00',
 'SeqReleaseDate': '2014/12/17 00:00',
 'AsmUpdateDate': '2016/11/29 00:00',
 'SubmissionDate': '2014/12/17 00:00',
 'LastUpdateDate': '2016/11/29 00:00',
 'SubmitterOrganization': 'Saccharomyces Genome Database',
 'RefSeq_category': 'reference genome',
 'PropertyList': ['full-genome-representation',
  'genbank_has_annotation',
  'has-chromosome',
  'has-mitochondrion',
  'has-synonym',
  'has_annotation',
  'latest',
  'latest_genbank',
  'latest_refseq',
  'reference',
  'refseq_has_annotation'],
 'FromType': 'None',
 'Synonym': {'Genbank': 'GCA_000146045.2',
  'RefSeq': 'GCF_000146045.2',
  'Similarity': 'different'},
 'ContigN50': '924431',
 'ScaffoldN50': '924431',
 'FtpPath_GenBank': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/146/045/GCA_000146045.2_R64',
 'FtpPath_RefSeq': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64',
 'FtpPath_Assembly_rpt': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_assembly_report.txt',
 'FtpPath_Stats_rpt': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_assembly_stats.txt',
 'FtpPath_Regions_rpt': 'None'}
    """

    assembly = NCBIAssembly(assembly_data, cache='/storage/fliu/data/biodb/ncbi/')
    print(assembly.cwd_ftp_path_gb, assembly.cwd_ftp_path_rs)
    genome_contigs = assembly.get_genomic_fna()
    genome_cds = assembly.get_protein_faa()

    from ftplib import FTP
    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    with (ftp):
        assembly.fetch_ncbi_ftp_data(ftp)
        # assembly.fetch_ncbi_ftp_data(ftp)