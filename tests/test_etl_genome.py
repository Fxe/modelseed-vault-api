"""
Tests for the NCBI genome ETL pipeline using E. coli K-12 MG1655 (GCF_000005845.2).
"""

from pathlib import Path

import pytest

from modelseed_vault.elt.extract.extract_ncbi_genome import ExtractGenomeNCBI
from modelseed_vault.elt.transform.genome.transform_genome import TransformGenome

DATA_DIR = Path(__file__).parent / 'data'
FNA = DATA_DIR / 'GCF_000005845.2_ASM584v2_genomic.fna.gz'
FAA = DATA_DIR / 'GCF_000005845.2_ASM584v2_protein.faa.gz'
GFF = DATA_DIR / 'GCF_000005845.2_ASM584v2_genomic.gff.gz'
GENOME_ID = 'GCF_000005845.2'


@pytest.fixture(scope='module')
def container():
    return ExtractGenomeNCBI().extract(FNA, FAA, GFF)


@pytest.fixture(scope='module')
def graph(container):
    return TransformGenome().transform(GENOME_ID, container)


# ── Extract tests ─────────────────────────────────────────────────────────────

def test_extract_contigs(container):
    assert container.contigs is not None
    contig_list = list(container.contigs.features)
    assert len(contig_list) >= 1
    accession = contig_list[0].id.split(' ')[0]
    assert accession == 'NC_000913.3'


def test_extract_proteins(container):
    assert container.proteins is not None
    protein_list = list(container.proteins.features)
    assert len(protein_list) > 100


def test_extract_features(container):
    assert container.features is not None
    assert len(container.features) > 100
    feat = container.features[0]
    assert 'contig_id' in feat
    assert 'protein_id' in feat
    assert 'locus_tag' in feat
    assert 'start' in feat and 'end' in feat and 'strand' in feat


def test_extract_feature_protein_ids_match(container):
    protein_ids = {f.id.split(' ')[0] for f in container.proteins.features}
    feature_protein_ids = {f['protein_id'] for f in container.features}
    overlap = protein_ids & feature_protein_ids
    assert len(overlap) > 100, f'Expected significant overlap, got {len(overlap)}'


# ── Transform tests ───────────────────────────────────────────────────────────

def test_transform_genome_node(graph):
    assert 'Genome' in graph.t_nodes
    assert f'Genome/{GENOME_ID}' in graph.t_nodes['Genome']


def test_transform_contig_nodes(graph):
    assert 'GenomicContig' in graph.t_nodes
    assert len(graph.t_nodes['GenomicContig']) >= 1


def test_transform_feature_nodes(graph):
    assert 'GenomicFeature' in graph.t_nodes
    assert len(graph.t_nodes['GenomicFeature']) > 100


def test_transform_protein_sequence_nodes(graph):
    assert 'ProteinSequence' in graph.t_nodes
    assert len(graph.t_nodes['ProteinSequence']) > 100


def test_transform_dna_sequence_nodes(graph):
    assert 'DNASequence' in graph.t_nodes
    assert len(graph.t_nodes['DNASequence']) > 100


def test_transform_has_contig_edges(graph):
    assert 'has_contig' in graph.t_edges
    assert len(graph.t_edges['has_contig']) >= 1


def test_transform_has_feature_edges(graph):
    assert 'has_feature' in graph.t_edges
    assert len(graph.t_edges['has_feature']) > 100


def test_transform_has_protein_sequence_edges(graph):
    assert 'has_protein_sequence' in graph.t_edges
    assert len(graph.t_edges['has_protein_sequence']) > 100


def test_transform_has_dna_sequence_edges(graph):
    assert 'has_dna_sequence' in graph.t_edges
    assert len(graph.t_edges['has_dna_sequence']) > 100


def test_transform_hash_node_keys_are_hex(graph):
    for node in graph.t_nodes['GenomicContig'].values():
        assert len(node.key) == 64  # SHA-256 hex digest
    for node in graph.t_nodes['ProteinSequence'].values():
        assert len(node.key) == 64
    for node in graph.t_nodes['DNASequence'].values():
        assert len(node.key) == 64


def test_transform_summary(graph, capsys):
    graph.summary()
    captured = capsys.readouterr()
    assert 'Genome' in captured.out
    assert 'GenomicContig' in captured.out
