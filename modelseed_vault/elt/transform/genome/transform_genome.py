import logging

from modelseed_vault.core.transform_graph import TransformGraph, Node, HashNode
from modelseed_vault.core.graph_ontology import (
    GENOME, GENOMIC_CONTIG, GENOMIC_FEATURE, LOCUS_TAG, LOCUS_TAG_OLD,
    PROTEIN_SEQUENCE, DNA_SEQUENCE,
    HAS_CONTIG, HAS_FEATURE, HAS_LOCUS_TAG,
    HAS_PROTEIN_SEQUENCE, HAS_DNA_SEQUENCE, HAS_PROTEIN_TRANSLATION, PARENT
)
from modelseed_vault.elt.extract.extract_ncbi_genome import ExtractContainerGenome

logger = logging.getLogger(__name__)

_COMPLEMENT = str.maketrans('ACGTacgt', 'TGCAtgca')


def _reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


class TransformGenome:

    def transform(self, genome_id: str, container: ExtractContainerGenome) -> TransformGraph:
        g = TransformGraph()
        unmatched = []

        # Genome node
        node_genome = g.add_transform_node(Node(genome_id, GENOME))

        # contig accession → MSFeature (seq = full nucleotide string)
        contig_lookup = {f.id.split(' ')[0]: f for f in container.contigs.features}

        # protein_id → MSFeature
        protein_lookup = {f.id.split(' ')[0]: f for f in container.proteins.features}

        # GenomicContig HashNodes — one per contig sequence
        contig_nodes = {}
        for accession, feat in contig_lookup.items():
            seq_str = str(feat.seq)
            node_contig = g.add_transform_node(
                HashNode(seq_str, GENOMIC_CONTIG, data={'accession': accession, 'length': len(seq_str)})
            )
            contig_nodes[accession] = node_contig
            g.add_transform_edge(node_genome, node_contig, HAS_CONTIG)

        # locus_tag → DNASequence node (populated when processing gene features)
        locus_tag_dna = {}
        # locus_tag → ProteinSequence node (populated when processing CDS features)
        locus_tag_protein = {}

        # Process all GFF features
        for feat in container.features:
            contig_id = feat['contig_id']
            feature_type = feat['feature_type']
            attrs = feat['attrs']

            # Use GFF positional fields as the unique identifier — locus_tag is
            # not unique across feature types (e.g. gene + CDS share the same tag)
            node_id = f"{genome_id}_{contig_id}_{feature_type}_{feat['start']}_{feat['end']}"

            # All 9 GFF columns + all attributes — enough to reconstruct the record
            node_data = {
                'contig_id': contig_id,
                'source': feat['source'],
                'feature_type': feature_type,
                'start': feat['start'],
                'end': feat['end'],
                'score': feat['score'],
                'strand': feat['strand'],
                'phase': feat['phase'],
            }
            node_data.update({f'gff_attr_{k}': v for k, v in attrs.items()})

            node_feature = g.add_transform_node(
                Node(node_id, GENOMIC_FEATURE, data=node_data)
            )

            # Connect feature to its contig
            contig_node = contig_nodes.get(contig_id)
            if contig_node:
                g.add_transform_edge(contig_node, node_feature, HAS_FEATURE)
            else:
                unmatched.append({'reason': 'contig_not_found', 'contig_id': contig_id, 'feature': feat})
                logger.warning('unmatched feature: contig %s not found (node_id=%s)', contig_id, node_id)

            # LocusTag node — one per unique locus_tag; all features sharing it get HAS_LOCUS_TAG
            locus_tag = attrs.get('locus_tag')
            node_locus = None
            if locus_tag:
                node_locus = g.add_transform_node(Node(locus_tag, LOCUS_TAG))
                g.add_transform_edge(node_feature, node_locus, HAS_LOCUS_TAG)
            locus_tag_old = attrs.get('old_locus_tag')
            if locus_tag_old:
                node_locus_old = g.add_transform_node(Node(locus_tag_old, LOCUS_TAG, [LOCUS_TAG_OLD]))
                if node_locus is not None:
                    g.add_transform_edge(node_locus, node_locus_old, PARENT)
                else:
                    g.add_transform_edge(node_feature, node_locus_old, HAS_LOCUS_TAG)

            # gene: extract DNASequence from the contig at the feature's coordinates
            if feature_type == 'gene':
                contig_feat = contig_lookup.get(contig_id)
                if contig_feat:
                    full_seq = str(contig_feat.seq)
                    dna_seq = full_seq[feat['start'] - 1: feat['end']]  # GFF is 1-based inclusive
                    if feat['strand'] == '-':
                        dna_seq = _reverse_complement(dna_seq)
                    if dna_seq:
                        node_dna = g.add_transform_node(
                            HashNode(dna_seq, DNA_SEQUENCE, data={'length': len(dna_seq)})
                        )
                        g.add_transform_edge(node_feature, node_dna, HAS_DNA_SEQUENCE)
                        if locus_tag:
                            locus_tag_dna[locus_tag] = node_dna
                else:
                    unmatched.append({'reason': 'contig_not_found_for_dna', 'contig_id': contig_id, 'feature': feat})
                    logger.warning('cannot extract DNASequence for gene: contig %s not found (node_id=%s)', contig_id, node_id)

            # CDS: connect ProteinSequence
            elif feature_type == 'CDS':
                protein_id = attrs.get('protein_id')
                prot_feat = protein_lookup.get(protein_id) if protein_id else None
                if prot_feat:
                    prot_seq = str(prot_feat.seq)
                    node_protein = g.add_transform_node(
                        HashNode(prot_seq, PROTEIN_SEQUENCE, data={'protein_id': protein_id, 'length': len(prot_seq)})
                    )
                    g.add_transform_edge(node_feature, node_protein, HAS_PROTEIN_SEQUENCE)
                    if locus_tag:
                        locus_tag_protein[locus_tag] = node_protein
                else:
                    unmatched.append({'reason': 'protein_not_found', 'protein_id': protein_id, 'feature': feat})
                    logger.warning('unmatched CDS: protein_id %s not found (node_id=%s)', protein_id, node_id)

        # Link DNASequence → ProteinSequence for matched locus_tags
        for locus_tag, node_dna in locus_tag_dna.items():
            node_protein = locus_tag_protein.get(locus_tag)
            if node_protein:
                g.add_transform_edge(node_dna, node_protein, HAS_PROTEIN_TRANSLATION)
            else:
                logger.warning('gene locus_tag %s has no matching CDS — DNASequence has no protein translation', locus_tag)

        for locus_tag in locus_tag_protein:
            if locus_tag not in locus_tag_dna:
                logger.warning('CDS locus_tag %s has no matching gene — ProteinSequence has no DNA source', locus_tag)

        g.unmatched = unmatched
        return g
