from modelseed_vault.core.graph_ontology import GENOME_SET, REFSEQ_GENOME, KBASE_GENE, KEGG_ORTHOLOGY


class ObjectFactory:

    def __init__(self):
        pass

    def build(self, node):
        labels = node.labels
        if GENOME_SET in labels:
            pass
        elif REFSEQ_GENOME in labels:
            pass
        elif KBASE_GENE in labels:
            pass
        elif KEGG_ORTHOLOGY in labels:
            pass
