

class ObjectFactory:

    def __init__(self):
        pass

    def build(self, node):
        labels = node.labels
        if 'GenomeSet' in labels:
            pass
        elif 'RefSeqGenome' in labels:
            pass
        elif 'KBaseGene' in labels:
            pass
        elif 'KeggOrthology' in labels:
            pass
