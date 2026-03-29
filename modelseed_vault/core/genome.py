from modelseed_vault.dao_neo4j import Neo4jDAO
from modelseedpy.core.msgenome import MSGenome, MSFeature

class GenomicsProtein:

    def __init__(self, protein_sequence: str):
        self.protein_sequence = protein_sequence


class GenomicsFeature(MSFeature):

    def __init__(self, feature_id: str, sequence: str,
                 start: int, end: int, strand: str, feature_type: str,
                 protein: GenomicsProtein):

        super().__init__(feature_id, sequence)
        self.start = start
        self.end = end
        self.strand = strand
        self.feature_type = feature_type
        self.protein = protein

    def get_ontology_term(self, ontology_term):
        pass

    def remove_ontology_term(self, ontology_term, value):
        pass

    def add_ontology_term(self, ontology_term, value):
        pass


class GenomicsGenome:

    def __init__(self, genome_sequence: str, neo4j_dao: Neo4jDAO):
        self.neo4j_dao = neo4j_dao
        self.features = None
        pass

    def get_feature(self, feature_id: str) -> GenomicsFeature:
        if self.features is None:
            self.features = {}
        if feature_id not in self.features:
            self.features[feature_id] = GenomicsFeature(feature_id)
        return self.features[feature_id]
    
    def load_features(self):
        pass

