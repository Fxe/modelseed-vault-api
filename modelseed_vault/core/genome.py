from modelseed_vault.dao_neo4j import Neo4jDAO
from modelseed_vault.core.hash import HashString
from modelseedpy.core.msgenome import MSGenome, MSFeature


class ProteinSequence(HashString):

    _STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
    _EXTENDED_ONLY = {"U", "O"}  # selenocysteine, pyrrolysine
    _AMBIGUOUS_AA = set("BJZX*")

    def __new__(cls, sequence, strip_ending_star=True):
        obj = super().__new__(cls, sequence, strip_ending_star=strip_ending_star)
        return obj

    @property
    def sequence(self) -> str:
        return str(self)

    def is_standard(self) -> bool:
        """True if the sequence contains only standard amino acids."""
        return bool(self) and set(self) <= self._STANDARD_AA

    def is_extended(self) -> bool:
        """True if all residues are standard or extended (U, O)."""
        allowed = self._STANDARD_AA | self._EXTENDED_ONLY
        return bool(self) and set(self) <= allowed

    def is_ambiguous(self) -> bool:
        """True if the sequence contains any ambiguous residue codes."""
        return any(res in self._AMBIGUOUS_AA for res in self)

    def is_valid(self) -> bool:
        """True if all residues are recognized AA codes (standard, extended, or ambiguous)."""
        allowed = self._STANDARD_AA | self._EXTENDED_ONLY | self._AMBIGUOUS_AA
        return bool(self) and set(self) <= allowed

    def z_compress(self):
        import zlib
        return zlib.compress(str(self).encode("utf-8"))


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

