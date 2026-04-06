# ── Node Labels ──────────────────────────────────────────────────────────────

# SBML
SBML_MODEL = 'SBMLModel'
SBML_COMPARTMENT = 'SBMLCompartment'
SBML_SPECIES = 'SBMLSpecies'
SBML_GENE = 'SBMLGene'
SBML_REACTION = 'SBMLReaction'
SBML_COMPLEX = 'SBMLComplex'

# COBRA JSON
COBRA_MODEL = 'COBRAModel'
COBRA_COMPARTMENT = 'COBRACompartment'
COBRA_METABOLITE = 'COBRAMetabolite'
COBRA_GENE = 'COBRAGene'
COBRA_REACTION = 'COBRAReaction'
COBRA_COMPLEX = 'COBRAComplex'

# Genome
GENOME = 'Genome'
GENOMIC_CONTIG = 'GenomicContig'
GENOMIC_FEATURE = 'GenomicFeature'
LOCUS_TAG = 'LocusTag'
LOCUS_TAG_OLD = 'OldLocusTag'
PROTEIN_SEQUENCE = 'ProteinSequence'
DNA_SEQUENCE = 'DNASequence'

# UniProt (static)
UNIPROTKB_ACCESSION = 'uniprotkb_accession'
UNIPROTKB_SUBCELL = 'uniprotkb_subcell'
ECO_TERM = 'eco_term'
RHEA_REACTION = 'RHEA_reaction'
CHEBI_TERM = 'ChEBI_term'
EC_NUMBER = 'EC_number'
RE_SEQ_PROTEIN = 're_seq_protein'
KEGG_GENE = 'kegg_gene'
ALPHAFOLDDB = 'alphafolddb'

# DAO / graph nodes
PROTEIN_ANNOTATION = 'ProteinAnnotation'
RAST = 'RAST'
METHOD = 'Method'
GENOME_SET = 'GenomeSet'
REFSEQ_GENOME = 'RefSeqGenome'
KBASE_GENE = 'KBaseGene'
KEGG_ORTHOLOGY = 'KeggOrthology'


def uniprotkb_collection(uniprot_type: str) -> str:
    return f"uniprotkb_{uniprot_type}"


# ── Edge / Relationship Types ─────────────────────────────────────────────────

# SBML
HAS_SBML_COMPARTMENT = 'has_compartment'
HAS_SBML_SPECIES = 'has_species'
HAS_SBML_REACTION = 'has_reaction'
HAS_SBML_GENE = 'has_gene'
HAS_REACTANT = 'has_reactant'
HAS_PRODUCT = 'has_product'
HAS_SBML_COMPLEX = 'has_gpr_complex'

# COBRA JSON
HAS_COMPARTMENT = 'has_compartment'
HAS_SPECIES = 'has_species'
HAS_REACTION = 'has_reaction'
HAS_GENE = 'has_gene'
HAS_STOICHIOMETRY_COEFFICIENT = 'has_stoichiometry_coefficient'
HAS_GENE_COMPLEX = 'has_gpr_complex'

# Genome
HAS_CONTIG = 'has_contig'
HAS_FEATURE = 'has_feature'
HAS_LOCUS_TAG = 'has_locus_tag'
HAS_PROTEIN_SEQUENCE = 'has_protein_sequence'
HAS_DNA_SEQUENCE = 'has_dna_sequence'
HAS_PROTEIN_TRANSLATION = 'has_protein_translation'

PARENT = 'parent'
CHILD = 'child'

# UniProt (dynamic — parameterized by uniprot_type)
def uniprotkb_has_ec(uniprot_type: str) -> str:
    return f"uniprotkb_{uniprot_type}_has_ec"


def uniprotkb_has_subcell(uniprot_type: str) -> str:
    return f"uniprotkb_{uniprot_type}_has_subcell"


def uniprotkb_has_accession(uniprot_type: str) -> str:
    return f"uniprotkb_{uniprot_type}_has_accession"


def uniprotkb_has_protein_sequence(uniprot_type: str) -> str:
    return f"uniprotkb_{uniprot_type}_has_protein_sequence"


def uniprotkb_has_reference_to_kegg_gene(uniprot_type: str) -> str:
    return f"uniprotkb_{uniprot_type}_has_reference_to_kegg_gene"


def uniprotkb_has_reference_to_alphafolddb(uniprot_type: str) -> str:
    return f"uniprotkb_{uniprot_type}_has_reference_to_alphafolddb"


def uniprotkb_has_cofactor_chebi_term(uniprot_type: str) -> str:
    return f"uniprotkb_{uniprot_type}_has_cofactor_chebi_term"


def uniprotkb_has_catalytic_activity_rhea_reaction(uniprot_type: str) -> str:
    return f"uniprotkb_{uniprot_type}_has_catalytic_activity_rhea_reaction"


def uniprotkb_has_catalytic_activity_ec_number(uniprot_type: str) -> str:
    return f"uniprotkb_{uniprot_type}_has_catalytic_activity_ec_number"


# Client / DAO
HAS_ANNOTATION = 'has_annotation'
HAS_GENOME = 'has_genome'
