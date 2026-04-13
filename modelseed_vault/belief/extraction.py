

PROFILE_COBRA_V0 = {
    'nodes': {
        'GenomicFeature': [True, True],
        'LocusTag': [True, True],
        'SBMLGene': [True, True],
        'SBMLComplex': [True, True],
        'SBMLReaction': [True, True],
        'SBMLSpecies': [True, False],
        'COBRAComplex': [True, True],
        'COBRAGene': [True, True],
        'COBRAReaction': [True, True],
        'COBRAMetabolite': [True, False]
        # 'SBMLSpecies': [True, False],}
    },
    'attributes': {
        'GenomicFeature': {
        #'source': 'function_source',
        'source': 'gff_source',
        'gff_attr_product': 'function',
        'gff_attr_locus_tag': 'locus_tag',
        'gff_attr_old_locus_tag': 'old_locus_tag',
        'feature_type': None
    },
    'COBRAReaction': {
        'name': None,
        'lower_bound': None,
        'upper_bound': None,
    },
    'COBRAMetabolite': {
        'name': None,
        'compartment': None,
        'formula': None,
    },
    'SBMLReaction': {
        'name': None,
        'id': None
    },
    'SBMLSpecies': {
        'name': None,
        'id': None,
        'compartment': None,
    },
    'ModelSEEDCompound': {
        'name': None,
        'formula': None,
    },
    'ModelSEEDReaction': {
        'name': None,
        'definition': None,
        'deltag': None,
        'deltag_eQuilibrator': None,
    }}
}


class ExtractionProfile:

    def __init__(self):
        pass
