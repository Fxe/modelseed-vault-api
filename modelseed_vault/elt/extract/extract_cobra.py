from modelseed_vault.elt.transform.cobra.parse import (
    parse_elements_with_provenance,
    parse_reactions_with_provenance,
    parse_model_tag,
    parse_parameters,
    parse_fbc_objectives,
    parse_fbc_gene_products,
    parse_unit_definitions,
    parse_groups,
    parse_gene_associations,
)
from pathlib import Path
import hashlib
import json
import re


def _split_or(expr: str) -> list:
    """Split a gene reaction rule expression by ' or ' at the top paren level."""
    parts, current, depth = [], [], 0
    i = 0
    while i < len(expr):
        if expr[i] == '(':
            depth += 1
            current.append(expr[i])
        elif expr[i] == ')':
            depth -= 1
            current.append(expr[i])
        elif depth == 0 and expr[i:i + 4] == ' or ':
            parts.append(''.join(current))
            current = []
            i += 4
            continue
        else:
            current.append(expr[i])
        i += 1
    if current:
        parts.append(''.join(current))
    return parts


def _parse_grr(rule: str) -> list:
    """
    Parse a COBRA gene_reaction_rule string into a list of complexes.
    Each complex is a list of gene ID strings (AND group).
    OR relationships produce separate complexes.
    """
    rule = rule.strip()
    if not rule:
        return []
    complexes = []
    for grp in _split_or(rule):
        grp = grp.strip().strip('()')
        genes = [g.strip() for g in re.split(r'\s+and\s+', grp) if g.strip()]
        if genes:
            complexes.append(genes)
    return complexes


class ExtractContainerCobra:

    def __init__(self,
                 model: dict,
                 compartments: dict,
                 species: dict,
                 reactions: dict,
                 parameters: list,
                 fbc_objectives: list,
                 fbc_gene_products: list,
                 unit_definitions: list,
                 groups: list):
        self.model = model
        self.compartments = compartments
        self.species = species
        self.reactions = reactions
        self.parameters = parameters
        self.fbc_objectives = fbc_objectives
        self.fbc_gene_products = fbc_gene_products
        self.unit_definitions = unit_definitions
        self.groups = groups


class ExtractCobraJson:

    def __init__(self):
        pass

    def extract(self, json_file: Path) -> ExtractContainerCobra:
        raw_bytes = Path(json_file).read_bytes()
        file_sha256 = hashlib.sha256(raw_bytes).hexdigest()
        d = json.loads(raw_bytes)

        model = {
            'id': d.get('id'),
            'name': d.get('name'),
            'version': d.get('version'),
        }

        compartments = {
            'file_sha256': file_sha256,
            'elements': [{'id': cid, 'name': cname} for cid, cname in d.get('compartments', {}).items()],
        }

        species = {
            'file_sha256': file_sha256,
            'elements': list(d.get('metabolites', [])),
        }

        reaction_elements = []
        for rxn in d.get('reactions', []):
            stoich: dict = rxn.get('metabolites', {})
            reactants = [[sid, str(abs(v))] for sid, v in stoich.items() if v < 0]
            products  = [[sid, str(v)]      for sid, v in stoich.items() if v > 0]
            rec = {k: v for k, v in rxn.items() if k != 'metabolites'}
            rec['_reactants'] = reactants
            rec['_products']  = products
            rec['_complexes'] = _parse_grr(rxn.get('gene_reaction_rule', ''))
            reaction_elements.append(rec)

        reactions = {
            'file_sha256': file_sha256,
            'elements': reaction_elements,
        }

        fbc_gene_products = [
            {'fbc:id': g['id'], 'fbc:label': g.get('name', g['id']), **g.get('annotation', {})}
            for g in d.get('genes', [])
        ]

        return ExtractContainerCobra(
            model=model,
            compartments=compartments,
            species=species,
            reactions=reactions,
            parameters=[],
            fbc_objectives=[],
            fbc_gene_products=fbc_gene_products,
            unit_definitions=[],
            groups=[],
        )


class ExtractCobraSBML:

    def __init__(self):
        pass

    def extract(self, sbml_file: Path) -> ExtractContainerCobra:
        p = str(sbml_file)
        rxn_xpath = "//*[local-name()='model']/*[local-name()='listOfReactions']/*[local-name()='reaction']"

        with open(p, 'rb') as fh:
            model = parse_model_tag(fh)
        with open(p, 'rb') as fh:
            compartments = parse_elements_with_provenance(fh, 'compartment')
        with open(p, 'rb') as fh:
            species = parse_elements_with_provenance(fh, 'species')
        with open(p, 'rb') as fh:
            reactions = parse_reactions_with_provenance(fh, xpath=rxn_xpath)
        with open(p, 'rb') as fh:
            gene_associations = parse_gene_associations(fh)
        for elem in reactions['elements']:
            elem['_complexes'] = gene_associations.get(elem.get('id', ''), [])
        with open(p, 'rb') as fh:
            parameters = parse_parameters(fh)
        with open(p, 'rb') as fh:
            fbc_objectives = parse_fbc_objectives(fh)
        with open(p, 'rb') as fh:
            fbc_gene_products = parse_fbc_gene_products(fh)
        with open(p, 'rb') as fh:
            unit_definitions = parse_unit_definitions(fh)
        with open(p, 'rb') as fh:
            groups = parse_groups(fh)

        return ExtractContainerCobra(
            model=model,
            compartments=compartments,
            species=species,
            reactions=reactions,
            parameters=parameters,
            fbc_objectives=fbc_objectives,
            fbc_gene_products=fbc_gene_products,
            unit_definitions=unit_definitions,
            groups=groups,
        )
