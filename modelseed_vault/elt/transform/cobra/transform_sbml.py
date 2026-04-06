from modelseed_vault.core.transform_graph import TransformGraph, Node
from modelseed_vault.core.graph_ontology import (
    SBML_MODEL, SBML_COMPARTMENT, SBML_SPECIES, SBML_GENE, SBML_REACTION, SBML_COMPLEX,
    HAS_SBML_COMPARTMENT, HAS_SBML_SPECIES, HAS_SBML_REACTION, HAS_SBML_GENE,
    HAS_REACTANT, HAS_PRODUCT, HAS_SBML_COMPLEX,
)
from modelseed_vault.utils import sha_hex
from modelseed_vault.elt.transform.cobra.parse import scan_gpr_nodes


def _clean(d: dict) -> dict:
    result = {k: v for k, v in d.items() if not k.startswith('_')}
    if '_raw_xml' in d:
        result['_raw_xml'] = d['_raw_xml']
    return result


def _build_graph(model_id: str, container) -> TransformGraph:
    """
    Build a TransformGraph from an ExtractContainerCobra.

    Nodes:  SBMLModel, SBMLCompartment, SBMLSpecies, SBMLReaction,
            SBMLComplex, SBMLGene
    Edges:  has_sbml_compartment, has_sbml_species, has_sbml_reaction,
            has_sbml_gene, has_reactant, has_product, has_sbml_complex
    """
    g = TransformGraph()

    # ── SBMLModel ────────────────────────────────────────────────────────────
    node_model = g.add_transform_node(Node(model_id, SBML_MODEL, data=_clean(container.model)))

    # ── SBMLCompartment ──────────────────────────────────────────────────────
    d_compartment = {}
    for elem in container.compartments['elements']:
        node = g.add_transform_node(
            Node(f"{model_id}:{elem['id']}", SBML_COMPARTMENT, data=_clean(elem))
        )
        d_compartment[elem['id']] = node
        g.add_transform_edge(node_model, node, HAS_SBML_COMPARTMENT)

    # ── SBMLSpecies ──────────────────────────────────────────────────────────
    d_species = {}
    for elem in container.species['elements']:
        node = g.add_transform_node(
            Node(f"{model_id}:{elem['id']}", SBML_SPECIES, data=_clean(elem))
        )
        d_species[elem['id']] = node
        g.add_transform_edge(node_model, node, HAS_SBML_SPECIES)
        comp_id = elem.get('compartment')
        if comp_id and comp_id in d_compartment:
            g.add_transform_edge(node, d_compartment[comp_id], HAS_SBML_COMPARTMENT)

    # ── SBMLGene ─────────────────────────────────────────────────────────────
    d_gene = {}
    for gp in container.fbc_gene_products:
        gene_id = gp.get('fbc:id') or gp.get('id', '')
        if not gene_id:
            continue
        node_gene = g.add_transform_node(
            Node(f"{model_id}:{gene_id}", SBML_GENE, data=gp)
        )
        d_gene[gene_id] = node_gene
        g.add_transform_edge(node_model, node_gene, HAS_SBML_GENE)

    # ── SBMLReaction + stoichiometry + gene complexes ─────────────────────
    for elem in container.reactions['elements']:
        rxn_id = elem['id']
        node_rxn = g.add_transform_node(
            Node(f"{model_id}:{rxn_id}", SBML_REACTION, data=_clean(elem))
        )
        g.add_transform_edge(node_model, node_rxn, HAS_SBML_REACTION)

        for sid, stoich in elem.get('_reactants', []):
            if sid in d_species:
                data = {'stoichiometry': stoich} if stoich is not None else {}
                g.add_transform_edge(node_rxn, d_species[sid], HAS_REACTANT, data)

        for sid, stoich in elem.get('_products', []):
            if sid in d_species:
                data = {'stoichiometry': stoich} if stoich is not None else {}
                g.add_transform_edge(node_rxn, d_species[sid], HAS_PRODUCT, data)

        complexes = elem.get('_complexes') or scan_gpr_nodes(elem.get('_raw_xml', ''))
        for complex_genes in complexes:
            genes = sorted(set(complex_genes))
            if not genes:
                continue
            complex_key = sha_hex("|".join(genes))
            node_complex = g.add_transform_node(
                Node(f"{model_id}:{complex_key}", SBML_COMPLEX, data={'genes': genes})
            )
            g.add_transform_edge(node_rxn, node_complex, HAS_SBML_COMPLEX)
            for gene_id in genes:
                if gene_id not in d_gene:
                    node_gene = g.add_transform_node(
                        Node(f"{model_id}:{gene_id}", SBML_GENE, data={'id': gene_id})
                    )
                    d_gene[gene_id] = node_gene
                    g.add_transform_edge(node_model, node_gene, HAS_SBML_GENE)
                g.add_transform_edge(node_complex, d_gene[gene_id], HAS_SBML_GENE)

    return g


class TransformCobraSBML:

    def transform(self, model_id: str, container) -> TransformGraph:
        return _build_graph(model_id, container)
