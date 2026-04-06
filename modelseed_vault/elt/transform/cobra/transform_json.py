from modelseed_vault.core.transform_graph import TransformGraph, Node
from modelseed_vault.core.graph_ontology import (
    COBRA_MODEL, COBRA_COMPARTMENT, COBRA_METABOLITE, COBRA_GENE, COBRA_REACTION, COBRA_COMPLEX,
    HAS_COMPARTMENT, HAS_SPECIES, HAS_REACTION, HAS_GENE,
    HAS_STOICHIOMETRY_COEFFICIENT, HAS_GENE_COMPLEX,
)
from modelseed_vault.utils import sha_hex


def _clean(d: dict) -> dict:
    return {k: v for k, v in d.items() if not k.startswith('_')}


class TransformCobraJson:

    def transform(self, model_id: str, container) -> TransformGraph:
        """
        Build a TransformGraph from an ExtractContainerCobra (JSON source).

        Uses COBRA-specific node labels to preserve provenance:

        Nodes:  COBRAModel, COBRACompartment, COBRAMetabolite, COBRAReaction,
                COBRAGene, COBRAComplex
        Edges:  has_compartment, has_species, has_reaction, has_gene,
                has_stoichiometry_coefficient, has_gene_complex, has_gene (complex→gene)
        """
        g = TransformGraph()

        # ── COBRAModel ───────────────────────────────────────────────────────
        node_model = g.add_transform_node(
            Node(model_id, COBRA_MODEL, data=_clean(container.model))
        )

        # ── COBRACompartment ─────────────────────────────────────────────────
        d_compartment = {}
        for elem in container.compartments['elements']:
            node = g.add_transform_node(
                Node(f"{model_id}:{elem['id']}", COBRA_COMPARTMENT, data=_clean(elem))
            )
            d_compartment[elem['id']] = node
            g.add_transform_edge(node_model, node, HAS_COMPARTMENT)

        # ── COBRAMetabolite ──────────────────────────────────────────────────
        d_metabolite = {}
        for elem in container.species['elements']:
            node = g.add_transform_node(
                Node(f"{model_id}:{elem['id']}", COBRA_METABOLITE, data=_clean(elem))
            )
            d_metabolite[elem['id']] = node
            g.add_transform_edge(node_model, node, HAS_SPECIES)
            comp_id = elem.get('compartment')
            if comp_id and comp_id in d_compartment:
                g.add_transform_edge(node, d_compartment[comp_id], HAS_COMPARTMENT)

        # ── COBRAGene ────────────────────────────────────────────────────────
        d_gene = {}
        for gp in container.fbc_gene_products:
            gene_id = gp.get('fbc:id') or gp.get('id', '')
            if not gene_id:
                continue
            node_gene = g.add_transform_node(
                Node(f"{model_id}:{gene_id}", COBRA_GENE, data=gp)
            )
            d_gene[gene_id] = node_gene
            g.add_transform_edge(node_model, node_gene, HAS_GENE)

        # ── COBRAReaction + stoichiometry + gene complexes ───────────────────
        for elem in container.reactions['elements']:
            rxn_id = elem['id']
            node_rxn = g.add_transform_node(
                Node(f"{model_id}:{rxn_id}", COBRA_REACTION, data=_clean(elem))
            )
            g.add_transform_edge(node_model, node_rxn, HAS_REACTION)

            # Reactants: coefficient is negative
            for sid, stoich in elem.get('_reactants', []):
                if sid in d_metabolite:
                    g.add_transform_edge(
                        node_rxn, d_metabolite[sid],
                        HAS_STOICHIOMETRY_COEFFICIENT,
                        {'coefficient': -float(stoich)},
                    )

            # Products: coefficient is positive
            for sid, stoich in elem.get('_products', []):
                if sid in d_metabolite:
                    g.add_transform_edge(
                        node_rxn, d_metabolite[sid],
                        HAS_STOICHIOMETRY_COEFFICIENT,
                        {'coefficient': float(stoich)},
                    )

            # Gene complexes (AND groups)
            for complex_genes in elem.get('_complexes', []):
                genes = sorted(set(complex_genes))
                if not genes:
                    continue
                complex_key = sha_hex("|".join(genes))
                node_complex = g.add_transform_node(
                    Node(f"{model_id}:{complex_key}", COBRA_COMPLEX, data={'genes': genes})
                )
                g.add_transform_edge(node_rxn, node_complex, HAS_GENE_COMPLEX)
                for gene_id in genes:
                    if gene_id in d_gene:
                        g.add_transform_edge(node_complex, d_gene[gene_id], HAS_GENE)

        return g
