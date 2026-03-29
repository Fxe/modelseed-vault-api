from modelseed_vault.core.node import Node
from modelseed_vault.elt.transform.cobra.parse import parse_species_list
from modelseed_vault.core.transform_graph import TransformGraph


class TransformCobraSBML:

    def __init__(self, model_id: str):
        self.model_id = model_id

        """
        d_node_compartment = {}
        for n in l_node_compartment:
            sid = n.data['id']
            if sid not in d_node_compartment:
                d_node_compartment[sid] = n
            else:
                raise ValueError('@')
        d_node_species = {}
        for n in l_node_species:
            sid = n.data['id']
            if sid not in d_node_species:
                d_node_species[sid] = n
            else:
                raise ValueError('@')
        """

    def transform_edges(self, node_model, nodes, reactions, wrapper):
        d_node_species = {n.data['id']: n for n in nodes if n.label == 'SBMLSpecies'}
        d_node_reaction = {n.data['id']: n for n in nodes if n.label == 'SBMLReaction'}
        d_node_compartment = {n.data['id']: n for n in nodes if n.label == 'SBMLCompartment'}

        res = {}
        for d_reaction in reactions['elements']:
            str_xml = d_reaction['_raw_xml']
            node_reaction = d_node_reaction[d_reaction['id']]
            edges = self.build_stoich_links(node_reaction.key, d_node_species, str_xml, wrapper)
            for k, v in edges.items():
                if k not in res:
                    res[k] = v
                else:
                    raise ValueError("dup edge")

        for node in d_node_compartment.values():
            t = (('SBMLModel', node_model.key), 'has_sbml_compartment', ('SBMLCompartment', node.key))
            if t not in res:
                res[t] = None
            else:
                raise ValueError("dup edge")

        for node in d_node_species.values():
            t = (('SBMLModel', node_model.key), 'has_sbml_species', ('SBMLSpecies', node.key))
            if t not in res:
                res[t] = None
            else:
                raise ValueError("dup edge")

        for node in d_node_reaction.values():
            t = (('SBMLModel', node_model.key), 'has_sbml_reaction', ('SBMLReaction', node.key))
            if t not in res:
                res[t] = None
            else:
                raise ValueError("dup edge")

        return res

    def transform(self, model_id: str, model_metadata: dict,
                  compartments, species, reactions, wrapper: str) -> TransformGraph:
        node_model = Node(model_id, "SBMLModel", model_metadata)

        l_node_species = [self._transform_species(model_id, x) for x in species['elements']]
        l_node_compartment = [self._transform_compartment(model_id, x) for x in compartments['elements']]
        l_node_reaction = [self._transform_reaction(model_id, x) for x in reactions['elements']]
        nodes = [node_model]
        nodes += l_node_species
        nodes += l_node_compartment
        nodes += l_node_reaction

        edges = self.transform_edges(node_model, nodes, reactions, wrapper)

        return nodes, edges

    @staticmethod
    def _transform(model_id: str, d: dict, label: str):
        node_id = f"{model_id}:{d['id']}"
        props = d.copy()
        return Node(node_id, label, props)

    @staticmethod
    def _transform_species(model_id: str, d: dict):
        return TransformCobraSBML._transform(model_id, d, 'SBMLSpecies')

    @staticmethod
    def _transform_compartment(model_id: str, d: dict):
        return TransformCobraSBML._transform(model_id, d, 'SBMLCompartment')

    @staticmethod
    def build_stoich_links(node_id, d_node_species: dict, xml_reaction: str, wrapper: str = None):
        s_tuples = {}

        d_reactants = parse_species_list(xml_reaction, list_tag='listOfReactants', wrapper=wrapper)
        d_products = parse_species_list(xml_reaction, list_tag='listOfProducts', wrapper=wrapper)

        for o in d_reactants:
            s_prop = dict(o)

            if 'species' not in o:
                raise ValueError("species key not found")

            s_prop.pop('species')
            species_id = o['species']
            node_species = d_node_species[species_id]

            k = (('SBMLReaction', node_id), 'has_reactant', ('SBMLSpecies', node_species.key))
            if k not in s_tuples:
                s_tuples[k] = s_prop
            else:
                raise ValueError('dup')
        for o in d_products:
            s_prop = dict(o)

            if 'species' not in o:
                raise ValueError("species key not found")

            s_prop.pop('species')
            species_id = o['species']
            node_species = d_node_species[species_id]

            k = (('SBMLReaction', node_id), 'has_product', ('SBMLSpecies', node_species.key))
            if k not in s_tuples:
                s_tuples[k] = s_prop
            else:
                raise ValueError('dup')

        return s_tuples

    @staticmethod
    def build_stoich_links_old(node_id: str, d_node_species, d):
        s_tuples = {}
        for (species_id, stoich_value) in d['_reactants']:
            s_prop = None
            if stoich_value is not None:
                s_prop = {'stoichiometry': stoich_value}
            node_species = d_node_species[species_id]
            k = (node_id, 'has_reactant', node_species.key)
            if k not in s_tuples:
                s_tuples[k] = s_prop
            else:
                raise ValueError('dup')
            # print('r', species_id, stoich_value)

        for (species_id, stoich_value) in d['_products']:
            s_prop = None
            if stoich_value is not None:
                s_prop = {'stoichiometry': stoich_value}
            node_species = d_node_species[species_id]
            k = (node_id, 'has_product', node_species.key)
            if k not in s_tuples:
                s_tuples[k] = s_prop
            else:
                raise ValueError('dup')
        return s_tuples

    @staticmethod
    def _transform_reaction(model_id: str, d: dict):
        node_id = f"{model_id}:{d['id']}"
        props = d.copy()
        props.pop('_raw_notes', None)
        props.pop('_raw_annotation', None)
        props.pop('_raw_stoichiometry', None)
        props.pop('_products', None)
        props.pop('_reactants', None)

        return Node(node_id, 'SBMLReaction', props)
