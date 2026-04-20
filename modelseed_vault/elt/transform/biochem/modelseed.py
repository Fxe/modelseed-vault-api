from modelseed_vault.core.transform_graph import Node, TransformGraph


def clean_none(cpd):
    return {k: v for k, v in cpd.items() if v is not None}


class TransformModelSEEDBiochem:

    @staticmethod
    def transform_compound(compounds) -> dict:
        nodes = {}
        for _cpd in compounds:
            cpd = _cpd.copy()
            if 'aliases' in cpd:
                del cpd['aliases']
            if 'notes' in cpd and type(cpd['notes']) == list:
                cpd['notes'] = '; '.join(cpd['notes'])
            node = Node(cpd['id'], "ModelSEEDCompound", data=clean_none(cpd))
            nodes[node.id] = node
        return nodes

    @staticmethod
    def transform_reactions(reactions):
        # compound_eid = {}
        edges = []
        nodes = []

        for _rxn in reactions:
            rxn_id = _rxn['id']
            rxn = _rxn.copy()
            del rxn['id']

            delete_keys = set()
            for k in rxn:
                if rxn[k] is None:
                    delete_keys.add(k)
            for k in delete_keys:
                del rxn[k]

            stoich = []
            if 'aliases' in rxn:
                del rxn['aliases']
            if 'notes' in rxn and type(rxn['notes']) == list:
                rxn['notes'] = '; '.join(rxn['notes'])
            if 'pathways' in rxn and type(rxn['pathways']) == list:
                rxn['pathways'] = '; '.join(rxn['pathways'])
            if 'ec_numbers' in rxn and type(rxn['ec_numbers']) == list:
                rxn['ec_numbers'] = '; '.join(rxn['ec_numbers'])
            if 'stoichiometry' in rxn:
                stoich = rxn['stoichiometry']
                del rxn['stoichiometry']
            if 'thermodynamics' in rxn:
                if 'Group contribution' in rxn['thermodynamics']:
                    rxn['deltag_group_contribution'] = rxn['thermodynamics']['Group contribution'][0]
                    rxn['deltagerr_group_contribution'] = rxn['thermodynamics']['Group contribution'][1]
                if 'eQuilibrator' in rxn['thermodynamics']:
                    rxn['deltag_eQuilibrator'] = rxn['thermodynamics']['eQuilibrator'][0]
                    rxn['deltagerr_eQuilibrator'] = rxn['thermodynamics']['eQuilibrator'][1]
                del rxn['thermodynamics']

            node = Node(rxn_id, "ModelSEEDReaction", data=clean_none(rxn))
            nodes.append(node)

            for o in stoich:
                edges.append([('ModelSEEDReaction', rxn_id), 'has_stoichiometry_coefficient',
                              ('ModelSEEDCompound', o['compound']), o])
        """
        for e in edges:
            if e[2][1] not in compound_eid:
                raise ValueError(f"missing eid {e[2][1]}")
                #compound_node = vault.get_node('ModelSEEDCompound', e[2][1])
                #if compound_node is not None:
                #    compound_eid[compound_node['entry']] = compound_node['elementId']
    
    
        for e in edges:
            reaction_node = vault.get_node('ModelSEEDReaction', e[0][1])
            if reaction_node is not None:
                vault.add_edge(reaction_node['elementId'], compound_eid[e[2][1]], e[1], e[3])
        """

        return nodes, edges

    def transform(self, compounds, reactions) -> TransformGraph:
        g = TransformGraph()

        compound_nodes = TransformModelSEEDBiochem.transform_compound(compounds)
        for node in compound_nodes.values():
            g.add_transform_node(node)

        nodes, edges = TransformModelSEEDBiochem.transform_reactions(reactions)
        for n in nodes:
            _found = g.get_node_by_key_label_pair(n.key, n.primary_label)
            if _found is None:
                g.add_transform_node(n)

        for e in edges:
            src = e[0]
            src_node = g.get_node_by_key_label_pair(src[1], src[0])
            if src_node is None:
                src_node = Node(src[1], src[0])
                g.add_transform_node(src_node)
            dst = e[2]
            dst_node = g.get_node_by_key_label_pair(dst[1], dst[0])
            if dst_node is None:
                dst_node = Node(dst[1], dst[0])
                g.add_transform_node(dst_node)
            g.add_transform_edge(src_node, dst_node, e[1], data=e[3])

        return g
