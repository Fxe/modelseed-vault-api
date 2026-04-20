from modelseed_vault.vault import Vault
from modelseed_vault.core.transform_graph import TransformGraph
from modelseed_vault.core.graph_ontology import GENOMIC_CONTIG, PROTEIN_SEQUENCE, DNA_SEQUENCE


class LoadNeo4j:

    def __init__(self, vault: Vault):
        self.vault = vault

    def load2(self, graph: TransformGraph):

        # find all nodes already in vault (these are proxies)
        refs = set()
        for node_type in graph.t_nodes:
            for node in graph.t_nodes[node_type].values():
                refs.add((node.primary_label, node.key))
        query_result = self.vault.query_eid(refs)

        node_to_eids = {(ref['type'], ref['key']): ref['elementId'] for ref in query_result if
                        ref['elementId'] is not None}

        # build load node payload for new nodes
        payload_nodes = []
        for node_type in graph.t_nodes:
            for node in graph.t_nodes[node_type].values():
                tkp = (node.primary_label, node.key)
                if tkp not in node_to_eids:
                    payload_nodes.append(node)

        bulk_add_node_result = self.vault.bulk_add_nodes2(payload_nodes)
        for k in bulk_add_node_result:
            primary_label, key = k.split('/')
            node_to_eids[(primary_label, key)] = bulk_add_node_result[k]

        payload_edges = []
        for edge_type in graph.t_edges:
            for o, u in graph.t_edges[edge_type].items():
                payload_edges.append([node_to_eids[(u.src.primary_label, u.src.key)],
                                      u.label,
                                      node_to_eids[(u.dst.primary_label, u.dst.key)],
                                      u.data])
                
        self.vault.bulk_add_edges2(payload_edges)

    def load(self, graph: TransformGraph):
        """Write a TransformGraph to Neo4j via the Vault REST API.

        Executes in two phases so that edge creation can reference the
        Neo4j-assigned element IDs (eid) produced during node creation.

        Phase 1 — nodes: iterates graph.t_nodes, POSTs each node via
        vault.add_node2, and stores the returned eid on node.eid.

        Phase 2 — edges: iterates graph.t_edges, POSTs each Edge via
        vault.add_edge using the node eids set in Phase 1, and stores the
        returned eid on edge.eid.

        Args:
            graph: The TransformGraph to persist.

        Raises:
            requests.HTTPError: If any Vault API call returns a non-2xx status.
            ValueError: If a node or edge response does not contain an eid.
        """
        # Phase 1 — load nodes
        node_eid = {}
        for label in graph.t_nodes:
            for node in graph.t_nodes[label].values():
                # GENOMIC_CONTIG DNA_SEQUENCE
                # if node.primary_label in {GENOMIC_CONTIG, PROTEIN_SEQUENCE, DNA_SEQUENCE}:
                if node.primary_label in {PROTEIN_SEQUENCE}:
                    eid = self.vault.add_protein(node)
                else:
                    res = self.vault.add_node2(node)
                    res.raise_for_status()
                    data = res.json()
                    if 'elementId' not in data:
                        raise ValueError(f"No eid in response for node {node.id}")
                    eid = data['elementId']
                if node.id not in node_eid:
                    node_eid[node.id] = eid
                else:
                    raise ValueError('dup should not happen')

        print('total_eid', len(node_eid))

        # Phase 2 — load edges
        for edge_label in graph.t_edges:
            for edge in graph.t_edges[edge_label].values():
                res = self.vault.add_edge(node_eid[edge.src.id],
                                          node_eid[edge.dst.id],
                                          edge_label, edge.data)
                res.raise_for_status()
                eid = res.content.decode('utf-8')
                edge.eid = eid
