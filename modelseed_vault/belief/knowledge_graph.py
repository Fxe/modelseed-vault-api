from modelseed_vault.core.transform_graph import Node
import networkx as nx
from modelseed_vault.vault import Vault


class KnowledgeGraphBuilder:
    def __init__(self, vault: Vault):
        self.vault = vault
        self.graph = nx.DiGraph()
        self.node_data = {}
        self.edge_data = {}
        self.visited = set()

    def add_node(self, node: Node):
        if node.id not in self.graph.nodes:
            self.graph.add_node(node.id)
        self.node_data[node.id] = node.data

    def add_edge(self, src, dst, vocab, eid, data):
        self.edge_data[eid] = data
        if src in self.graph.nodes and dst in self.graph.nodes:
            self.graph.add_edge(src, dst, ontology=vocab, element_id=eid)
        else:
            print(f"!! not found {src} {dst}")

    @staticmethod
    def _make_node(raw):
        return Node(raw['entry'], raw['labels'][0], data=raw)

    def expand(self, node: Node, child=True, parent=True):
        self.visited.add(node.id)
        ret = {}

        # children
        if child:
            for edge, raw in self.vault.get_node_child(node):
                child = KnowledgeGraphBuilder._make_node(raw)
                self.add_node(child)
                self.add_edge(node.id, child.id, edge['t'], edge['elementId'], edge['properties'])
                # print(f"{node.id} -[{edge['t']}]-> {child.id}")

                ret.setdefault(child.primary_label, {})
                if child.id not in self.visited:
                    ret[child.primary_label][child.key] = child

        # parents
        if parent:
            for edge, raw in self.vault.get_node_parent(node):
                parent = KnowledgeGraphBuilder._make_node(raw)
                self.add_node(parent)
                self.add_edge(parent.id, node.id, edge['t'], edge['elementId'], edge['properties'])
                # print(f"{node.id} <-[{edge['t']}]- {parent.id}")

                ret.setdefault(parent.primary_label, {})
                if parent.id not in self.visited:
                    ret[parent.primary_label][parent.key] = parent

        return ret

    def expand_recursive(self, node, allowed, max_depth=3, _depth=0):
        """
        Recursively expand a node through the knowledge graph.

        Args:
            node:       starting Node object
            allowed:    set of primary_labels eligible for further expansion
            max_depth:  stop recursing after this many levels
            _depth:     (internal) current recursion depth

        Returns:
            dict  –  {depth: {primary_label: {key: Node, ...}, ...}, ...}
        """
        if _depth >= max_depth:
            return {}

        exp_c, exp_p = allowed[node.primary_label]
        exp_nodes = self.expand(node, exp_c, exp_p)

        # for k, n in exp_nodes.items():
        #    print(k, nodes)

        next_nodes = [
            n
            for label, nodes in exp_nodes.items()
            if label in allowed
            for n in nodes.values()
        ]

        result = {_depth: exp_nodes}

        if not next_nodes:
            return result

        for child in next_nodes:
            child_result = self.expand_recursive(
                child, allowed, max_depth, _depth + 1
            )
            result.update(child_result)

        return result

    def build(self, label, entry, allowed, max_depth=3):
        """
        Convenience method: fetch a node and fully expand it.

        Args:
            label:      node label (e.g. 'GenomicFeature')
            entry:      node entry identifier
            allowed:    set of primary_labels to recurse into
            max_depth:  maximum expansion depth

        Returns:
            dict of results keyed by depth
        """
        raw = self.vault.get_node(label, entry)
        root = self._make_node(raw)
        self.add_node(root)
        return self.expand_recursive(root, allowed, max_depth)

    def reset(self):
        """Clear the graph and visited state for a fresh run."""
        self.graph.clear()
        self.node_data.clear()
        self.visited.clear()
