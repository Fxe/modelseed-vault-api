from networkx import DiGraph, compose
from modelseed_vault.utils import sha_hex


class Node:
    def __init__(self, key: str, label: str, data=None):
        if not key:
            raise ValueError('empty key')
        if not label:
            raise ValueError('empty label')
        self._key = key.strip()
        self.label = label.strip()
        self.data = data if data else {}

    @property
    def key(self):
        return self._key.replace(" ", "_")

    @property
    def id(self):
        return f"{self.label}/{self.key}"

    def __hash__(self):
        return hash(self.id)

    def __eq__(self, other):
        return isinstance(other, Node) and self.id == other.id

    def to_json(self):
        out = {
            '_key': self._key
        }
        out.update(self.data)

        return out


class NodeHash(Node):

    def __init__(self, key, label, data=None):
        self._value = key
        super().__init__(sha_hex(key), label, data)

    def to_json(self):
        out = super().to_json()
        out['_value'] = self._value
        return out


class TransformGraph(DiGraph):
    def __init__(self, incoming_graph_data=None):
        super().__init__(incoming_graph_data)
        self.t_nodes = {}
        self.t_edges = {}

    def concat(self, graph):
        res = compose(self, graph)
        res.t_nodes.update(self.t_nodes)
        for klass in graph.t_nodes:
            if klass not in res.t_nodes:
                res.t_nodes[klass] = {}
                res.t_nodes[klass].update(graph.t_nodes[klass])
            else:
                for k, v in graph.t_nodes[klass].items():
                    res.t_nodes[klass][k] = v
        res.t_edges.update(self.t_edges)
        for klass in graph.t_edges:
            if klass not in res.t_edges:
                res.t_edges[klass] = {}
                res.t_edges[klass].update(graph.t_edges[klass])
            else:
                for k, v in graph.t_edges[klass].items():
                    res.t_edges[klass][k] = v
        #res.t_nodes.update(graph.t_nodes)
        return res

    def add_transform_edge(self, src, dst, label, data=None):
        l1 = list(filter(lambda x: x.id == src, self.nodes))
        l2 = list(filter(lambda x: x.id == dst, self.nodes))
        if len(l1) == 1 and len(l2) == 1:
            if label not in self.t_edges:
                self.t_edges[label] = {}
            self.t_edges[label][(l1[0], l2[0])] = data
            self.add_edge(l1[0], l2[0], data=data if data else {})

    def add_transform_edge2(self, src, dst, label, data=None):
        if src in self.nodes and dst in self.nodes:
            if label not in self.t_edges:
                self.t_edges[label] = {}
            self.t_edges[label][(src, dst)] = data
            self.add_edge(src, dst, data=data if data else {})

    def add_transform_node(self, node_id, label, data=None):
        node = Node(node_id, label, data)

        return self.add_transform_node2(node)

    def add_transform_node2(self, node: Node):
        if node.label not in self.t_nodes:
            self.t_nodes[node.label] = {}
        if node.id not in self.t_nodes[node.label]:
            self.t_nodes[node.label][node.id] = node
            self.add_node(node)
        else:
            return self.t_nodes[node.label][node.id]
            #raise Exception('dup')

        return node

    def summary(self):
        for k in self.t_nodes:
            print('N', k, len(self.t_nodes[k]))
        for k in self.t_edges:
            print('E', k, len(self.t_edges[k]))
