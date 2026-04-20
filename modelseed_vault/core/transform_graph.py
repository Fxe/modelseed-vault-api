from networkx import DiGraph, compose
from modelseed_vault.utils import sha_hex


class Node:
    """A labeled graph node identified by a (primary_label, key) pair.

    A node may carry multiple Neo4j labels, but exactly one — ``primary_label``
    — has the uniqueness constraint used for key-based lookups.

    Attributes:
        primary_label: The label that carries the unique key constraint
            (e.g. 'Genome'). Used to construct ``id`` and for Vault API calls.
        labels: Full list of Neo4j labels on this node (e.g.
            ``['Label1', 'Label2', 'Label3']``), where ``primary_label`` is the
            one with the constraint. ``None`` when the node has only one label.
        data: Arbitrary properties stored on the node.
        eid: Neo4j element ID — surrogate key assigned by Neo4j after the node
            is persisted. None until the node has been written to the database.
    """

    def __init__(self, key: str, primary_label: str, labels=None, data=None, eid: str = None):
        """
        Args:
            key: Domain identifier for the node. Spaces are normalised to
                underscores via the ``key`` property.
            primary_label: The label with the unique key constraint.
            labels: Full list of Neo4j labels on this node. Pass ``None`` when
                the node has only ``primary_label``.
            data: Optional dict of node properties.
            eid: Neo4j element ID. Leave as ``None`` when constructing a node
                before it has been persisted.
        """
        if not key:
            raise ValueError('empty key')
        if not primary_label:
            raise ValueError('empty primary label')
        self._key = key.strip()
        self.primary_label = primary_label.strip()
        self.labels = labels
        self.data = data if data else {}
        self.eid = eid  # surrogate key assigned by Neo4j

    @property
    def key(self) -> str:
        """Domain key with spaces replaced by underscores."""
        return self._key.replace(" ", "_")

    @property
    def id(self) -> str:
        """Composite identifier: ``<label>/<key>``."""
        return f"{self.primary_label}/{self.key}"

    def __hash__(self):
        return hash(self.id)

    def __eq__(self, other):
        return isinstance(other, Node) and self.id == other.id

    def to_json(self) -> dict:
        """Serialise the node to a dict suitable for the Vault API.

        Returns:
            Dict containing ``_key``, ``_primary_label``, ``_labels`` (the
            full label list, or ``None`` for single-label nodes), and all
            entries from ``data``.
        """
        out = {
            '_key': self._key,
            '_primary_label': self.primary_label,
            '_labels': self.labels,
        }
        out.update(self.data)
        return out


class Edge:
    """A directed edge between two nodes in a TransformGraph.

    Attributes:
        src: Source node.
        dst: Destination node.
        label: Edge type (e.g. 'has_genome', 'has_reactant').
        data: Arbitrary properties stored on the edge.
        eid: Neo4j element ID — surrogate key assigned by Neo4j after the edge
            is persisted. None until the edge has been written to the database.
    """

    def __init__(self, src: Node, dst: Node, label: str, data=None, eid: str = None):
        """
        Args:
            src: Source node.
            dst: Destination node.
            label: Edge type label.
            data: Optional dict of edge properties.
            eid: Neo4j element ID. Leave as ``None`` before the edge is
                persisted.
        """
        self.src = src
        self.dst = dst
        self.label = label
        self.data = data if data else {}
        self.eid = eid  # surrogate key assigned by Neo4j

    @property
    def id(self) -> str:
        """Composite identifier: ``<src.id>/<label>/<dst.id>``."""
        return f"{self.src.id}/{self.label}/{self.dst.id}"

    def __hash__(self):
        return hash(self.id)

    def __eq__(self, other):
        return isinstance(other, Edge) and self.id == other.id


class HashNode(Node):
    """A Node whose key is the SHA-256 hex digest of the original value.

    Neo4j enforces a size limit on key constraints. HashNode works around this
    by storing the SHA-256 hash as the ``key`` while preserving the original
    value in ``_value``.
    """

    def __init__(self, key: str, primary_label: str, labels=None, data=None):
        """
        Args:
            key: Original string value. The stored key will be its SHA-256
                hex digest.
            primary_label: The label with the unique key constraint.
            labels: Full list of Neo4j labels. Pass ``None`` for single-label
                nodes.
            data: Optional dict of node properties.
        """
        self._value = key
        super().__init__(sha_hex(key), primary_label, labels, data)

    def to_json(self) -> dict:
        """Serialise the node, including the original ``_value`` field."""
        out = super().to_json()
        out['_value'] = self._value
        return out


class TransformGraph(DiGraph):
    """An in-memory subgraph of Neo4j backed by a networkx DiGraph.

    Nodes and edges are indexed by label/type in ``t_nodes`` and ``t_edges``
    so that callers can iterate over all nodes or edges of a given type without
    scanning the full graph. Use ``concat`` to merge two TransformGraphs before
    writing to the Vault.

    Attributes:
        t_nodes: ``{primary_label: {node.id: Node}}`` — nodes indexed by
            primary label.
        t_edges: ``{edge_label: {edge.id: Edge}}`` — edges indexed by type.
    """

    def __init__(self, incoming_graph_data=None):
        super().__init__(incoming_graph_data)
        self.t_nodes = {}
        self.t_edges = {}

    def get_node_by_key_label_pair(self, key, label) -> Node:
        return self.t_nodes.get(label, {}).get(f"{label}/{key}")

    def concat(self, graph: "TransformGraph") -> "TransformGraph":
        """Merge another TransformGraph into this one and return the result.

        Uses networkx ``compose`` for the underlying DiGraph merge, then
        reconciles ``t_nodes`` and ``t_edges`` from both graphs. Duplicate
        keys from ``graph`` overwrite entries from ``self``.

        Args:
            graph: The TransformGraph to merge in.

        Returns:
            A new TransformGraph containing nodes and edges from both graphs.
        """
        res = compose(self, graph)
        res.t_nodes.update(self.t_nodes)
        for klass in graph.t_nodes:
            if klass not in res.t_nodes:
                res.t_nodes[klass] = {}
            for k, v in graph.t_nodes[klass].items():
                res.t_nodes[klass][k] = v
        res.t_edges.update(self.t_edges)
        for klass in graph.t_edges:
            if klass not in res.t_edges:
                res.t_edges[klass] = {}
            for k, v in graph.t_edges[klass].items():
                res.t_edges[klass][k] = v
        return res

    def add_transform_node(self, node: Node) -> Node:
        """Add a node to the graph if it does not already exist.

        If a node with the same ``id`` is already present, the existing node
        is returned without modification.

        Args:
            node: The Node to add.

        Returns:
            The node stored in the graph (either ``node`` or the pre-existing
            one if a duplicate was detected).
        """
        if node.primary_label not in self.t_nodes:
            self.t_nodes[node.primary_label] = {}
        if node.id not in self.t_nodes[node.primary_label]:
            self.t_nodes[node.primary_label][node.id] = node
            self.add_node(node)
        else:
            return self.t_nodes[node.primary_label][node.id]
        return node

    def add_transform_edge(self, src: Node, dst: Node, label: str, data=None) -> Edge:
        """Add a directed edge between two existing nodes.

        The edge is only added when both ``src`` and ``dst`` are already
        present in the graph. Callers should add nodes via
        ``add_transform_node`` before calling this method.

        Args:
            src: Source node.
            dst: Destination node.
            label: Edge type label (e.g. 'has_genome', 'has_reactant').
            data: Optional dict of edge properties.

        Returns:
            The Edge object stored in the graph, or ``None`` if either node
            was not found.
        """
        if src in self.nodes and dst in self.nodes:
            edge = Edge(src, dst, label, data)
            if label not in self.t_edges:
                self.t_edges[label] = {}
            self.t_edges[label][edge.id] = edge
            self.add_edge(src, dst, data=edge.data)
            return edge
        return None

    def summary(self):
        """Print a count of nodes and edges grouped by label/type."""
        for k in self.t_nodes:
            print('N', k, len(self.t_nodes[k]))
        for k in self.t_edges:
            print('E', k, len(self.t_edges[k]))
