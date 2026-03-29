from neo4j import GraphDatabase
import uuid
import json


class Neo4jDAO:

    def __init__(self, uri: str, user: str, password: str):
        self.driver = GraphDatabase.driver(uri, auth=(user, password))

    def close(self):
        self.driver.close()

    @staticmethod
    def _get_node(tx, node_id, node_label):
        # print(type(tx), node_id, node_label)
        result = tx.run(f"""
            MATCH (n:{node_label})
            WHERE elementId(n) = $node_id
            RETURN n
            """, node_id=node_id
        )
        return result.single()

    def create_node(self, labels, node_key):
        _param_label = ':'.join(labels)
        cypher_merge_node = f"""
        MERGE (n:{_param_label} {{key: '{node_key}'}})
        ON CREATE
          SET n.created_at = timestamp(), n.updated_at = timestamp()
        ON MATCH
          SET n.updated_at = timestamp()
        RETURN n
        """
        return self.driver.execute_query(cypher_merge_node)

    def _create_edge(self, rel_type, node_from_label, node_from_key, node_to_label, node_to_key, ext_args):
        cypher_merge_r = f"""
        MATCH (n_from:{node_from_label} {{key: '{node_from_key}'}})
        MATCH (n_to:{node_to_label} {{key: '{node_to_key}'}})
        MERGE (n_from)-[r:{rel_type}]->(n_to)
        ON CREATE
          SET r.created_at = timestamp(), r.updated_at = timestamp() {ext_args}
        ON MATCH
          SET r.updated_at = timestamp()    
        RETURN r
        """
        res = self.driver.execute_query(cypher_merge_r)
        return res

    def create_edge(self, rel_type, node_from, node_to, ext_args):
        query = f"""
        MATCH (n_from:{node_from.key_label} {{key: '{node_from.key}'}})
        MATCH (n_to:{node_to.key_label} {{key: '{node_to.key}'}})
        MERGE (n_from)-[r:{rel_type}]->(n_to)
        ON CREATE
          SET r.created_at = timestamp(), r.updated_at = timestamp() {ext_args}
        ON MATCH
          SET r.updated_at = timestamp()    
        RETURN r
        """

        res = self.driver.execute_query(query)
        return res

    @staticmethod
    def _get_node_by_key(tx, node_id, node_label):
        # print(type(tx), node_id, node_label)
        result = tx.run(f"""
            MATCH (n:{node_label})
            WHERE n.key = $node_id
            RETURN n
            """, node_id=node_id
        )
        return result.single()

    @staticmethod
    def _query(tx, query):
        result = tx.run(query)
        return [o for o in result]

    def get_node(self, node_id: str, node_label: str) -> dict:
        with self.driver.session(database="neo4j") as session:
            return session.execute_read(Neo4jDAO._get_node, node_id, node_label)

    def get_node_by_key(self, node_key: str, node_label: str) -> dict:
        with self.driver.session(database="neo4j") as session:
            return session.execute_read(Neo4jDAO._get_node_by_key, node_key, node_label)

    def get_or_create_node_by_key(self, node_key: str, node_label: str) -> dict:
        res = self.get_node_by_key(node_key, node_label)
        if res is None:
            return self._create_node([])
        return res

    def query(self, query: str):
        with self.driver.session(database="neo4j") as session:
            return session.execute_read(Neo4jDAO._query, query)


class GraphNode:

    def __init__(self, node_id, node_key, key_label, labels, dao):
        self.id = node_id
        self.key = node_key
        self.key_label = key_label
        self.labels = labels
        self.dao = dao
        self.created_at = None
        self.updated_at = None

    @staticmethod
    def from_neo4j_node(node):
        return GraphNode(node.element_id, node.get('key'), node.labels)

    def __str__(self):
        _str_label = ':'.join(self.labels)
        return f'<GraphNode:{_str_label} - {self.key} [{self.id}]'

    def __repr__(self):
        _str_label = ':'.join(self.labels)
        return f'<GraphNode:{_str_label} - {self.key} [{self.id}]'

    @staticmethod
    def _to_prop(k, v):
        if k == 'created_at':
            return 'created_at: timestamp()'
        elif type(v) == str:
            return f"{k}: '{v}'"
        else:
            return f"{k}: {v}"

    def delete(self):
        if self.id is not None:
            query = f"MATCH (n) WHERE elementId(n) = '{self.id}' DELETE n;"
            with self.dao.driver.session() as session:
                res = session.run(query)
                return res

    def fetch(self):
        if self.id is None:
            query = f"MATCH (n:{self.key_label}) WHERE n.key = '{self.key}' RETURN n;"
            with self.dao.driver.session() as session:
                single = session.run(query).single()
                if single:
                    self.id = single['n'].element_id
                    return self
                return None

    def create(self):
        if self.id is None:
            response = self.dao.create_node(self.labels, self.key)
            node = response[0][0]['n']
            self.id = node.element_id
            return self


class GraphNodeRastFunction(GraphNode):

    def __init__(self, node_id, rast_function, dao: Neo4jDAO):
        super().__init__(node_id, rast_function, 'ProteinAnnotation', ['ProteinAnnotation', 'RAST'], dao)


class GraphNodeProtein(GraphNode):

    def __init__(self, node_id, rast_function, dao: Neo4jDAO):
        super().__init__(node_id, rast_function, 'ProteinSequence', ['ProteinSequence'], dao)


class GraphNodeRastExecution(GraphNode):

    def __init__(self, node_id, node_key, host_name, analysis_execs, dao):
        super().__init__(node_id, node_key, 'Method', ['Method'], dao)
        self.host_name = host_name
        self.analysis_execs = analysis_execs

    @staticmethod
    def from_rast(analysis_events):
        host_name = None
        analysis_execs = []
        for event in analysis_events:
            if host_name is None:
                host_name = event['hostname']
            elif host_name != event['hostname']:
                raise ValueError('! events with different host')

            analysis_exec = {}
            for k in event:
                if k != 'hostname':
                    v = event[k]
                    if k == 'parameters':
                        v = ' '.join([str(x) for x in event[k]])
                    analysis_exec[k] = v
            analysis_execs.append(analysis_exec)
        return GraphNodeRastExecution(None, f'RAST_API_{str(uuid.uuid4())}', host_name, analysis_execs)

    def create(self):
        label_str = ':'.join(self.labels)
        props = {
            'key': self.key,
            'host_name': self.host_name,
            'analysis_execs': [json.dumps(x) for x in self.analysis_execs],
            'created_at': None
        }
        prop_str = ', '.join([self._to_prop(k, v) for k, v in props.items()])
        query = f"""
        CREATE (n:{label_str} {{{prop_str}}}) RETURN n
        """
        with self.dao.driver.session() as session:
            return session.run(query)
