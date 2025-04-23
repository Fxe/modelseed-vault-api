from modelseed_annotation_api.dao_neo4j import Neo4jDAO
from modelseed_annotation_api.dao_neo4j import GraphNodeRastFunction, GraphNodeProtein

"""
Main client class for interacting with ModelSEED annotation services.
"""


class ModelSEEDAnnotationClient:
    """
    Client for interacting with ModelSEED annotation services.
    
    This class provides methods to interact with various ModelSEED annotation
    services and endpoints.
    """
    
    def __init__(self, dao_vault):
        """
        Initialize the ModelSEED annotation client.
        """
        self.dao_vault = dao_vault

    def get_genome_set(self, key=None, node_id=None):
        query_cypher = f'MATCH (n:GenomeSet)-[r:has_genome]-(g:RefSeqGenome)  WHERE n.key = "{key}" RETURN g;'
        if node_id:
            query_cypher = f'MATCH (n:GenomeSet)-[r:has_genome]-(g:RefSeqGenome)  ' \
                           f'WHERE elementId(n) = "{node_id}" RETURN g;'

        node_elements = self.dao_vault.query(query_cypher)
        node_set = self.dao_vault.get_node(node_id, "GenomeSet")

        return node_set, node_elements

    def get_genome(self):
        pass

    def get_rast_annotation(self, k_v, method_id):
        for h, rast_o in k_v.items():
            rast_function = rast_o.get('function')
            if rast_function is None:
                rast_function = 'NONE'
                ext = ", r.method = '{method_id}'"
            else:
                import json
                ext = f", r.quality = '{json.dumps(rast_o['quality'])}', r.method = '{method_id}', r.rast_method = '{rast_o['annotations'][0][1]}'"

            rast_function = rast_function.strip()
            node_protein = GraphNodeProtein(None, h, self.dao_vault).fetch()
            if node_protein is None:
                node_protein = GraphNodeProtein(None, h, self.dao_vault).create()

            node_rast = GraphNodeRastFunction(None, rast_function, self.dao_vault).fetch()
            if node_rast is None:
                node_rast = GraphNodeRastFunction(None, rast_function, self.dao_vault).create()

            self.dao_vault.create_edge('has_annotation', node_protein, node_rast, ext)

    def get_protein_ontology(self, protein_id: str) -> dict:
        """
        Get the annotation for a given protein ID.
        """
        raise NotImplementedError("This method will be implemented in future versions")
    
    def get_protein_cluster(self, cluster_id: str, cluster_type: str) -> dict:
        """
        Get the annotation for a given protein cluster ID.
        """
        
        raise NotImplementedError("This method will be implemented in future versions")

    def annotate_sequence(self, sequence: str) -> dict:
        """
        Annotate a given sequence using ModelSEED services.
        
        Args:
            sequence (str): The sequence to annotate
            
        Returns:
            dict: Annotation results
        """
        
        raise NotImplementedError("This method will be implemented in future versions")
    
    def get_annotation(self, annotation_id: str) -> dict:
        """
        Retrieve a specific annotation by ID.
        
        Args:
            annotation_id (str): The ID of the annotation to retrieve
            
        Returns:
            dict: The annotation data
        """
        raise NotImplementedError("This method will be implemented in future versions") 