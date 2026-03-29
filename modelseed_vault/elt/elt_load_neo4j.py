from modelseed_vault.vault import Vault
from modelseed_vault.core.transform_graph import TransformGraph


class LoadNeo4j:

    def __init__(self, vault: Vault):
        self.vault = vault


    def load(self, graph: TransformGraph):
        pass