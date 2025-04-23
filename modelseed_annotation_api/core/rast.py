from modelseedpy.core.msgenome import normalize_role
from modelseedpy import RastClient


class RastAnnotation:

    def __init__(self, annotation):
        self.annotation = annotation
        self.functions = {}


class RastFunction:

    def __init__(self, rast_function):
        self.search_string = normalize_role(rast_function)
        self.rast_function = rast_function


class VaultRastClient(RastClient):

    def __init__(self):
        super().__init__()

    def annotate_genome(self, genome, split_terms=True):
        pass
