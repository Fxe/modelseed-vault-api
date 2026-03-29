from modelseed_vault.dao_neo4j import Neo4jDAO


class AnnotationFunction:
    
    def __init__(self, function_id, value):
        """
        Initialize an AnnotationFunction object.

        Args:
            function_id (str): The ID of the function
            value (str): The value of the function
        """
        self.id = function_id
        self._value = value
        self.search_value = value
        self.synonyms = set()
        self.sub_functions = set()
        self.function_group = set()
        self.source = set()

    @staticmethod
    def from_json(data):
        """
        Create an AnnotationFunction object from a JSON dictionary.

        Args:
            data (dict): A dictionary containing the function data

        Returns:
            AnnotationFunction: An AnnotationFunction object
        """
        annotation_function = AnnotationFunction(data['id'], data['value'])
        annotation_function.search_value = data['search_value']
        annotation_function.synonyms |= set(data['synonyms'])
        annotation_function.function_group |= set(data['function_group'])
        annotation_function.source |= set(data['source'])
        for o in data['sub_functions']:
            sub_function = AnnotationFunction.from_json(o)
            annotation_function.sub_functions.add(sub_function)
        return annotation_function

    @property
    def value(self):
        return self._value

    def get_data(self):
        return {
            'id': self.id,
            'value': self._value,
            'search_value': self.search_value,
            'synonyms': list(sorted(self.synonyms)),
            'function_group': list(sorted(self.function_group)),
            'sub_functions': list(map(lambda x: x.get_data(), self.sub_functions)),
            'source': list(sorted(self.source))
        }
