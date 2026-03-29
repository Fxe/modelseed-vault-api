from modelseed_vault.elt.transform.cobra.parse import parse_elements_with_provenance, parse_reactions_with_provenance
from pathlib import Path


class ExtractCobraJson:

    def __init__(self):
        pass

    def extract(self, json_file: Path):
        return None


class ExtractCobraSBML:

    def __init__(self):
        pass

    def extract(self, sbml_file: Path):
        with open(str(sbml_file), 'rb') as fh:
            species = parse_elements_with_provenance(fh, 'species')
        with open(str(sbml_file), 'rb') as fh:
            compartment = parse_elements_with_provenance(fh, 'compartment')
        xpath = "//*[local-name()='model']/*[local-name()='listOfReactions']/*[local-name()='reaction']"
        with open(str(sbml_file), 'rb') as fh:
            reaction = parse_reactions_with_provenance(fh, xpath=xpath)

        return {'compartments': compartment, 'species': species, 'reactions': reaction}
