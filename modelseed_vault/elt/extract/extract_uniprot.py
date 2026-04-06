from pathlib import Path
from modelseed_vault.elt.extract.parser_uniprot import UniprotProtParser
import gzip

class ExtractContainerUniprot:

    def __init__(self, kind:str):
        self.kind = None # swissprot or uniprot

class ExtractUniprot:

    def __init__(self):
        pass

    def extract(self, xml_path: Path) -> ExtractContainerUniprot:
        parser = UniprotProtParser()
        parser.set_xml_ns('{https://uniprot.org/uniprot}')
        # parser.max_entries = 1
        """
        if xml_path # endswith gz
            gzip.open

        for o in parser.parse(fh):
            entry_list.append(o)

        for o in entry_list:
            o
        
        """