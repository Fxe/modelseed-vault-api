import json
from pathlib import Path


class ExtractModelSEEDBiochem:

    def __init__(self, root: Path):
        self.root = root
        self.path_biochem = self.root / "Biochemistry"
        if not self.path_biochem.exists() or not self.path_biochem.is_dir():
            raise ValueError(f'invalid root {self.root}')

    def extract(self):
        compounds = []
        for f in self.path_biochem.iterdir():
            if f.name.startswith('compound_') and f.name.endswith('.json'):
                with open(str(f)) as fh:
                    compounds += json.load(fh)

        reactions = []
        for f in self.path_biochem.iterdir():
            if f.name.startswith('reaction_') and f.name.endswith('.json'):
                with open(str(f)) as fh:
                    reactions += json.load(fh)

        return compounds, reactions
