import csv
from collections import defaultdict
from modelseed_vault.core.hash import _hash_string
from modelseedpy import RastClient
import json


def read_genes_tab(path):
    """Read an LBNL FitnessBrowser genes.tab file into a list of dicts."""
    genes = []
    with open(path, "r") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            genes.append({
                "locusId": row["locusId"],
                "sysName": row.get("sysName", ""),
                "scaffoldId": row["scaffoldId"],
                "begin": int(row["begin"]),
                "end": int(row["end"]),
                "strand": row["strand"],
                "type": row.get("type", "1"),
                "desc": row.get("desc", ""),
            })
    return genes


def read_proteins(path):
    """
    Read a proteins.faa FASTA file.
    Returns dict mapping locusId -> sequence string.
    Header format: >genome:locusId  description
    """
    proteins = {}
    current_id = None
    seq_parts = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    proteins[current_id] = "".join(seq_parts)
                header = line[1:].split()[0]  # e.g. "acidovorax_3H11:Ac3H11_1"
                current_id = header.split(":", 1)[-1]  # e.g. "Ac3H11_1"
                seq_parts = []
            else:
                seq_parts.append(line)
        if current_id is not None:
            proteins[current_id] = "".join(seq_parts)
    return proteins


def generate_adjacent_edges_rast_ec_annotation(proteins: dict, rast: RastClient):
    result = rast.annotate_protein_sequences(proteins)
    """
    example result:
     
    {'protein1': 'Transcriptional regulator, AraC family',
      'protein2': 'Gamma-glutamyltranspeptidase (EC 2.3.2.2) @ Glutathione hydrolase (EC 3.4.19.13)'}
    """
    import json
    with open('S:/python3/modelseed-vault-api/modelseed_vault/data/function_to_ec.json', 'r') as fh:
        ec_data = json.load(fh)

    for protein_id, annotation in result.items():
        annotation_ec = ec_data.get(_hash_string(annotation))



def generate_adjacent_edges(genes):
    """
    Given a list of gene dicts (from read_genes_tab), return a list of
    ADJACENT_TO edge dicts matching the vault chain JSON schema.

    Edges are generated for each consecutive gene pair on the same scaffold,
    sorted by begin position.
    """
    by_scaffold = defaultdict(list)
    for g in genes:
        by_scaffold[g["scaffoldId"]].append(g)

    edges = []
    for scaffold_id, scaffold_genes in sorted(by_scaffold.items()):
        scaffold_genes.sort(key=lambda g: g["begin"])
        contig = scaffold_id

        for i in range(len(scaffold_genes) - 1):
            src = scaffold_genes[i]
            tgt = scaffold_genes[i + 1]
            gap = max(0, tgt["begin"] - src["end"] - 1)
            edges.append({
                "source": {"entry": src["locusId"], "elementId": None},
                "contig": contig,
                "type": "ADJACENT_TO",
                "target": {"entry": tgt["locusId"], "elementId": None},
                "properties": {
                    "source_start": src["begin"],
                    "source_end": src["end"],
                    "source_strand": src["strand"],
                    "gap": gap,
                    "target_start": tgt["begin"],
                    "target_end": tgt["end"],
                    "target_strand": tgt["strand"],
                },
            })

    return edges
