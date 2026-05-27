from modelseed_vault.core.transform_graph import TransformGraph, Node
from modelseed_vault.core.genome import classify_protein_sequence, ProteinSequence, SeqType
from pathlib import Path
import polars as pl
from tqdm import tqdm
import json
import pyarrow as pa
import pyarrow.parquet as pq
from dataclasses import dataclass, field, asdict
from typing import Optional
import re


@dataclass
class DomainStats:
    """Statistics for a single region (full sequence or best domain)."""
    evalue: float
    score: float
    bias: float


@dataclass
class DomainEstimation:
    """Domain number estimation counts from HMMER tblout."""
    exp: float
    reg: int
    clu: int
    ov: int
    env: int
    dom: int
    rep: int
    inc: int


@dataclass
class KoFamHit:
    """A single hmmsearch hit against a KO profile."""
    target_name: str
    target_accession: str
    query_name: str
    query_accession: str
    full_sequence: DomainStats
    best_domain: DomainStats
    domain_estimation: DomainEstimation
    description: str = ""

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class KoFamBlock:
    """All hits and metadata for one KO profile search."""
    ko: str
    hits: list[KoFamHit] = field(default_factory=list)
    metadata: dict = field(default_factory=dict)

    @property
    def n_hits(self) -> int:
        return len(self.hits)


class KoFamParser:
    """Parser for KofamScan per-KO tabular (hmmsearch --tblout) output files.

    The input file contains one or more concatenated hmmsearch blocks.
    Each block is preceded by a bare KO identifier line (e.g. "K00001")
    and terminated by a "# [ok]" marker.
    """

    # Lines we skip inside a block: header/comment lines that aren't metadata
    _HEADER_COMMENT = re.compile(r"^#\s*(target name|-+|$)")
    # KO identifier lines are bare tokens like "K00001" (no '#', no whitespace split needed)
    _KO_LINE = re.compile(r"^(K\d{5})\s*$")
    # Metadata lines look like "# Program:  hmmsearch"
    _METADATA_LINE = re.compile(r"^#\s*([A-Z][A-Za-z ]+?):\s+(.*?)\s*$")

    def parse(self, file: Path) -> dict[str, KoFamBlock]:
        """Parse a KofamScan tabular output file.

        Returns a dict mapping KO identifier -> KoFamBlock (hits + metadata).
        If the same KO appears multiple times, later blocks overwrite earlier ones.
        """
        results: dict[str, KoFamBlock] = {}
        current: Optional[KoFamBlock] = None

        with open(file) as f:
            for raw_line in f:
                line = raw_line.rstrip("\n")
                if not line.strip():
                    continue

                # New KO block
                ko_match = self._KO_LINE.match(line)
                if ko_match:
                    current = KoFamBlock(ko=ko_match.group(1))
                    results[current.ko] = current
                    continue

                if current is None:
                    # Stray content before any KO header — skip
                    continue

                # End-of-block marker
                if line.strip() == "# [ok]":
                    current = None
                    continue

                # Metadata lines inside the trailing comment section
                meta_match = self._METADATA_LINE.match(line)
                if meta_match:
                    key, value = meta_match.group(1).strip(), meta_match.group(2)
                    current.metadata[key] = value
                    continue

                # Skip other comment/header lines
                if line.startswith("#"):
                    continue

                # Data row
                hit = self._parse_hit(line, current.ko)
                if hit is not None:
                    current.hits.append(hit)

        return results

    @staticmethod
    def _parse_hit(line: str, ko: str) -> Optional[KoFamHit]:
        """Parse a single data row from hmmsearch --tblout.

        Column layout (whitespace-separated):
          0  target_name
          1  target_accession
          2  query_name
          3  query_accession
          4  full_seq_evalue    5  full_seq_score    6  full_seq_bias
          7  best_dom_evalue    8  best_dom_score    9  best_dom_bias
         10  exp  11 reg  12 clu  13 ov  14 env  15 dom  16 rep  17 inc
         18+ description of target (may contain spaces)
        """
        # Split into max 19 fields so the description keeps internal whitespace
        fields = line.split(maxsplit=18)
        if len(fields) < 18:
            return None  # malformed row

        try:
            return KoFamHit(
                target_name=fields[0],
                target_accession=fields[1],
                query_name=fields[2],
                query_accession=fields[3],
                full_sequence=DomainStats(
                    evalue=float(fields[4]),
                    score=float(fields[5]),
                    bias=float(fields[6]),
                ),
                best_domain=DomainStats(
                    evalue=float(fields[7]),
                    score=float(fields[8]),
                    bias=float(fields[9]),
                ),
                domain_estimation=DomainEstimation(
                    exp=float(fields[10]),
                    reg=int(fields[11]),
                    clu=int(fields[12]),
                    ov=int(fields[13]),
                    env=int(fields[14]),
                    dom=int(fields[15]),
                    rep=int(fields[16]),
                    inc=int(fields[17]),
                ),
                description=fields[18] if len(fields) > 18 else "",
            )
        except (ValueError, IndexError):
            return None


def fix_label(lkp):
    label, key = lkp.split('/', 1)
    if label == 'GOTerm':
        label = 'OntologyGO'
    return f"{label}/{key}"


class KBERDL:

    def __init__(self, root_reference: Path, block: int):
        self.ldf_ccol = pl.scan_parquet(f'{root_reference}/{block}/contig_collection.parquet')
        self.ldf_contig = pl.scan_parquet(f'{root_reference}/{block}/contig.parquet')
        self.ldf_ccol_x_contig = pl.scan_parquet(f'{root_reference}/{block}/contig_x_contig_collection.parquet')
        self.ldf_contig_x_feature = pl.scan_parquet(f'{root_reference}/{block}/contig_x_feature.parquet')
        self.ldf_feature_x_protein = pl.scan_parquet(f'{root_reference}/{block}/feature_x_protein.parquet')
        self.ldf_protein = pl.scan_parquet(f'{root_reference}/{block}/protein.parquet')
        self.ldf_feature = pl.scan_parquet(f'{root_reference}/{block}/feature.parquet')
        self.ldf_name = pl.read_parquet(f'{root_reference}/{block}/name.parquet')
        self.load = None
        self.vault = None

    def register_labels(self):
        for label in {'Genome',
                      'IdentifierGenome', 'IdentifierContig',
                      'GenomicContig', 'GenomicFeature', 'LocusTag', 'ProteinSequence'}:
            res = self.vault.register(label)
            print(label, res)

    def etl_genomes(self):
        batch_size = 20000
        offset = 0

        total = self.ldf_ccol.select(pl.len()).collect().item()

        with tqdm(total=total, desc="Processing batches", unit="rows") as pbar:
            while True:
                batch = (
                    self.ldf_ccol
                    .slice(offset, batch_size)
                    .join(
                        self.ldf_name.lazy(),
                        left_on="contig_collection_id",
                        right_on="entity_id",
                        how="inner"
                    )
                    .collect()
                )

                g = TransformGraph()
                for row in batch.iter_rows(named=True):
                    node_genome_data = {
                        'ncbi_taxon_id': row['ncbi_taxon_id'],
                        'gtdb_taxon_id': row['gtdb_taxon_id'],
                        'contig_bp': row['contig_bp'],
                    }
                    node_genome = Node(row['hash'], "Genome", data=node_genome_data)
                    node_identifier1 = Node(row['name'], "IdentifierGenome")
                    node_identifier2 = Node(row['contig_collection_id'], "IdentifierGenome", labels=['KBERDL'])
                    g.add_transform_node(node_genome)
                    g.add_transform_node(node_identifier1)
                    g.add_transform_node(node_identifier2)
                    g.add_transform_edge(node_genome, node_identifier1, "has_identifier", data={'source': 'KBERDL'})
                    g.add_transform_edge(node_genome, node_identifier2, "has_identifier", data={'source': 'KBERDL'})
                if len(g.t_nodes) > 0:
                    self.load.load2(g)

                if batch.is_empty():
                    break

                offset += batch_size
                pbar.update(len(batch))

    def etl_contigs(self):
        batch_size = 20000
        offset = 0

        total = self.ldf_contig.select(pl.len()).collect().item()

        with tqdm(total=total, desc="Processing batches", unit="rows") as pbar:
            while True:
                batch = (
                    self.ldf_contig
                    .slice(offset, batch_size)
                    .join(
                        self.ldf_name.lazy(),
                        left_on="contig_id",
                        right_on="entity_id",
                        how="inner"
                    )
                    .collect()
                )

                g = TransformGraph()
                for row in batch.iter_rows(named=True):
                    node_genome_data = {
                        'hash': row['gc_content'],
                        'gc_content': row['gc_content'],
                        'length': row['length'],
                    }
                    node_contig = Node(row['contig_id'], "GenomicContig", data=node_genome_data)
                    node_identifier1 = Node(row['name'], "IdentifierContig")
                    node_identifier2 = Node(row['contig_id'], "IdentifierContig", labels=['KBERDL'])
                    g.add_transform_node(node_contig)
                    g.add_transform_node(node_identifier1)
                    g.add_transform_node(node_identifier2)
                    g.add_transform_edge(node_contig, node_identifier1, "has_identifier", data={'source': 'KBERDL'})
                    g.add_transform_edge(node_contig, node_identifier2, "has_identifier", data={'source': 'KBERDL'})
                if len(g.t_nodes) > 0:
                    self.load.load2(g)

                if batch.is_empty():
                    break

                offset += batch_size
                pbar.update(len(batch))

    def etl_genome_contig_edges(self):
        batch_size = 100000
        offset = 0

        total = self.ldf_ccol_x_contig.select(pl.len()).collect().item()

        with tqdm(total=total, desc="Processing batches", unit="rows") as pbar:
            while True:
                batch = (
                    self.ldf_ccol_x_contig
                    .slice(offset, batch_size)
                    .join(
                        self.ldf_ccol.lazy(),
                        left_on="contig_collection_id",
                        right_on="contig_collection_id",
                        how="inner"
                    )
                    .select("contig_collection_id", "hash", "contig_id")
                    .collect()
                )

                g = TransformGraph()
                for row in batch.iter_rows(named=True):
                    node_contig = Node(row['contig_id'], "GenomicContig")
                    node_genome = Node(row['hash'], "Genome")
                    g.add_transform_node(node_contig)
                    g.add_transform_node(node_genome)
                    g.add_transform_edge(node_genome, node_contig, "has_contig", data={'source': 'KBERDL'})
                if len(g.t_nodes) > 0:
                    self.load.load2(g)

                if batch.is_empty():
                    break

                offset += batch_size
                pbar.update(len(batch))

    def etl_features(self):
        batch_size = 1000
        offset = 0

        total = self.ldf_feature.select(pl.len()).collect().item()
        with tqdm(total=total, desc="Processing batches", unit="rows") as pbar:
            while True:
                batch = (
                    self.ldf_feature
                    .slice(offset, batch_size)
                    .collect()
                )

                nodes = []
                for row in batch.iter_rows(named=True):
                    node_data = {
                        'source': 'KBERDL',
                        'feature_type': row['type'],
                        'start': row['start'],
                        'end': row['end'],
                        'strand': row['strand'],
                    }
                    node_feature = Node(row['feature_id'], "GenomicFeature", labels=['GenomicFeatureCDS'],
                                        data=node_data)
                    nodes.append(node_feature)
                if len(nodes) > 0:
                    self.vault.bulk_add_nodes2(nodes)

                if batch.is_empty():
                    break

                offset += batch_size
                pbar.update(len(batch))

    def etl_contig_features_edge(self):
        batch_size = 100000
        offset = 0

        total = self.ldf_contig_x_feature.select(pl.len()).collect().item()

        with tqdm(total=total, desc="Processing batches", unit="rows") as pbar:
            while True:
                batch = (
                    self.ldf_contig_x_feature
                    .slice(offset, batch_size)
                    .collect()
                )

                if batch.is_empty():
                    break

                contig_refs = [['GenomicContig', x] for x in set(batch['contig_id'])]
                feature_refs = [['GenomicFeature', x] for x in set(batch['feature_id'])]
                query_features = self.vault.query_eid(feature_refs)
                query_contigs = self.vault.query_eid(contig_refs)
                contig_eid = {o['key']: o['elementId'] for o in query_contigs}
                feature_eid = {o['key']: o['elementId'] for o in query_features}
                payload = []
                for row in batch.iter_rows(named=True):
                    payload.append([
                        contig_eid[row['contig_id']],
                        'has_feature',
                        feature_eid[row['feature_id']],
                        {'source': 'KBERDL'}
                    ])
                if len(payload) > 0:
                    self.vault.bulk_add_edges2(payload)

                offset += batch_size
                pbar.update(len(batch))

    def etl_protein(self, ldf_protein):
        """
        ldf_protein: NR Protein Parquet
        """
        batch_size = 20000
        offset = 0

        total = ldf_protein.select(pl.len()).collect().item()

        with tqdm(total=total, desc="Processing batches", unit="rows") as pbar:
            while True:
                batch = (
                    ldf_protein
                    .slice(offset, batch_size)
                    .select(["protein_id", "hash", "length"])
                    .collect()
                )

                nodes = []
                for row in batch.iter_rows(named=True):
                    node_data = {
                        'cdm_protein_id': row['protein_id'],
                        'length': row['length'],
                    }
                    node = Node(row['hash'], "ProteinSequence", data=node_data)
                    nodes.append(node)

                if len(nodes) > 0:
                    self.vault.bulk_add_nodes2(nodes)
                if batch.is_empty():
                    break

                offset += batch_size
                pbar.update(len(batch))

    def etl_feature_proteins_edges(self):
        batch_size = 100000
        offset = 0

        total = self.ldf_feature_x_protein.select(pl.len()).collect().item()

        with tqdm(total=total, desc="Processing batches feature_x_protein", unit="rows") as pbar:
            while True:
                batch = (
                    self.ldf_feature_x_protein
                    .slice(offset, batch_size)
                    .join(
                        self.ldf_protein.lazy(),
                        left_on="protein_id",
                        right_on="protein_id",
                        how="inner"
                    )
                    .select(["feature_id", "hash"])
                    .collect()
                )
                if batch.is_empty():
                    break
                protein_refs = [['ProteinSequence', x] for x in set(batch['hash'])]
                feature_refs = [['GenomicFeature', x] for x in set(batch['feature_id'])]
                query_features = self.vault.query_eid(feature_refs)
                query_proteins = self.vault.query_eid(protein_refs)
                protein_eid = {o['key']: o['elementId'] for o in query_proteins}
                feature_eid = {o['key']: o['elementId'] for o in query_features}
                payload = []
                for row in batch.iter_rows(named=True):
                    payload.append([
                        feature_eid[row['feature_id']],
                        'has_protein_sequence',
                        protein_eid[row['hash']],
                        {'source': 'KBERDL'}
                    ])
                if len(payload) > 0:
                    self.vault.bulk_add_edges2(payload)

                offset += batch_size
                pbar.update(len(batch))

    def etl_clusters(self, base_dir=Path('/home/fliu/GSP2026/reference_data/berdl_db/ke-pangenomes')):
        parquet_dir = Path(f"{base_dir}/parquet/table_gene_cluster_V1.0")
        lf_gene_cluster = pl.scan_parquet(str(parquet_dir / "*.parquet"))

        batch_size = 200000
        offset = 0
        total = lf_gene_cluster.select(pl.len()).collect().item()
        with tqdm(total=total, desc=f"Processing BERLD KE clusters", unit="rows") as pbar:
            while True:
                batch = (
                    lf_gene_cluster
                    .slice(offset, batch_size)
                    .collect()
                )
                if batch.is_empty():
                    break
                payload = []
                for row in batch.iter_rows(named=True):
                    node = Node(row['gene_cluster_id'], 'ClusterGenomicFeature',
                                labels=['ClusterPangenome'],
                                data={
                                    'mOTUpan_core': row['is_core'] == '1',
                                    'gtdb_species_clade_id': row['gtdb_species_clade_id']
                                })
                    payload.append(node)
                if len(payload) > 0:
                    self.vault.bulk_add_nodes2(payload)

                offset += batch_size
                pbar.update(len(batch))
        #Processing BERLD KE clusters: 100%|██████████| 105720840/105720840 [1:01:38<00:00, 28581.36rows/s]


class KBERLKEGG:

    def __init__(self, vault):
        self.vault = vault
        pass

    def aaa(self):
        ldf_bakta_ko = pl.read_parquet("/home/fliu/scratch/data/CDM/bakta_nodes_KEGGOrtholog.parquet")
        ldf_ko_annotation = pl.read_parquet("/home/fliu/scratch/data/CDM/ko_to_annotation.parquet")
        kos1 = set(ldf_ko_annotation['KEGGOrtholog'])
        kos2 = set(ldf_bakta_ko['key'])
        all_kos = kos1 | kos2

        payload = []
        for ko_id in all_kos:
            node = Node(ko_id, 'KEGGOrtholog')
            payload.append(node)
        len(payload)
        self.vault.bulk_add_nodes2(payload)

        functional_annotation = set(ldf_ko_annotation['FunctionalAnnotation'])
        import hashlib
        def _hash_value(value: str) -> str:
            return hashlib.sha256(value.encode()).hexdigest()

        def get_product(_product):
            node_product = Node(_hash_value(_product), 'FunctionalAnnotation', data={'_value': _product})
            return node_product

        func_to_node = {}
        for func_a in functional_annotation:
            func_to_node[func_a] = get_product(func_a)
        h_to_node = {n.key: n for n in func_to_node.values()}
        refs = [[x.primary_label, x.key] for x in func_to_node.values()]
        query_func_a = self.vault.query_eid(refs)

        payload = []
        for o in query_func_a:
            if o['elementId'] is None:
                payload.append(h_to_node[o['key']])
        self.vault.bulk_add_nodes2(payload)

        ko_eids = self.vault.query_eid([['KEGGOrtholog', x] for x in ldf_ko_annotation['KEGGOrtholog']])
        ko_to_eids = {o['key']: o['elementId'] for o in ko_eids}
        fa_eids = self.vault.query_eid(
            [['FunctionalAnnotation', func_to_node[x].key] for x in ldf_ko_annotation['FunctionalAnnotation']])
        fa_to_eids = {o['key']: o['elementId'] for o in fa_eids}

        payload = []
        for row in ldf_ko_annotation.iter_rows(named=True):
            edge = [
                ko_to_eids[row['KEGGOrtholog']],
                'has_annotation_event',
                fa_to_eids[func_to_node[row['FunctionalAnnotation']].key],
                {"agent": "KEGG", "version": row['version']}
            ]
            payload.append(edge)
        self.vault.bulk_add_edges2(payload)

        return ko_to_eids

    def bbb(self, ko_to_eids: dict):
        df_ko_gene = pl.read_parquet("/home/fliu/scratch/data/CDM/ko_to_gene.parquet")
        df_bakta_gene = pl.read_parquet("/home/fliu/scratch/data/CDM/bakta_nodes_IdentifierGene.parquet")
        igene1 = set(df_ko_gene['IdentifierGene'])
        igene2 = set(df_bakta_gene['key'])
        all_genes = igene1 | igene2

        payload = []
        for _id in all_genes:
            node = Node(_id, 'IdentifierGene')
            payload.append(node)
        len(payload)
        gene_eids = self.vault.bulk_add_nodes2(payload)
        # pass
        payload = []
        for row in df_ko_gene.iter_rows(named=True):
            identifier = row['IdentifierGene']
            identifier = identifier.replace(' ', '_')
            edge = [
                ko_to_eids[row['KEGGOrtholog']],
                'has_annotation_event',
                gene_eids[f"IdentifierGene/{identifier}"],
                {"agent": "KEGG", "version": row['version']}
            ]
            payload.append(edge)
        self.vault.bulk_add_edges2(payload)


class KBERDLBakta:

    def __init__(self):
        # ldf_bakta_knowledge = pl.scan_parquet('/home/fliu/scratch/data/CDM/bakta_edges/*.parquet')
        self.ldf_bakta_knowledge = None
        self.vault = None

    def very_dumb_method_that_reads_entire_collection(self, collection_bakta):
        from collections import Counter
        i = 0
        cursor = collection_bakta.find()
        ct_psc = Counter()
        ct_pscc = Counter()
        ct_ips = Counter()
        ct_ups = Counter()
        ct_doc = Counter()
        all_nodes = {}
        all_edges = []
        for doc in tqdm(cursor):
            ct_doc.update({k: 1 for k in doc})
            ct_psc.update({k: 1 for k in doc.get('psc', {})})
            ct_pscc.update({k: 1 for k in doc.get('pscc', {})})
            ct_ips.update({k: 1 for k in doc.get('ips', {})})
            ct_ups.update({k: 1 for k in doc.get('ups', {})})
            nodes, edges = process_bakta(doc)
            all_nodes.update(nodes)
            all_edges += edges
            i += 1
            # if i >= 1000000:
            #    break
        from collections import Counter
        ct = Counter()
        for k in all_nodes:
            ct.update({k.split('/')[0]: 1})
        ct

    def validate_ec(ec: str) -> bool:
        import re
        if not re.match(r'^[1-7]\.([\dn]+|-)\.([\dn]+|-)\.([\dn]+|-)$', ec):
            return False
        parts = ec.split('.')
        seen_dash = False
        for part in parts:
            if part == '-':
                seen_dash = True
            elif seen_dash:
                return False
        return True

    @staticmethod
    def _hash_value(value: str) -> str:
        import hashlib
        return hashlib.sha256(value.encode()).hexdigest()

    @staticmethod
    def get_product(d):
        _product = d.get('product')
        if _product:
            node_product = Node(_hash_value(_product), 'FunctionalAnnotation', data={'_value': _product})
            return node_product

    @staticmethod
    def process_bakta(d):
        def _hash_value(value: str) -> str:
            import hashlib
            return hashlib.sha256(value.encode()).hexdigest()
        def build_edge_data(d):
            _data = {"agent": "tool", "method": "bakta"}
            _allowed_types = (str, bool, int, float)
            for k, v in d.items():
                if v is not None:
                    if isinstance(v, list):
                        _data[f"_bakta_{k}"] = '; '.join(str(i) for i in v)
                    elif isinstance(v, _allowed_types):
                        _data[f"_bakta_{k}"] = v
                    else:
                        raise TypeError(f"unsupported type for key '{k}': {type(v).__name__}")
            return _data


        nodes = {}
        edges = []

        node_protein = Node(d['_id'], 'ProteinSequence')
        nodes[node_protein.id] = node_protein

        node_product = KBERDLBakta.get_product(d)
        if node_product:
            nodes[node_product.id] = node_product
            edges.append([node_protein.id,
                          'has_annotation_event',
                          node_product.id,
                          {"agent": "tool", "method": "bakta"}])
        else:
            _product = 'none'
            node_product = Node(_hash_value(_product), 'FunctionalAnnotation', data={'_value': _product})
            nodes[node_product.id] = node_product
            edges.append([node_protein.id,
                          'has_annotation_event',
                          node_product.id,
                          {"agent": "tool", "method": "bakta"}])

        if 'ups' in d:
            node_ups = Node(d['ups']['uniparc_id'], 'UniParc')
            nodes[node_ups.id] = node_ups
            edges.append([node_protein.id,
                          'has_annotation_event',
                          node_ups.id,
                          {"agent": "tool", "method": "bakta"}])
        if 'ips' in d:
            node_ips_u100 = Node(d['ips']['uniref100_id'], 'UniProtUniRef100')
            nodes[node_ips_u100.id] = node_ips_u100
            edges.append([node_protein.id,
                          'has_annotation_event',
                          node_ips_u100.id,
                          {"agent": "tool", "method": "bakta"}])
            # node_product = get_product(doc)
            # if node_product:
            #    nodes[node_product.id] = node_product
        if 'psc' in d:
            psc = d['psc']
            _data = build_edge_data(psc)
            node_psc_u90 = Node(psc['uniref90_id'], 'UniProtUniRef90')
            nodes[node_psc_u90.id] = node_psc_u90
            edges.append([node_protein.id,
                          'has_annotation_event',
                          node_psc_u90.id,
                          _data])
            # node_product = get_product(doc)
            # if node_product:
            #    nodes[node_product.id] = node_product
            if 'ec_ids' in psc:
                for _ec in psc['ec_ids']:
                    node_ec = Node(_ec, 'ECNumber')
                    nodes[node_ec.id] = node_ec
                    edges.append([node_protein.id,
                                  'has_annotation_event',
                                  node_ec.id,
                                  {"agent": "tool", "method": "bakta"}])
            kegg_orthology_id = psc.get('kegg_orthology_id')
            if kegg_orthology_id is not None:
                node_ko = Node(kegg_orthology_id, 'KEGGOrtholog')
                nodes[node_ko.id] = node_ko
                edges.append([node_protein.id,
                              'has_annotation_event',
                              node_ko.id,
                              {"agent": "tool", "method": "bakta"}])
            # go_ids
            go_ids = psc.get('go_ids')
            if go_ids is not None:
                for _go_id in go_ids:
                    node_go = Node(_go_id, 'OntologyGO')
                    nodes[node_go.id] = node_go
                    edges.append([node_protein.id,
                                  'has_annotation_event',
                                  node_go.id,
                                  {"agent": "tool", "method": "bakta"}])
        if 'pscc' in d:
            pscc = d['pscc']
            _data = build_edge_data(pscc)
            node_pscc_u50 = Node(pscc['uniref50_id'], 'UniProtUniRef50')
            nodes[node_pscc_u50.id] = node_pscc_u50
            # node_product = get_product(doc)
            # if node_product:
            #    nodes[node_product.id] = node_product
            edges.append([node_protein.id,
                          'has_annotation_event',
                          node_pscc_u50.id,
                          _data])

        if 'expert' in d:
            for _expert_a in d['expert']:
                e_product = _expert_a['product']
                node_product = Node(_hash_value(e_product), 'FunctionalAnnotation', data={'_value': e_product})
                nodes[node_product.id] = node_product
                _data = build_edge_data(_expert_a)
                edges.append([node_protein.id,
                              'has_annotation_event',
                              node_product.id,
                              _data])
                # print(_data)

        genes = d.get('genes')
        if genes is not None:
            for _gene_str in genes:
                node_gene = Node(_gene_str, 'IdentifierGene')
                nodes[node_gene.id] = node_gene
                edges.append([node_protein.id,
                              'has_annotation_event',
                              node_gene.id,
                              {"agent": "tool", "method": "bakta"}])
        return nodes, edges

    @staticmethod
    def load_bakta_output_json_to_mongo(bakta_output_file, collection_bakta):
        with open(bakta_output_file, 'r') as fh:
            bakta_results = json.load(fh)

            id_to_result = {o['id']: o for o in bakta_results['features']}
            cursor = collection_bakta.find(
                {"_id": {"$in": list(id_to_result)}},
                {"_id": 1})
            found_ids = {doc['_id'] for doc in cursor}
            from pymongo import InsertOne
            operations = []
            for h in id_to_result:
                if h not in found_ids:
                    d = dict(id_to_result[h])
                    del d['id']
                    d['_id'] = h
                    operations.append(InsertOne(d))
            print(len(operations))
            if len(operations) > 0:
                collection_bakta.bulk_write(operations, ordered=False)

    """
    def aaa(self, doc):
        if 'aa' in doc:
            del doc['aa']
        if 'description' in doc:
            del doc['description']
        if 'locus' in doc:
            del doc['locus']

    def enforce_schema_null_missing(df: pl.DataFrame, schema: dict[str, pl.DataType]) -> pl.DataFrame:
        # 1) add any missing columns as typed nulls
        for name, dtype in schema.items():
            if name not in df.columns:
                df = df.with_columns(pl.lit(None, dtype=dtype).alias(name))

        # 2) (optional) drop extras + order columns consistently
        df = df.select(list(schema.keys()))

        # 3) cast columns to the target dtypes
        return df.cast(schema, strict=False)
    """

    def edges(self):
        batch_size = 200000
        offset = 0
        not_found = set()
        total = self.ldf_bakta_knowledge.select(pl.len()).collect().item()
        with tqdm(total=total, desc=f"Processing Bakta Edges", unit="rows") as pbar:
            while True:
                batch = (
                    self.ldf_bakta_knowledge
                    .slice(offset, batch_size)
                    .collect()
                )
                if batch.is_empty():
                    break
                u_refs = set()
                for lkp in batch['src']:
                    label, key = lkp.split('/', 1)
                    u_refs.add((label, key))
                for lkp in batch['dst']:
                    label, key = lkp.split('/', 1)
                    if label == 'GOTerm':
                        label = 'OntologyGO'
                    u_refs.add((label, key))
                eids = self.vault.query_eid(u_refs)
                lkp_to_eid = {}
                for o in eids:
                    if o['elementId']:
                        lkp_to_eid[f"{o['type']}/{o['key']}"] = o['elementId']
                payload = []
                for row in batch.iter_rows(named=True):
                    lkp_src = fix_label(row['src'])
                    lkp_dst = fix_label(row['dst'])
                    src_eid = lkp_to_eid.get(lkp_src)
                    dst_eid = lkp_to_eid.get(lkp_dst)
                    if src_eid is None:
                        not_found.add(row['src'])
                    if dst_eid is None:
                        not_found.add(row['dst'])
                    if src_eid is not None and dst_eid is not None:
                        edge = [src_eid, row['label'], dst_eid, json.loads(row['data'])]
                        payload.append(edge)
                if len(payload) > 0:
                    self.vault.bulk_add_edges2(payload)

                offset += batch_size
                pbar.update(len(batch))


class KBERDLKoFam:

    def __init__(self, path_out: Path):
        self.path_out = path_out

    @staticmethod
    def process_dir(out_dir, path_out: Path):
        if not (out_dir.name.startswith('ke_nr_') and out_dir.is_dir()):
            return None

        result = out_dir / 'output'
        tabular = out_dir / 'tabular' / 'tabular.txt'
        outfile_name = path_out / f"{out_dir.name}.parquet"

        if not outfile_name.exists() and result.exists() and tabular.exists():
            df_out = pl.read_csv(result, has_header=False, separator='\t', new_columns=["h", "KO"])
            parser = KoFamParser()
            tabular_data = parser.parse(tabular)
            return KBERDLKoFam.write_to_parquet(df_out, tabular_data, outfile_name)

        return None

    def run_executor(self, base: Path, threads=8):
        from concurrent.futures import ThreadPoolExecutor, as_completed

        dirs = list(base.iterdir())

        with ThreadPoolExecutor(max_workers=threads) as executor:
            futures = {executor.submit(KBERDLKoFam.process_dir, out_dir, self.path_out): out_dir for out_dir in dirs}
            for future in tqdm(as_completed(futures), total=len(dirs)):
                try:
                    future.result()
                except Exception as e:
                    print(f"Error processing {futures[future]}: {e}")

    @staticmethod
    def write_to_parquet(df_out, tabular_data, outfile: Path):
        hit_to_ko = {}
        for ko in tabular_data:
            for hit in tabular_data[ko].hits:
                if hit.target_name not in hit_to_ko:
                    hit_to_ko[hit.target_name] = {}
                if ko not in hit_to_ko[hit.target_name]:
                    hit_to_ko[hit.target_name][ko] = hit
                else:
                    raise ValueError('!')

        data = {
            'src': [],
            'dst': [],
            'type': [],
            'data': [],
        }
        for h, ko in df_out.iter_rows():
            if ko is not None:
                if len(ko) == 6:
                    # hits = tabular_data[ko].hits
                    # hit = [hit for hit in hits if hit.target_name == h][0]
                    hit = hit_to_ko[h][ko]
                    _data = {
                        'full_sequence_evalue': hit.full_sequence.evalue,
                        'full_sequence_score': hit.full_sequence.score,
                        'full_sequence_bias': hit.full_sequence.bias,
                        'best_domain_evalue': hit.best_domain.evalue,
                        'best_domain_score': hit.best_domain.score,
                        'best_domain_bias': hit.best_domain.bias,

                        'domain_est_exp': hit.domain_estimation.exp,
                        'domain_est_reg': hit.domain_estimation.reg,
                        'domain_est_clu': hit.domain_estimation.clu,
                        'domain_est_ov': hit.domain_estimation.ov,
                        'domain_est_env': hit.domain_estimation.env,
                        'domain_est_dom': hit.domain_estimation.dom,
                        'domain_est_rep': hit.domain_estimation.rep,
                        'domain_est_inc': hit.domain_estimation.inc,
                    }
                    data["src"].append(h)
                    data["dst"].append(ko)
                    data["type"].append('has_annotation_event')
                    data["data"].append(_data)
                else:
                    raise ValueError('!')
        table = pa.table(data)
        pq.write_table(table, outfile)


class KBERDLRast:

    def __init__(self, client):
        self.client = client
        self.collection_rast = None

    def run(self, ldf):
        batch_size = 20000
        offset = 0
        total = ldf.select(pl.len()).collect().item()
        total_processed = 0
        baddies = {}
        with tqdm(total=total, desc="Processing batches", unit="rows") as pbar:
            while True:
                batch = ldf.slice(offset, batch_size).collect()

                if batch.height == 0:
                    break

                proteins = {}
                for row in batch.iter_rows(named=True):
                    seq = row['sequence']
                    stype = classify_protein_sequence(row['sequence'])
                    if stype == SeqType.BASIC or stype == SeqType.EXTENDED:
                        protein = ProteinSequence(row['sequence'])
                        h = protein.hash_value
                        if h not in proteins:
                            proteins[h] = seq
                    else:
                        baddies[row['hash']] = seq
                    total_processed += 1

                batch_h = set(proteins)
                cursor = self.collection_rast.find(
                    {"_id": {"$in": list(batch_h)}},
                    {"_id": 1})
                found_ids = {doc['_id'] for doc in cursor}

                proteins_filter = {k: seq for k, seq in proteins.items() if k not in found_ids}

                self.aaa(proteins_filter)

                offset += batch_size
                pbar.update(len(batch))

    def aaa(self, proteins: dict):
        p_features = [{"id": k, "protein_translation": seq} for k, seq in proteins.items()]
        annotation_rast = self.client.f(p_features)
        h_to_result = {o['id']: o for o in annotation_rast[0]['features']}

        self.cache_to_mongo(h_to_result)

    @staticmethod
    def process_doc(result):
        doc = {}
        for k in result:
            if k == 'id':
                doc['_id'] = result[k]
            else:
                if k != 'protein_translation':
                    doc[k] = result[k]
        return doc

    def cache_to_mongo(self, h_to_result):
        from pymongo import ReplaceOne

        operations = []

        for h in h_to_result:
            result = h_to_result.get(h)
            if result is not None:
                operations.append(ReplaceOne(
                    {"_id": h},  # match condition
                    KBERDLRast.process_doc(result),  # full replacement document
                    upsert=True  # insert if not exists
                ))
            else:
                operations.append(ReplaceOne(
                    {"_id": h},  # match condition
                    {"_id": h, 'function': None},  # full replacement document
                    upsert=True  # insert if not exists
                ))
        return self.collection_rast.bulk_write(operations, ordered=False)


class KBERDLGo:

    def __init__(self):
        self.vault = None

    def etl(self, df_go_ontology):
        from modelseed_vault.core.transform_graph import Node
        payload = []
        ns_map = {
            'molecular_function': 'OntologyGOMolecularFunction',
            'biological_process': 'OntologyGOBiologicalProcess',
            'cellular_component': 'OntologyGOCellularComponent',
        }
        for row in df_go_ontology.iter_rows(named=True):
            node = Node(row['go_id'].replace('_', ':'), 'OntologyGO')
            labels = []
            if row['meta_deprecated']:
                labels.append('Obsolete')
            if row['meta_bp_namespace']:
                labels.append(ns_map[row['meta_bp_namespace']])
            if len(labels) > 0:
                node.labels = labels
            _data = {}
            if row['label']:
                _data['label'] = row['label']
            if row['meta_comments']:
                _data['comments'] = row['meta_comments']
            if row['meta_bp_created_by']:
                _data['meta_bp_created_by'] = row['meta_bp_created_by']
            if row['meta_definition']:
                _data['definition'] = row['meta_definition']
            if len(_data) > 0:
                node.data = _data
            payload.append(node)
        self.vault.bulk_add_nodes2(payload)
