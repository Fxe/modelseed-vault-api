from __future__ import annotations

from typing import Iterable, Optional
import json
from pathlib import Path
import networkx as nx
from collections import defaultdict


def collect_triple(center: str, edges: list) -> tuple:
    upstream = None    # closest neighbor with end <= center_start
    downstream = None  # closest neighbor with start >= center_end
    upstream_dist = float('inf')
    downstream_dist = float('inf')
    center_start = center_end = None

    for e in edges:
        src = e['source']['entry']
        tgt = e['target']['entry']
        p = e['properties']

        if src == center:
            # center is source; neighbor is target
            cs, ce = p['source_start'], p['source_end']
            ns, ne = p['target_start'], p['target_end']
            neighbor = tgt
        elif tgt == center:
            # center is target; neighbor is source
            cs, ce = p['target_start'], p['target_end']
            ns, ne = p['source_start'], p['source_end']
            neighbor = src
        else:
            continue

        center_start, center_end = cs, ce

        if ne <= cs:  # neighbor is upstream
            d = cs - ne
            if d < upstream_dist:
                upstream_dist = d
                upstream = neighbor
        elif ns >= ce:  # neighbor is downstream
            d = ns - ce
            if d < downstream_dist:
                downstream_dist = d
                downstream = neighbor

    return upstream, center, downstream


class QuintuplePairBuilder:
    """
    Build fc_adj_pairs and fc_adj_5tuples from per-genome ec + chain files using
    contig-level window reconstruction.

    Unlike TriplePairBuilder (which uses per-feature nearest-neighbour lookups),
    this class uses build_windows_from_edges to reconstruct full contig gene order
    and slides radius-1 (triple) and radius-2 (quintuple) windows in the correct
    5'→3' transcriptional direction.

    Parameters
    ----------
    gap_filter : int or None
        Maximum intergenic gap (bp) to allow inside a window.  Set ~200 for
        strict operon membership; None to keep all adjacent pairs.
    strand_filter : bool
        If True, only emit windows where all genes share the same strand.

    Attributes
    ----------
    fc_adj_pairs : dict[(f_up, f_center, f_down), list[genome_id]]
        Radius-1 windows in transcriptional order.
    fc_adj_5tuples : dict[(f1, f_x, f_center, f_y, f2), list[genome_id]]
        Radius-2 windows in transcriptional order.
    """

    def __init__(self, gap_filter=None, strand_filter=True):
        self.gap_filter = gap_filter
        self.strand_filter = strand_filter
        self.fc_adj_pairs: dict = {}
        self.fc_adj_5tuples: dict = {}

    def add_genome(self, genome_id, feature_chain, feature_ec) -> None:
        """Ingest one genome's ec + chain files."""
        genome_h = genome_id

        edges = [o for o in feature_chain if o['type'] == 'ADJACENT_TO']
        triples, quintuples = build_windows_from_edges(
            edges,
            gap_filter=self.gap_filter,
            strand_filter=self.strand_filter,
        )

        def _fn(fid):
            v = feature_ec.get(fid)
            return v[1] if v else None

        for triple in triples:
            key = tuple(_fn(fid) for fid in triple)
            if None in key:
                continue
            if key in self.fc_adj_pairs:
                self.fc_adj_pairs[key].append(genome_h)
            else:
                self.fc_adj_pairs[key] = [genome_h]

        for quint in quintuples:
            key = tuple(_fn(fid) for fid in quint)
            if None in key:
                continue
            if key in self.fc_adj_5tuples:
                self.fc_adj_5tuples[key].append(genome_h)
            else:
                self.fc_adj_5tuples[key] = [genome_h]

    # ------------------------------------------------------------------
    def add(self, file_ec: Path, file_chain: Path) -> None:
        """Ingest one genome's ec + chain files."""
        genome_h = file_ec.stem
        with open(file_ec) as fh:
            feature_ec = json.load(fh)
        with open(file_chain) as fh:
            feature_chain = json.load(fh)

        self.add_genome(genome_h, feature_chain, feature_ec)


class TriplePairBuilder:

    def __init__(self):
        self.fc_adj_pairs = {}

    def add(self, file_ec: Path, file_chain: Path):
        genome_h = file_ec.stem
        with open(file_ec, 'r') as fh:
            feature_ec = json.load(fh)
        with open(file_chain, 'r') as fh:
            feature_chain = json.load(fh)
        feature_capture = set()
        feature_triple = {}
        for o in feature_chain:
            if o['type'] == 'ADJACENT_TO':
                a = o['source']['entry']
                b = o['target']['entry']
                if a not in feature_triple:
                    feature_triple[a] = [o]
                else:
                    feature_triple[a].append(o)
                if b not in feature_triple:
                    feature_triple[b] = [o]
                else:
                    feature_triple[b].append(o)
        for feature_id in feature_triple:
            triple = collect_triple(feature_id, feature_triple[feature_id])
            if None not in triple:
                t0 = feature_ec.get(triple[0], [None, None])[1]
                t1 = feature_ec.get(triple[1], [None, None])[1]
                t2 = feature_ec.get(triple[2], [None, None])[1]
                f_triple = (t0, t1, t2)
                f_triple_rev = (t2, t1, t0)
                if f_triple not in self.fc_adj_pairs and f_triple_rev not in self.fc_adj_pairs:
                    self.fc_adj_pairs[f_triple] = [genome_h]
                elif f_triple in self.fc_adj_pairs:
                    self.fc_adj_pairs[f_triple].append(genome_h)
                elif f_triple_rev in self.fc_adj_pairs:
                    self.fc_adj_pairs[f_triple_rev].append(genome_h)
                else:
                    raise ValueError('!')


class NeighborhoodEstimator:
    """
    Hierarchical estimator for P(X, Y | F) where F is a center gene function
    and X, Y are its upstream and downstream neighbors on a contig.

    Hierarchy (in order of preference):
        1. Radius-2: P(X, Y | F, F1, F2) -- requires 5-tuple support
        2. Radius-2 marginal: P(X, Y | F)  -- joint, from 5-tuples
        3. Radius-1: P(X | F) and P(Y | F) -- independent, from triples
        4. Uniform fallback if F is unseen.

    Conventions:
        - Triples are (upstream, center, downstream), ordered along the contig.
        - 5-tuples are (F1, X, F, Y, F2), ordered along the contig.
        - "Distinct genomes" is the unit of support: a genome contributes
          at most once per window key.
    """

    def __init__(
        self,
        total_genomes: int,
        min_support_r2: int = 20,
        min_support_r1: int = 5,
    ):
        self.total_genomes = total_genomes
        self.min_support_r2 = min_support_r2
        self.min_support_r1 = min_support_r1

        # Radius-1 storage: (center) -> {(X, Y): genome_set}
        self._r1_joint: dict[str, dict[tuple[str, str], set]] = defaultdict(lambda: defaultdict(set))
        # Radius-1 marginals: (center, 'up') -> {X: genome_set}, (center, 'down') -> {Y: genome_set}
        self._r1_up: dict[str, dict[str, set]] = defaultdict(lambda: defaultdict(set))
        self._r1_down: dict[str, dict[str, set]] = defaultdict(lambda: defaultdict(set))
        # Genomes where each center was seen
        self._r1_center_genomes: dict[str, set] = defaultdict(set)

        # Radius-2 storage: (F, F1, F2) -> {(X, Y): genome_set}
        self._r2_joint: dict[tuple[str, str, str], dict[tuple[str, str], set]] = defaultdict(lambda: defaultdict(set))
        # Genomes per (F, F1, F2) flanking context
        self._r2_context_genomes: dict[tuple[str, str, str], set] = defaultdict(set)

        # Radius-2 marginal (no flanking conditioning): F -> {(X, Y): genome_set}
        self._r2_marginal: dict[str, dict[tuple[str, str], set]] = defaultdict(lambda: defaultdict(set))
        self._r2_center_genomes: dict[str, set] = defaultdict(set)

    # ---------- ingestion ----------

    def add_triples(self, fc_adj_pairs: dict[tuple[str, str, str], Iterable[str]]) -> None:
        """Ingest radius-1 windows: keys are (upstream, center, downstream)."""
        for (up, center, down), genomes in fc_adj_pairs.items():
            gset = set(genomes)
            self._r1_joint[center][(up, down)].update(gset)
            self._r1_up[center][up].update(gset)
            self._r1_down[center][down].update(gset)
            self._r1_center_genomes[center].update(gset)

    def add_quintuples(self, fc_adj_5tuples: dict[tuple[str, str, str, str, str], Iterable[str]]) -> None:
        """Ingest radius-2 windows: keys are (F1, X, F, Y, F2)."""
        for (f1, x, center, y, f2), genomes in fc_adj_5tuples.items():
            gset = set(genomes)
            self._r2_joint[(center, f1, f2)][(x, y)].update(gset)
            self._r2_context_genomes[(center, f1, f2)].update(gset)
            self._r2_marginal[center][(x, y)].update(gset)
            self._r2_center_genomes[center].update(gset)

    # ---------- core query ----------

    def joint(
        self,
        center: str,
        f1: Optional[str] = None,
        f2: Optional[str] = None,
        top_k: Optional[int] = None,
    ) -> list[tuple[tuple[str, str], float, int, str]]:
        """
        Return P(X, Y | F [, F1, F2]) as a ranked list.

        Each row: ((X, Y), probability, support, source_tier)
        source_tier ∈ {'r2_conditional', 'r2_marginal', 'r1_independent', 'unseen'}
        """
        # Tier 1: full radius-2 conditional, if flanks provided and support is sufficient.
        if f1 is not None and f2 is not None:
            ctx = (center, f1, f2)
            denom = len(self._r2_context_genomes.get(ctx, ()))
            if denom >= self.min_support_r2:
                rows = [
                    (xy, len(g) / denom, len(g), 'r2_conditional')
                    for xy, g in self._r2_joint[ctx].items()
                ]
                rows.sort(key=lambda r: -r[1])
                return rows[:top_k] if top_k else rows

        # Tier 2: radius-2 marginal joint P(X, Y | F).
        denom = len(self._r2_center_genomes.get(center, ()))
        if denom >= self.min_support_r2:
            rows = [
                (xy, len(g) / denom, len(g), 'r2_marginal')
                for xy, g in self._r2_marginal[center].items()
            ]
            rows.sort(key=lambda r: -r[1])
            return rows[:top_k] if top_k else rows

        # Tier 3: radius-1 independent factorization.
        denom = len(self._r1_center_genomes.get(center, ()))
        if denom >= self.min_support_r1:
            up_probs = {x: len(g) / denom for x, g in self._r1_up[center].items()}
            down_probs = {y: len(g) / denom for y, g in self._r1_down[center].items()}
            rows = []
            for x, px in up_probs.items():
                for y, py in down_probs.items():
                    # Joint support approximated by intersection (more accurate than px*py*denom)
                    g_up = self._r1_up[center][x]
                    g_dn = self._r1_down[center][y]
                    support = len(g_up & g_dn)
                    rows.append(((x, y), px * py, support, 'r1_independent'))
            rows.sort(key=lambda r: -r[1])
            return rows[:top_k] if top_k else rows

        # Tier 4: unseen.
        return [(('?', '?'), 0.0, 0, 'unseen')]

    def marginal_upstream(self, center: str, top_k: Optional[int] = None):
        """P(X | F) from radius-1 data."""
        denom = len(self._r1_center_genomes.get(center, ()))
        if denom < self.min_support_r1:
            return []
        rows = [(x, len(g) / denom, len(g)) for x, g in self._r1_up[center].items()]
        rows.sort(key=lambda r: -r[1])
        return rows[:top_k] if top_k else rows

    def marginal_downstream(self, center: str, top_k: Optional[int] = None):
        """P(Y | F) from radius-1 data."""
        denom = len(self._r1_center_genomes.get(center, ()))
        if denom < self.min_support_r1:
            return []
        rows = [(y, len(g) / denom, len(g)) for y, g in self._r1_down[center].items()]
        rows.sort(key=lambda r: -r[1])
        return rows[:top_k] if top_k else rows

    # ---------- diagnostics ----------

    def support_summary(self, center: str) -> dict:
        """How much evidence backs predictions for this center?"""
        return {
            'r1_center_genomes': len(self._r1_center_genomes.get(center, ())),
            'r1_distinct_pairs': len(self._r1_joint.get(center, {})),
            'r2_center_genomes': len(self._r2_center_genomes.get(center, ())),
            'r2_distinct_pairs': len(self._r2_marginal.get(center, {})),
            'r2_distinct_contexts': sum(
                1 for (c, _, _) in self._r2_context_genomes if c == center
            ),
        }


def build_windows_from_edges(edges, gap_filter=None, strand_filter=True):
    """
    From a list of ADJACENT_TO edges, reconstruct gene order per contig
    and emit radius-1 (triples) and radius-2 (quintuples) windows.

    Parameters
    ----------
    edges : iterable of edge dicts (like your input list).
    gap_filter : optional int, max allowed gap to consider adjacency 'tight'.
    strand_filter : if True, only emit windows where all genes are on the same strand.

    Returns
    -------
    (triples, quintuples) where each is a list of feature-id tuples ordered
    along the contig in transcriptional direction (5' -> 3').
    """
    # Group by contig and reconstruct linear order
    by_contig = defaultdict(list)
    for e in edges:
        by_contig[e['contig']].append(e)

    triples, quintuples = [], []

    for contig, contig_edges in by_contig.items():
        # Build a per-feature record: (feat_id, start, end, strand)
        feats = {}
        for e in contig_edges:
            p = e['properties']
            feats[e['source']['entry']] = (p['source_start'], p['source_end'], p['source_strand'])
            feats[e['target']['entry']] = (p['target_start'], p['target_end'], p['target_strand'])

        # Sort by start coordinate
        ordered = sorted(feats.items(), key=lambda kv: kv[1][0])

        # Optionally split into strand-homogeneous runs
        runs = []
        if strand_filter:
            current = []
            current_strand = None
            for fid, (s, e_, strand) in ordered:
                if strand != current_strand and current:
                    runs.append(current)
                    current = []
                current.append((fid, s, e_, strand))
                current_strand = strand
            if current:
                runs.append(current)
        else:
            runs = [[(fid, s, e_, st) for fid, (s, e_, st) in ordered]]

        # Optionally filter by gap
        if gap_filter is not None:
            new_runs = []
            for run in runs:
                cur = [run[0]]
                for i in range(1, len(run)):
                    gap = run[i][1] - run[i - 1][2]
                    if gap <= gap_filter:
                        cur.append(run[i])
                    else:
                        if len(cur) >= 3:
                            new_runs.append(cur)
                        cur = [run[i]]
                if len(cur) >= 3:
                    new_runs.append(cur)
            runs = new_runs

        # Slide windows. Reverse run if on minus strand so order is 5'->3'.
        for run in runs:
            strand = run[0][3]
            ids = [r[0] for r in run]
            if strand == '-':
                ids = ids[::-1]

            for i in range(len(ids) - 2):
                triples.append(tuple(ids[i:i + 3]))
            for i in range(len(ids) - 4):
                quintuples.append(tuple(ids[i:i + 5]))

    return triples, quintuples


def transition_probs(G: nx.DiGraph, function: str, top_k: int | None = None,
                     min_support: int = 0) -> list[tuple[str, float, int]]:
    """
    Return outgoing transition probabilities from `function`.

    Parameters
    ----------
    G : nx.DiGraph
        Graph built by build_transition_graph.
    function : str
        Source node.
    top_k : int or None
        If given, return only the top-k successors by probability.
    min_support : int
        Drop edges supported by fewer than this many distinct genomes.

    Returns
    -------
    list of (target_function, probability, support), sorted by probability desc.
    """
    if function not in G:
        return []

    rows = [
        (tgt, data['weight'], data['support'])
        for _, tgt, data in G.out_edges(function, data=True)
        if data['support'] >= min_support
    ]
    rows.sort(key=lambda r: r[1], reverse=True)
    return rows[:top_k] if top_k else rows


def build_transition_graph(fc_adj_pairs: dict, total_genomes: int) -> nx.DiGraph:
    """
    Build a probability transition graph from functional triples.

    Each triple (f1, f2, f3) contributes two directed edges:
        f1 -> f2  and  f2 -> f3
    Edge weight = P(target | source) estimated as:
        (# distinct genomes supporting source->target) / (# distinct genomes supporting source)

    Parameters
    ----------
    fc_adj_pairs : dict
        Maps (f1, f2, f3) -> list of genome IDs (may contain duplicates).
    total_genomes : int
        Total number of genomes (kept on the graph as metadata; not used in
        the conditional-probability calculation).

    Returns
    -------
    nx.DiGraph
        Nodes are function strings. Edges have:
          - 'weight'        : P(target | source)
          - 'support'       : # distinct genomes with source -> target
          - 'prevalence'    : support / total_genomes  (marginal frequency)
    """
    # Step 1: aggregate distinct genomes per directed pair (source -> target)
    pair_genomes: dict[tuple[str, str], set] = {}
    # And distinct genomes per source node (denominator for P(target|source))
    source_genomes: dict[str, set] = {}

    for (f1, f2, f3), genomes in fc_adj_pairs.items():
        gset = set(genomes)  # dedupe within this triple
        for src, tgt in ((f1, f2), (f2, f3)):
            pair_genomes.setdefault((src, tgt), set()).update(gset)
            source_genomes.setdefault(src, set()).update(gset)

    # Step 2: build the graph
    G = nx.DiGraph()
    G.graph['total_genomes'] = total_genomes

    for (src, tgt), gset in pair_genomes.items():
        support = len(gset)
        src_total = len(source_genomes[src])
        prob = support / src_total if src_total else 0.0
        G.add_edge(
            src, tgt,
            weight=prob,
            support=support,
            prevalence=support / total_genomes if total_genomes else 0.0,
        )

    # Attach per-node marginal (distinct genomes that ever have this function as a source)
    for node, gset in source_genomes.items():
        if node in G:
            G.nodes[node]['source_genome_count'] = len(gset)

    return G


def build_transition_graph2(fc_adj_pairs: dict, total_genomes: int) -> nx.DiGraph:
    """
    Build a directional transition graph from (upstream, center, downstream) triples.

    Each triple contributes ONE forward edge per genome:
        center -> downstream

    This gives a proper Markov kernel: outgoing weights from each node sum to 1
    (over genomes where the center has a tracked downstream neighbor).

    Weights are:
        P(downstream = T | center = F)
            = |distinct genomes with F->T as center->downstream|
            / |distinct genomes where F appears as center of any triple|

    Parameters
    ----------
    fc_adj_pairs : dict
        Maps (upstream, center, downstream) -> list of genome IDs.
        The middle element is treated as the source; the third as the successor.
    total_genomes : int
        Corpus size (stored as metadata; also used for edge 'prevalence').

    Returns
    -------
    nx.DiGraph with edge attributes:
        - 'weight'     : P(downstream | center)
        - 'support'    : # distinct genomes with center -> downstream
        - 'prevalence' : support / total_genomes
      and node attribute:
        - 'center_genome_count' : # distinct genomes where node was a center
    """
    pair_genomes: dict[tuple[str, str], set] = {}
    center_genomes: dict[str, set] = {}

    for (upstream, center, downstream), genomes in fc_adj_pairs.items():
        gset = set(genomes)
        center_genomes.setdefault(center, set()).update(gset)
        pair_genomes.setdefault((center, downstream), set()).update(gset)

    G = nx.DiGraph()
    G.graph['total_genomes'] = total_genomes

    for (src, tgt), gset in pair_genomes.items():
        support = len(gset)
        denom = len(center_genomes[src])
        G.add_edge(
            src, tgt,
            weight=support / denom if denom else 0.0,
            support=support,
            prevalence=support / total_genomes if total_genomes else 0.0,
        )

    for node, gset in center_genomes.items():
        if node in G:
            G.nodes[node]['center_genome_count'] = len(gset)

    return G

