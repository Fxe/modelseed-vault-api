# Gene Neighborhood Probability Estimator — Design Context

This document captures the logic, design decisions, and pitfalls from building a probabilistic estimator over gene functional neighborhoods. Use as context for future Claude Code sessions.

## Problem Statement

Given a center gene function `F`, estimate the likelihood of its neighboring gene functions:

```
X — F — Y
```

Where X is the upstream neighbor and Y is the downstream neighbor (in 5'→3' transcriptional direction). The estimator should predict P(X, Y | F) and ideally condition on wider context P(X, Y | F, F1, F2) when data supports it.

## Input Data Structure

The graph is built from `ADJACENT_TO` edges from a bioinformatics knowledge graph (Neo4j-backed). Each edge has the form:

```python
{
  'type': 'ADJACENT_TO',
  'source': {'entry': 'CDM:feat-...', 'elementId': '...'},
  'target': {'entry': 'CDM:feat-...', 'elementId': '...'},
  'properties': {
    'gap': 957,
    'source_strand': '-',
    'target_strand': '-',
    'source_start': 22571, 'source_end': 23266,
    'target_start': 24224, 'target_end': 24859,
  },
  'contig': 'CDM:cntg-...'
}
```

From these edges, we derive:

- **Triples**: `(upstream, center, downstream)` ordered along the contig in transcriptional direction.
- **Quintuples**: `(F1, X, F, Y, F2)` — same idea, radius 2.

These are aggregated across many genomes into:

- `fc_adj_pairs: dict[(f1, f2, f3) -> list[genome_ids]]`
- `fc_adj_5tuples: dict[(f1, f2, f3, f4, f5) -> list[genome_ids]]`

A genome may appear multiple times in a list; for probability estimation we **dedupe per key** (count each genome at most once per window).

## The `collect_triple` Helper

For a given center feature, find its closest upstream and downstream neighbors by coordinate distance. Features can appear as either source or target of an `ADJACENT_TO` edge — must normalize both directions.

```python
def collect_triple(center: str, edges: list) -> tuple:
    upstream = downstream = None
    upstream_dist = downstream_dist = float('inf')

    for e in edges:
        src, tgt = e['source']['entry'], e['target']['entry']
        p = e['properties']

        if src == center:
            cs, ce = p['source_start'], p['source_end']
            ns, ne = p['target_start'], p['target_end']
            neighbor = tgt
        elif tgt == center:
            cs, ce = p['target_start'], p['target_end']
            ns, ne = p['source_start'], p['source_end']
            neighbor = src
        else:
            continue

        if ne <= cs:  # upstream by coordinate
            d = cs - ne
            if d < upstream_dist:
                upstream_dist, upstream = d, neighbor
        elif ns >= ce:  # downstream by coordinate
            d = ns - ce
            if d < downstream_dist:
                downstream_dist, downstream = d, neighbor

    return (upstream, center, downstream)
```

**Important caveat:** This is *coordinate*-based, not biologically directional. For minus-strand genes, true biological upstream is at higher coordinates. Flip ordering based on strand when building windows for transcriptional analysis.

## Key Pitfall: Naive Markov Construction

### The bad first attempt

Initial graph construction decomposed each triple `(f1, f2, f3)` into TWO edges: `f1 → f2` and `f2 → f3`. The probability was defined as:

```
P(target | source) = |genomes(src→tgt)| / |genomes where src appears as a source|
```

### Why it broke

For F = `DAHP synthase`, three different "successors" all had ~0.998 probability:

```
('putative secreted protein', 0.9988, 14223)
('Chorismate mutase I / Cyclohexadienyl dehydrogenase', 0.9981, 14214)
('Phosphoenolpyruvate synthase regulatory protein', 0.9908, 14109)
```

These are **not mutually exclusive successors**. The same ~14,200 genomes contain F adjacent to all three because they're all in the same operon (aromatic amino acid biosynthesis cluster). The graph correctly reported co-occurrence frequencies, but they were misinterpreted as Markov transition probabilities.

The decomposition let a single genome contribute to multiple outgoing edges from F simultaneously — once for each triple F participated in. **The outgoing weights are independent conditional co-occurrence frequencies, not a probability distribution.**

### The fix: directional triple decomposition

Each triple `(upstream, center, downstream)` contributes EXACTLY ONE forward edge: `center → downstream`. Now each genome supporting F as a center contributes to exactly one outgoing F→? edge.

```python
def build_transition_graph(fc_adj_pairs: dict, total_genomes: int) -> nx.DiGraph:
    pair_genomes = defaultdict(set)
    center_genomes = defaultdict(set)

    for (upstream, center, downstream), genomes in fc_adj_pairs.items():
        gset = set(genomes)
        center_genomes[center].update(gset)
        pair_genomes[(center, downstream)].update(gset)

    G = nx.DiGraph()
    G.graph['total_genomes'] = total_genomes

    for (src, tgt), gset in pair_genomes.items():
        support = len(gset)
        denom = len(center_genomes[src])
        G.add_edge(src, tgt,
                   weight=support / denom if denom else 0.0,
                   support=support,
                   prevalence=support / total_genomes)
    return G
```

Outgoing weights now sum to ~1 for any node that appears as a center. This is a proper Markov kernel.

## Two Estimator Approaches

### Approach 1: Radius-1 (direct X — F — Y)

Use triples directly. Estimate `P(X | F)` and `P(Y | F)` as marginals. Optionally factor `P(X, Y | F) ≈ P(X | F) · P(Y | F)`.

**Pros:** Maximum data coverage, simple, interpretable.
**Cons:** Independence assumption is wrong (X and Y are correlated through the operon). Cannot disambiguate F's functional context when F is promiscuous.

### Approach 2: Radius-2 (F1 — X — F — Y — F2)

Use quintuples. Condition on flanking pair `(F1, F2)` for the strongest signal.

**Pros:** Captures operon-level context. Mutual constraint from flanking anchors disambiguates which functional module F is in.
**Cons:** Data sparsity. The 5-tuple space is combinatorially larger; same 14k genomes split across many distinct windows. Risk of confident estimates from few observations.

### Biological grounding

Approach 2 makes biological sense because:
- Operons typically span 2–7 genes; a 5-gene window captures most operon contexts.
- Two flanking genes provide much stronger evidence of conserved gene-cluster membership than one.
- Signal decays past ~5–7 genes; coincidental adjacency noise dominates beyond.

### Recommendation: Hierarchical estimator with backoff

Same logic as Kneser-Ney smoothing in n-gram language models. Gene functions are "tokens," contigs are "sentences."

```
P(X, Y | F) = {
  empirical 5-window cond.   if N(F1, _, F, _, F2) ≥ threshold_r2
  empirical 5-window marg.   else if N(_, _, F, _, _) ≥ threshold_r2
  empirical 3-window         else if N(_, F, _) ≥ threshold_r1
  unseen                     otherwise
}
```

## The Full Estimator

```python
import networkx as nx
from collections import defaultdict
from typing import Iterable, Optional


class NeighborhoodEstimator:
    """
    Hierarchical estimator for P(X, Y | F) with backoff:
      Tier 1: r2_conditional — P(X, Y | F, F1, F2)
      Tier 2: r2_marginal    — P(X, Y | F) from quintuples
      Tier 3: r1_independent — P(X | F) · P(Y | F) from triples
      Tier 4: unseen
    """

    def __init__(self, total_genomes: int,
                 min_support_r2: int = 20,
                 min_support_r1: int = 5):
        self.total_genomes = total_genomes
        self.min_support_r2 = min_support_r2
        self.min_support_r1 = min_support_r1

        self._r1_joint = defaultdict(lambda: defaultdict(set))
        self._r1_up = defaultdict(lambda: defaultdict(set))
        self._r1_down = defaultdict(lambda: defaultdict(set))
        self._r1_center_genomes = defaultdict(set)

        self._r2_joint = defaultdict(lambda: defaultdict(set))
        self._r2_context_genomes = defaultdict(set)
        self._r2_marginal = defaultdict(lambda: defaultdict(set))
        self._r2_center_genomes = defaultdict(set)

    def add_triples(self, fc_adj_pairs):
        for (up, center, down), genomes in fc_adj_pairs.items():
            gset = set(genomes)
            self._r1_joint[center][(up, down)].update(gset)
            self._r1_up[center][up].update(gset)
            self._r1_down[center][down].update(gset)
            self._r1_center_genomes[center].update(gset)

    def add_quintuples(self, fc_adj_5tuples):
        for (f1, x, center, y, f2), genomes in fc_adj_5tuples.items():
            gset = set(genomes)
            self._r2_joint[(center, f1, f2)][(x, y)].update(gset)
            self._r2_context_genomes[(center, f1, f2)].update(gset)
            self._r2_marginal[center][(x, y)].update(gset)
            self._r2_center_genomes[center].update(gset)

    def joint(self, center, f1=None, f2=None, top_k=None):
        # Tier 1
        if f1 is not None and f2 is not None:
            ctx = (center, f1, f2)
            denom = len(self._r2_context_genomes.get(ctx, ()))
            if denom >= self.min_support_r2:
                rows = [(xy, len(g) / denom, len(g), 'r2_conditional')
                        for xy, g in self._r2_joint[ctx].items()]
                rows.sort(key=lambda r: -r[1])
                return rows[:top_k] if top_k else rows

        # Tier 2
        denom = len(self._r2_center_genomes.get(center, ()))
        if denom >= self.min_support_r2:
            rows = [(xy, len(g) / denom, len(g), 'r2_marginal')
                    for xy, g in self._r2_marginal[center].items()]
            rows.sort(key=lambda r: -r[1])
            return rows[:top_k] if top_k else rows

        # Tier 3
        denom = len(self._r1_center_genomes.get(center, ()))
        if denom >= self.min_support_r1:
            rows = []
            for x, g_up in self._r1_up[center].items():
                px = len(g_up) / denom
                for y, g_dn in self._r1_down[center].items():
                    py = len(g_dn) / denom
                    support = len(g_up & g_dn)  # true co-occurrence
                    rows.append(((x, y), px * py, support, 'r1_independent'))
            rows.sort(key=lambda r: -r[1])
            return rows[:top_k] if top_k else rows

        # Tier 4
        return [(('?', '?'), 0.0, 0, 'unseen')]

    def marginal_upstream(self, center, top_k=None):
        denom = len(self._r1_center_genomes.get(center, ()))
        if denom < self.min_support_r1:
            return []
        rows = [(x, len(g) / denom, len(g)) for x, g in self._r1_up[center].items()]
        rows.sort(key=lambda r: -r[1])
        return rows[:top_k] if top_k else rows

    def marginal_downstream(self, center, top_k=None):
        denom = len(self._r1_center_genomes.get(center, ()))
        if denom < self.min_support_r1:
            return []
        rows = [(y, len(g) / denom, len(g)) for y, g in self._r1_down[center].items()]
        rows.sort(key=lambda r: -r[1])
        return rows[:top_k] if top_k else rows

    def support_summary(self, center):
        return {
            'r1_center_genomes': len(self._r1_center_genomes.get(center, ())),
            'r1_distinct_pairs': len(self._r1_joint.get(center, {})),
            'r2_center_genomes': len(self._r2_center_genomes.get(center, ())),
            'r2_distinct_pairs': len(self._r2_marginal.get(center, {})),
            'r2_distinct_contexts': sum(1 for (c, _, _) in self._r2_context_genomes if c == center),
        }
```

## Window Construction from Raw Edges

```python
from collections import defaultdict

def build_windows_from_edges(edges, gap_filter=None, strand_filter=True):
    """
    Reconstruct gene order per contig and emit triples and quintuples in 5'→3' direction.

    gap_filter: max bp gap to consider adjacency tight (operon threshold ~200bp).
    strand_filter: if True, only emit windows where all genes share a strand.
    """
    by_contig = defaultdict(list)
    for e in edges:
        by_contig[e['contig']].append(e)

    triples, quintuples = [], []

    for contig, contig_edges in by_contig.items():
        feats = {}
        for e in contig_edges:
            p = e['properties']
            feats[e['source']['entry']] = (p['source_start'], p['source_end'], p['source_strand'])
            feats[e['target']['entry']] = (p['target_start'], p['target_end'], p['target_strand'])

        ordered = sorted(feats.items(), key=lambda kv: kv[1][0])

        # Split into strand-homogeneous runs if requested
        runs = []
        if strand_filter:
            current, current_strand = [], None
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

        # Apply gap filter
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

        for run in runs:
            strand = run[0][3]
            ids = [r[0] for r in run]
            if strand == '-':
                ids = ids[::-1]  # flip to 5'→3'

            for i in range(len(ids) - 2):
                triples.append(tuple(ids[i:i + 3]))
            for i in range(len(ids) - 4):
                quintuples.append(tuple(ids[i:i + 5]))

    return triples, quintuples
```

After this, map feature IDs to functional annotations and aggregate per-window genome lists to build `fc_adj_pairs` and `fc_adj_5tuples`.

## Usage Pattern

```python
est = NeighborhoodEstimator(total_genomes=27254, min_support_r2=20, min_support_r1=5)
est.add_triples(fc_adj_pairs)
est.add_quintuples(fc_adj_5tuples)

F = '2-keto-3-deoxy-D-arabino-heptulosonate-7-phosphate synthase I alpha (EC 2.5.1.54)'

# Marginal joint, auto-tier
for (x, y), p, n, tier in est.joint(F, top_k=10):
    print(f"[{tier}] p={p:.4f} n={n} X={x[:40]} Y={y[:40]}")

# Conditioned on flanks
est.joint(F,
          f1='Phosphoenolpyruvate synthase regulatory protein',
          f2='putative secreted protein',
          top_k=10)

print(est.support_summary(F))
```

## Design Decisions Worth Remembering

**Probability definition for tier 3.** The `r1_independent` tier reports `P(X|F) · P(Y|F)` as the probability but the true intersection `|genomes(X,_) ∩ genomes(_,Y)|` as the support. If actual support is much less than `P(X|F) · P(Y|F) · denom`, X and Y are anti-correlated given F — a useful diagnostic.

**Strand handling.** Lives in window construction, not in the estimator. Once windows are emitted in transcriptional order, the estimator is direction-agnostic.

**Gap filter.** Optional but recommended for operon-level prediction (~200bp threshold). Trades ~30–50% of data for cleaner signal. Leave off for general syntenic context.

**Support thresholds.** Default `min_support_r2=20`, `min_support_r1=5`. Tune empirically using held-out log-likelihood, not by gut.

**Co-occurrence ≠ Markov.** The original mistake: decomposing triples into two edges and treating outgoing weights as a distribution. Outgoing weights only sum to 1 if each genome contributes to exactly one outgoing edge per source. Use directional decomposition (one edge per triple) for proper kernels.

## Open Items / Possible Extensions

- **Held-out log-likelihood evaluator** for empirical threshold tuning. Hold out 10% of genomes, fit on the rest, report mean log P(observed_X, observed_Y | F) per held-out window.
- **Smoothing.** Add-α or full Kneser-Ney style backoff weights rather than the hard tier switch.
- **Asymmetric radii.** Use radius-2 upstream + radius-1 downstream when one side is sparse.
- **Strand-aware sign in the original graph.** The provided `ADJACENT_TO` edges encode coordinate adjacency, not transcriptional. Re-running the pipeline with strand-normalized triples may shift estimates noticeably.
- **Confidence intervals on probabilities.** Wilson or Jeffreys interval on `support / denom` to flag predictions backed by few genomes.

## Biological Context for the DAHP Synthase Example

The DAHP synthase / chorismate mutase / PEP synthase regulator cluster reflects the conserved aromatic amino acid biosynthesis locus (aroF–tyrA–ppsR region in many bacteria). The "putative secreted protein" co-occurring with this cluster across ~14k genomes is likely a conserved hypothetical embedded in the locus; a candidate for functional reannotation.
