"""
feature_chain.py
================

Build position-based edges between genomic features.

Two edge types are produced, both strand-agnostic (strand is carried as a
property rather than used to partition the chain):

  ADJACENT_TO  -- a "nearest-downstream" DAG. Feature A points to every
                  feature B that starts strictly after A ends and whose
                  start is minimal among all such candidates. Ties on the
                  minimal start all receive an edge, which is how the
                  non-linear  A->B, A->C  fan-out arises. Re-convergence
                  (B->D, C->D) falls out naturally when several features
                  share the same nearest successor.

  OVERLAPS     -- an undirected edge between any two features whose
                  [start, end] intervals intersect. Emitted once per pair.

Design notes
------------
* Coordinates are treated as a closed interval [start, end]. Two intervals
  overlap iff  A.start <= B.end  and  B.start <= A.end.
* An overlapping feature is deliberately NOT an adjacency successor, so the
  two edge sets are disjoint in meaning: ADJACENT_TO is "what comes next
  with a clean gap", OVERLAPS is everything that intersects.
* Edges are keyed by the feature `entry` field. `elementId` is also carried
  on each endpoint in case you need the Neo4j-native id downstream.
* Features should belong to the same contig/sequence before being passed in
  here -- this function does not know about contig membership. Group by
  contig upstream and call once per group (see `chain_by_contig`).
"""

from __future__ import annotations

from itertools import groupby
from typing import Any, Callable, Dict, Iterable, List


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def _node_ref(feature: Dict[str, Any]) -> Dict[str, Any]:
    """Compact endpoint reference carried on every edge."""
    return {
        "entry": feature["entry"],
        "elementId": feature.get("elementId"),
    }


def _coords(feature: Dict[str, Any]) -> tuple[int, int]:
    """Return (start, end) from a feature's properties, normalized so
    start <= end regardless of how the source recorded them."""
    props = feature["properties"]
    s, e = int(props["start"]), int(props["end"])
    return (s, e) if s <= e else (e, s)


def _strand(feature: Dict[str, Any]) -> Any:
    return feature["properties"].get("strand")


# --------------------------------------------------------------------------
# adjacency
# --------------------------------------------------------------------------

def adjacency_edges(features: Iterable[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """ADJACENT_TO edges: each feature -> its nearest non-overlapping
    downstream feature(s).

    For feature A, candidates are all features B with B.start > A.end.
    Among those, the minimal B.start is selected; every B at that minimal
    start gets an edge. This yields a DAG, not a linked list: a feature can
    fan out to several successors (A->B, A->C) and several features can
    converge on one (B->D, C->D).
    """
    feats = list(features)

    # Pre-compute coords once; sort by (start, end) for deterministic output.
    indexed = sorted(
        ((f, *_coords(f)) for f in feats),
        key=lambda t: (t[1], t[2]),
    )

    edges: List[Dict[str, Any]] = []

    for src, _src_start, src_end in indexed:
        # Candidates start strictly after this feature ends.
        candidates = [
            (f, s, e) for (f, s, e) in indexed if s > src_end
        ]
        if not candidates:
            continue

        min_start = min(s for (_f, s, _e) in candidates)
        successors = [t for t in candidates if t[1] == min_start]

        for dst, dst_start, dst_end in successors:
            src_start_c, _ = _coords(src)
            edges.append({
                "type": "ADJACENT_TO",
                "source": _node_ref(src),
                "target": _node_ref(dst),
                "properties": {
                    "gap": dst_start - src_end - 1,   # bases between them
                    "source_strand": _strand(src),
                    "target_strand": _strand(dst),
                    "source_start": src_start_c,
                    "source_end": src_end,
                    "target_start": dst_start,
                    "target_end": dst_end,
                },
            })

    return edges


# --------------------------------------------------------------------------
# overlap
# --------------------------------------------------------------------------

def overlap_edges(features: Iterable[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """OVERLAPS edges: one undirected edge per pair of intersecting
    intervals. Uses a sweep so it is O(n log n + k) rather than O(n^2),
    where k is the number of overlapping pairs.
    """
    indexed = sorted(
        ((f, *_coords(f)) for f in features),
        key=lambda t: (t[1], t[2]),
    )

    edges: List[Dict[str, Any]] = []
    active: List[tuple] = []  # features whose interval may still overlap

    for cur, cur_start, cur_end in indexed:
        # Drop anything that ends before the current feature starts.
        active = [t for t in active if t[2] >= cur_start]
        # Everything still active overlaps `cur` (closed-interval test).
        for other, other_start, other_end in active:
            ov_start = max(cur_start, other_start)
            ov_end = min(cur_end, other_end)
            edges.append({
                "type": "OVERLAPS",
                "source": _node_ref(other),   # ordered: earlier start first
                "target": _node_ref(cur),
                "properties": {
                    "overlap_start": ov_start,
                    "overlap_end": ov_end,
                    "overlap_length": ov_end - ov_start + 1,
                    "source_strand": _strand(other),
                    "target_strand": _strand(cur),
                    "same_strand": _strand(other) == _strand(cur),
                },
            })
        active.append((cur, cur_start, cur_end))

    return edges


# --------------------------------------------------------------------------
# combined entry points
# --------------------------------------------------------------------------

def feature_edges(features: Iterable[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Both edge types for a single set of (same-contig) features."""
    feats = list(features)
    return adjacency_edges(feats) + overlap_edges(feats)


def chain_by_contig(
    features: Iterable[Dict[str, Any]],
    contig_of: Callable[[Dict[str, Any]], Any],
) -> List[Dict[str, Any]]:
    """Convenience wrapper: group features by contig, then build edges
    within each group only -- so adjacency/overlap never cross a contig
    boundary.

    NOTE on the supplied GenomicFeature record: it carries
    `entry`, `type`, `labels`, `properties` (feature_type/strand/start/
    end/source), `elementId`, `createdAt`, `updatedAt` -- and NO contig
    field. `source` is 'KBERDL', which is a provenance tag, not a contig.

    So contig membership must come from outside this record. Pass a
    `contig_of` callable that resolves it however your graph models it,
    e.g.:

        # if contig membership is a (LocusTag/Contig)-[:HAS_FEATURE]-> edge
        # resolved upstream into a dict:
        chain_by_contig(feats, lambda f: feature_to_contig[f["entry"]])

        # or if you later add it to properties:
        chain_by_contig(feats, lambda f: f["properties"]["contig"])

    Features that resolve to the same key are chained together; different
    keys are kept separate.
    """
    feats = list(features)
    feats.sort(key=lambda f: (str(contig_of(f)), *_coords(f)))

    all_edges: List[Dict[str, Any]] = []
    for _contig, group in groupby(feats, key=contig_of):
        all_edges.extend(feature_edges(list(group)))
    return all_edges
