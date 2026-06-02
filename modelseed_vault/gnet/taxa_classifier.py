"""
taxa_classifier.py — Hierarchical conformity-based taxonomic labelling of gene
neighborhoods.

Given a hierarchy of fitted ``ConformityScorer`` objects (one per taxonomic rank,
ordered broad → narrow, e.g. family → genus → species), this module labels each
window (triple or quintuple) of a query genome by the **taxonomic scope at which
its gene neighborhood is statistically conserved** — producing a per-gene
"synteny-conservation-depth" track.

What this measures (and what it does NOT)
-----------------------------------------
The estimators are *nested marginalisations of the same data* (species genomes ⊂
genus ⊂ family), so this is **not** model selection between independent hypotheses
and it is **not** a readout of evolutionary origin / phylogenetic ownership. It
answers a narrower, defensible question:

    "At what taxonomic scope does conditioning on a finer clade label stop
     improving the prediction of this window's neighborhood?"

A window that the *family* model already explains is part of the conserved
backbone (CORE); one that only a *genus*/*species* model explains is
clade-characteristic/lineage-specific; one no competent model explains is
UNEXPLAINED (HGT / novel biology / assembly or annotation artifact).

Soundness controls implemented here
-----------------------------------
1. **Support-aware abstention.** A level may only judge a window if it is
   *competent* for the window's center function — defined as having a calibrated
   per-center null (``center in scorer._per_center_null``, i.e. ≥5 seen training
   windows for that center). A level with too little data **abstains** rather than
   voting "not conserved". This prevents conflating *under-sampled* with
   *not-conserved* — the central confound when species corpora are small.
2. **Calibrated cross-level comparison.** Decisions use the per-center z-score
   (each level normalised by its *own* null), never raw log-probabilities across
   differently-sized corpora.
3. **Explicit evidence.** Every label carries the per-level (z, prob, tier,
   support, competent, conformant) so the caller can audit *why* a region was
   labelled and never has to trust a bare class.

Leave-one-out (avoiding leakage) is the caller's responsibility: build the
hierarchy with the query genome (and, for a strict test, its species) held out.

See ``gnet_taxa_classifier.ipynb`` for an end-to-end demo with negative controls.
"""
from __future__ import annotations

import json
import math
import random
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

from .conformity import ConformityScorer

# Label used for the broadest level in the hierarchy when it is the deepest
# (i.e. broadest) level that already explains a window.
CORE_LABEL = "CORE"
# Window explained by no competent level (a competent level rejected it).
UNEXPLAINED = "UNEXPLAINED"
# No level had enough data to judge the window's center function.
UNJUDGED = "UNJUDGED"


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class LevelScore:
    """One taxonomic level's verdict on a single window."""
    rank: str
    zscore: float
    prob: float
    tier: str
    support: int        # # training genomes backing this window's center at this level
    competent: bool     # level has a calibrated null for the center → may judge
    conformant: bool    # competent AND zscore >= conformity_threshold


@dataclass
class WindowTaxaLabel:
    """Taxonomic-scope label for one window."""
    window: tuple
    center: str
    levels: list[LevelScore]            # broad → narrow
    label: str                          # CORE | <rank> | UNEXPLAINED | UNJUDGED
    assigned_rank: Optional[str]        # broadest conformant rank, or None
    gains: dict[str, float]             # rank → z gain vs previous competent level

    @property
    def is_explained(self) -> bool:
        return self.assigned_rank is not None


@dataclass
class GenomeTaxaReport:
    """Per-genome conservation-depth report."""
    genome_id: str
    window_labels: list[WindowTaxaLabel]
    rank_order: list[str]

    @property
    def n_windows(self) -> int:
        return len(self.window_labels)

    def label_counts(self) -> Counter:
        return Counter(w.label for w in self.window_labels)

    def fraction(self, label: str) -> float:
        n = self.n_windows
        return self.label_counts().get(label, 0) / n if n else float("nan")

    def track(self) -> list[str]:
        """Per-window label in genomic (within-contig sliding) order."""
        return [w.label for w in self.window_labels]

    def windows_with_label(self, label: str) -> list[WindowTaxaLabel]:
        return [w for w in self.window_labels if w.label == label]

    def top_by_gain(self, rank: str, n: int = 10) -> list[WindowTaxaLabel]:
        """Windows most strongly attributed to `rank` (largest conformity gain)."""
        cand = [w for w in self.window_labels if w.assigned_rank == rank]
        cand.sort(key=lambda w: w.gains.get(rank, float("-inf")), reverse=True)
        return cand[:n]

    def summary(self) -> dict:
        c = self.label_counts()
        n = max(1, self.n_windows)
        out = {"genome_id": self.genome_id, "n_windows": self.n_windows}
        for lab in [CORE_LABEL] + [r for r in self.rank_order[1:]] + [UNEXPLAINED, UNJUDGED]:
            out[f"frac_{lab}"] = c.get(lab, 0) / n
        return out


# ---------------------------------------------------------------------------
# Classifier
# ---------------------------------------------------------------------------

class ConformityTaxaClassifier:
    """
    Label windows by the taxonomic scope at which their neighborhood is conserved.

    Parameters
    ----------
    levels : sequence of (rank_name, ConformityScorer), ordered BROAD → NARROW
        e.g. [('family', fam_scorer), ('genus', genus_scorer), ('species', sp_scorer)].
        Each scorer must have had ``fit_null`` called on its training windows.
    conformity_threshold : float
        A window is "conformant" at a level if its per-center z-score ≥ this value
        (higher z = more conformant). Default -1.0 ≈ "not anomalous".
    parsimony_margin : float
        A finer (narrower) level only "claims" a window if its calibrated z-score
        exceeds the best broader level's by at least this margin. This is an
        Occam guard: do not attribute a region to a narrower clade unless finer
        conditioning *demonstrably improves* the neighborhood's predictability.
        Larger → more conservative (more CORE); ``inf`` recovers a pure
        broadest-conformant backoff; ``0`` recovers a pure best-fit arg-max.
        Default 1.0 (≈ one standard deviation of improvement).
    core_label : str
        Label given to windows best explained at the broadest level.

    Decision rule (support-aware, calibrated, parsimonious)
    -------------------------------------------------------
    For each window, every level reports a calibrated per-center z-score, but only
    *competent* levels (those with a fitted per-center null) may vote. Among the
    competent **and** conformant levels we take the best-fitting z, then assign the
    window to the **broadest** level whose z is within ``parsimony_margin`` of that
    best. So a window is called:
        - 'species'      if species conditioning beats genus/family by the margin
        - 'genus'        if genus conditioning beats family by the margin
        - CORE           if the broadest level already fits within the margin
        - UNEXPLAINED    if some level is competent but none is conformant
        - UNJUDGED       if no level is competent (no statistical power)

    The output is therefore "the *narrowest* taxonomic scope at which the
    neighborhood becomes meaningfully better explained" — i.e. conservation depth,
    not phylogenetic origin.
    """

    def __init__(
        self,
        levels: Sequence[tuple[str, ConformityScorer]],
        conformity_threshold: float = -1.0,
        parsimony_margin: float = 1.0,
        core_label: str = CORE_LABEL,
    ):
        if not levels:
            raise ValueError("levels must be non-empty")
        self.levels = list(levels)
        self.rank_order = [r for r, _ in self.levels]
        self.tau = conformity_threshold
        self.parsimony_margin = parsimony_margin
        self.core_label = core_label

    # ------------------------------------------------------------------ helpers

    @staticmethod
    def _support(scorer: ConformityScorer, center: str) -> int:
        est = scorer.est
        return max(
            len(est._r1_center_genomes.get(center, ())),
            len(est._r2_center_genomes.get(center, ())),
        )

    @staticmethod
    def _competent(scorer: ConformityScorer, center: str) -> bool:
        """A level may judge a center iff it has a calibrated per-center null.

        ``fit_null`` only records a per-center null when ≥5 training windows were
        seen for that center, so this ties competence directly to having enough
        data to calibrate a z-score.
        """
        return center in scorer._per_center_null

    # ------------------------------------------------------------------ core

    def classify_window(self, window: tuple) -> WindowTaxaLabel:
        center = window[len(window) // 2]

        level_scores: list[LevelScore] = []
        for rank, scorer in self.levels:
            ws = scorer.score_window(window)
            comp = self._competent(scorer, center)
            conf = comp and (not math.isnan(ws.zscore)) and ws.zscore >= self.tau
            level_scores.append(LevelScore(
                rank=rank, zscore=ws.zscore, prob=ws.prob, tier=ws.tier,
                support=self._support(scorer, center),
                competent=comp, conformant=conf,
            ))

        # conformity-gain decomposition: z gain over the previous *competent* level
        gains: dict[str, float] = {}
        prev_z: Optional[float] = None
        for ls in level_scores:
            if ls.competent and not math.isnan(ls.zscore):
                gains[ls.rank] = ls.zscore if prev_z is None else (ls.zscore - prev_z)
                prev_z = ls.zscore

        # parsimonious best-fit: best z among competent+conformant levels, then the
        # broadest level within `parsimony_margin` of that best.
        conformant = [ls for ls in level_scores if ls.conformant]
        if conformant:
            best_z = max(ls.zscore for ls in conformant)
            winner = next(ls for ls in conformant
                          if ls.zscore >= best_z - self.parsimony_margin)  # broad → narrow
            assigned_rank = winner.rank
            label = self.core_label if assigned_rank == self.rank_order[0] else assigned_rank
        elif any(ls.competent for ls in level_scores):
            assigned_rank, label = None, UNEXPLAINED
        else:
            assigned_rank, label = None, UNJUDGED

        return WindowTaxaLabel(
            window=window, center=center, levels=level_scores,
            label=label, assigned_rank=assigned_rank, gains=gains,
        )

    def classify_genome(self, genome_id: str, windows: list[tuple]) -> GenomeTaxaReport:
        return GenomeTaxaReport(
            genome_id=genome_id,
            window_labels=[self.classify_window(w) for w in windows],
            rank_order=list(self.rank_order),
        )


# ---------------------------------------------------------------------------
# Shuffle negative control
# ---------------------------------------------------------------------------

def contig_function_lists(file_ec: Path, file_chain: Path) -> list[list[Optional[str]]]:
    """Reconstruct per-contig gene-function order (5'→3' by start coordinate).

    Returns one list of function names (or None for unannotated genes) per contig.
    Used to build a synteny-destroying null: shuffle each contig's functions and
    re-slide windows (matched composition and annotation density, scrambled order).
    """
    with open(file_ec) as fh:
        feature_ec = json.load(fh)
    with open(file_chain) as fh:
        feature_chain = json.load(fh)

    edges = [o for o in feature_chain if o.get("type") == "ADJACENT_TO"]
    by_contig: dict[str, list] = defaultdict(list)
    for e in edges:
        by_contig[e["contig"]].append(e)

    def fn(fid: str) -> Optional[str]:
        v = feature_ec.get(fid)
        return v[1] if v else None

    runs: list[list[Optional[str]]] = []
    for _, ce in by_contig.items():
        feats: dict[str, int] = {}
        for e in ce:
            p = e["properties"]
            feats[e["source"]["entry"]] = p["source_start"]
            feats[e["target"]["entry"]] = p["target_start"]
        ordered = sorted(feats.items(), key=lambda kv: kv[1])
        runs.append([fn(fid) for fid, _ in ordered])
    return runs


def windows_from_runs(runs: list[list[Optional[str]]], radius: int = 2) -> list[tuple]:
    """Slide windows of size ``2*radius+1`` over function runs, dropping any with None."""
    width = 2 * radius + 1
    out: list[tuple] = []
    for run in runs:
        for i in range(len(run) - width + 1):
            w = tuple(run[i:i + width])
            if None not in w:
                out.append(w)
    return out


def shuffled_quintuples(
    file_ec: Path, file_chain: Path, seed: int = 0
) -> list[tuple]:
    """Synteny-destroying null: per-contig shuffle of gene functions, then quintuples."""
    rng = random.Random(seed)
    runs = contig_function_lists(file_ec, file_chain)
    shuffled = [rng.sample(run, len(run)) for run in runs]
    return windows_from_runs(shuffled, radius=2)
