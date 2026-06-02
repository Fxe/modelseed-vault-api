"""
conformity.py — Genome conformity scoring against a NeighborhoodEstimator.

Given a fitted NeighborhoodEstimator (built from a reference corpus), score
every window in a new genome and produce:

  - per-window log-probability and z-score
  - hotspot regions of consecutive anomalous windows
  - genome-level summary statistics

Null model
----------
fit_null() scores all training windows, collecting log P values per center
function.  Only seen windows (P > 0) contribute to the null; unseen windows
are always flagged as anomalous regardless of z-score threshold.

Z-score convention: higher = more conformant.
  z >> 0  →  neighborhood is more predictable than average (conserved)
  z ~  0  →  normal
  z << 0  →  anomalous (annotation error candidate, HGT, or OOD clade)
"""
from __future__ import annotations

import json
import math
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from .genome_net import NeighborhoodEstimator, build_windows_from_edges

# Log-probability assigned to windows with P = 0.
FLOOR_LOG_PROB: float = -20.0


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class WindowScore:
    """Score for a single genomic window (triple or quintuple)."""
    window: tuple           # function-name tuple, length 3 or 5
    center: str             # center function (index 1 for triples, 2 for quintuples)
    prob: float             # P from estimator  (0.0 if unseen)
    log_prob: float         # log P; equals FLOOR_LOG_PROB when prob == 0
    tier: str               # 'r2_conditional' | 'r2_marginal' | 'r1_joint' |
                            # 'r1_independent' | 'unseen'
    zscore: float           # (log_prob − null_mean) / null_std; NaN if no null

    @property
    def surprise(self) -> float:
        """Positive score: higher = more surprising."""
        return -self.log_prob

    @property
    def is_unseen(self) -> bool:
        return self.tier == 'unseen'


@dataclass
class Hotspot:
    """A run of consecutive anomalous windows."""
    window_start: int       # index into GenomeReport.window_scores (inclusive)
    window_end: int         # inclusive
    length: int             # number of windows in the run
    mean_zscore: float
    min_zscore: float       # most anomalous point in the run
    centers: list[str]      # center functions of each window in the run


@dataclass
class GenomeReport:
    """Full conformity report for one genome."""
    genome_id: str
    window_scores: list[WindowScore]
    hotspots: list[Hotspot]

    # ---- summary statistics ----
    global_log_prob: float      # mean log P over all windows
    global_zscore: float        # mean z-score (NaN if null not fitted)
    fraction_unseen: float      # fraction of windows with tier == 'unseen'
    fraction_r2: float          # fraction of windows resolved at r2 tier
    n_windows: int

    # ---- convenience accessors ----
    def worst_windows(self, n: int = 10) -> list[WindowScore]:
        """Most anomalous windows (lowest z-score first)."""
        return sorted(
            self.window_scores,
            key=lambda s: s.zscore if not math.isnan(s.zscore) else -1e9
        )[:n]

    def best_windows(self, n: int = 10) -> list[WindowScore]:
        """Most conformant windows (highest z-score first)."""
        return sorted(
            self.window_scores,
            key=lambda s: s.zscore if not math.isnan(s.zscore) else -1e9,
            reverse=True
        )[:n]

    def tier_counts(self) -> dict[str, int]:
        counts: dict[str, int] = defaultdict(int)
        for ws in self.window_scores:
            counts[ws.tier] += 1
        return dict(counts)


# ---------------------------------------------------------------------------
# Utility: extract windows from raw genome files
# ---------------------------------------------------------------------------

def extract_windows(
    file_ec: Path,
    file_chain: Path,
    gap_filter: Optional[int] = None,
    strand_filter: bool = True,
) -> tuple[list[tuple], list[tuple]]:
    """
    Load one genome's ec + chain files and return function-name windows.

    Parameters
    ----------
    file_ec : Path
        JSON mapping feature_id → [genome_hash, function_name, ec_numbers].
    file_chain : Path
        JSON list of graph edges (ADJACENT_TO, OVERLAPS, …).
    gap_filter : int or None
        Maximum intergenic gap (bp) allowed inside a window.
    strand_filter : bool
        If True, emit only same-strand windows.

    Returns
    -------
    (triples, quintuples)
        Each is a list of function-name tuples (length 3 or 5).
        Windows containing any unannotated gene (None function) are dropped.
    """
    triples_wf, quintuples_wf = extract_windows_with_features(
        file_ec, file_chain, gap_filter=gap_filter, strand_filter=strand_filter
    )
    return [w for w, _ in triples_wf], [w for w, _ in quintuples_wf]  # type: ignore[return-value]


def extract_windows_with_features(
    file_ec: Path,
    file_chain: Path,
    gap_filter: Optional[int] = None,
    strand_filter: bool = True,
) -> tuple[list[tuple[tuple, str]], list[tuple[tuple, str]]]:
    """
    Like extract_windows but also returns the center feature ID for each window.

    Returns
    -------
    (triples, quintuples)
        Each is a list of (fn_window_tuple, center_feature_id) pairs.
        The center feature ID is the middle element of the raw feature-ID window
        (index 1 for triples, index 2 for quintuples).
        Windows containing any unannotated gene are dropped.
    """
    with open(file_ec) as fh:
        feature_ec = json.load(fh)
    with open(file_chain) as fh:
        feature_chain = json.load(fh)

    edges = [o for o in feature_chain if o.get('type') == 'ADJACENT_TO']
    raw_triples, raw_quintuples = build_windows_from_edges(
        edges, gap_filter=gap_filter, strand_filter=strand_filter
    )

    def _fn(fid: str) -> Optional[str]:
        v = feature_ec.get(fid)
        return v[1] if v else None

    triples: list[tuple[tuple, str]] = []
    for t in raw_triples:
        if all(_fn(fid) is not None for fid in t):
            triples.append((tuple(_fn(fid) for fid in t), t[1]))  # type: ignore[arg-type]

    quintuples: list[tuple[tuple, str]] = []
    for q in raw_quintuples:
        if all(_fn(fid) is not None for fid in q):
            quintuples.append((tuple(_fn(fid) for fid in q), q[2]))  # type: ignore[arg-type]

    return triples, quintuples


# ---------------------------------------------------------------------------
# ConformityScorer
# ---------------------------------------------------------------------------

class ConformityScorer:
    """
    Score genomic windows against a fitted NeighborhoodEstimator.

    Workflow
    --------
    >>> scorer = ConformityScorer(estimator)
    >>> scorer.fit_null(builder.fc_adj_pairs, builder.fc_adj_5tuples)
    >>> report  = scorer.score_genome('genome_A', triples + quintuples)

    Parameters
    ----------
    estimator : NeighborhoodEstimator
        Pre-fitted estimator.
    floor_log_prob : float
        Log-probability assigned to unseen windows (P = 0).
    hotspot_zscore_threshold : float
        Windows with z-score < this value are considered anomalous.
    hotspot_min_length : int
        Minimum number of consecutive anomalous windows to report as hotspot.
    """

    def __init__(
        self,
        estimator: NeighborhoodEstimator,
        floor_log_prob: float = FLOOR_LOG_PROB,
        hotspot_zscore_threshold: float = -1.5,
        hotspot_min_length: int = 2,
    ):
        self.est = estimator
        self.floor = floor_log_prob
        self.hotspot_threshold = hotspot_zscore_threshold
        self.hotspot_min_length = hotspot_min_length

        # (mean_log_prob, std_log_prob) — built by fit_null()
        self._global_null: Optional[tuple[float, float]] = None
        self._per_center_null: dict[str, tuple[float, float]] = {}

    # ------------------------------------------------------------------ lookup

    def lookup_prob(self, window: tuple) -> tuple[float, str]:
        """
        Return (probability, tier) for a specific window without ranking.

        Directly reads the estimator's internal data structures for O(1) lookup.
        Returns (0.0, 'unseen') when the window has no model support.
        """
        n = len(window)
        if n == 3:
            return self._lookup_triple(window[0], window[1], window[2])
        if n == 5:
            return self._lookup_quintuple(*window)
        raise ValueError(f'Window length must be 3 or 5, got {n}')

    def _lookup_triple(
        self, up: str, center: str, down: str
    ) -> tuple[float, str]:
        est = self.est
        center_g = est._r1_center_genomes.get(center)
        if not center_g or len(center_g) < est.min_support_r1:
            return 0.0, 'unseen'
        denom = len(center_g)

        # Exact co-occurrence
        joint_g = est._r1_joint.get(center, {}).get((up, down))
        if joint_g:
            return len(joint_g) / denom, 'r1_joint'

        # Independent factorization
        up_g = est._r1_up.get(center, {}).get(up)
        dn_g = est._r1_down.get(center, {}).get(down)
        if up_g and dn_g:
            return (len(up_g) / denom) * (len(dn_g) / denom), 'r1_independent'

        return 0.0, 'unseen'

    def _lookup_quintuple(
        self, f1: str, x: str, center: str, y: str, f2: str
    ) -> tuple[float, str]:
        est = self.est

        # Tier 1 — r2_conditional
        ctx = (center, f1, f2)
        ctx_g = est._r2_context_genomes.get(ctx)
        if ctx_g and len(ctx_g) >= est.min_support_r2:
            pair_g = est._r2_joint.get(ctx, {}).get((x, y))
            if pair_g:
                return len(pair_g) / len(ctx_g), 'r2_conditional'

        # Tier 2 — r2_marginal
        r2_g = est._r2_center_genomes.get(center)
        if r2_g and len(r2_g) >= est.min_support_r2:
            pair_g = est._r2_marginal.get(center, {}).get((x, y))
            if pair_g:
                return len(pair_g) / len(r2_g), 'r2_marginal'

        # Fall back to triple (center window)
        return self._lookup_triple(x, center, y)

    # ------------------------------------------------------------------ null

    def fit_null(
        self,
        fc_adj_pairs: Optional[dict] = None,
        fc_adj_5tuples: Optional[dict] = None,
    ) -> ConformityScorer:
        """
        Build per-center and global log-P null distributions from training windows.

        Only windows where P > 0 contribute (unseen windows are excluded so
        z-scores reflect deviations from typical *seen* windows).

        Parameters
        ----------
        fc_adj_pairs : dict, optional
            Training triples dict from QuintuplePairBuilder.fc_adj_pairs.
        fc_adj_5tuples : dict, optional
            Training quintuples dict from QuintuplePairBuilder.fc_adj_5tuples.
        """
        per_center: dict[str, list[float]] = defaultdict(list)
        all_lp: list[float] = []

        def _collect(window: tuple) -> None:
            p, _ = self.lookup_prob(window)
            if p <= 0:
                return
            lp = max(math.log(p), self.floor)
            center = window[len(window) // 2]
            per_center[center].append(lp)
            all_lp.append(lp)

        # Prefer quintuples (more informative tier coverage)
        if fc_adj_5tuples:
            for key in fc_adj_5tuples:
                _collect(key)
        if fc_adj_pairs:
            for key in fc_adj_pairs:
                _collect(key)

        if not all_lp:
            return self

        # Global null
        mean_g = sum(all_lp) / len(all_lp)
        std_g = math.sqrt(
            sum((x - mean_g) ** 2 for x in all_lp) / len(all_lp)
        ) + 1e-9
        self._global_null = (mean_g, std_g)

        # Per-center null (require ≥ 5 samples for stability)
        for center, lps in per_center.items():
            if len(lps) < 5:
                continue
            m = sum(lps) / len(lps)
            s = math.sqrt(sum((x - m) ** 2 for x in lps) / len(lps)) + 1e-9
            self._per_center_null[center] = (m, s)

        return self

    def null_stats(self) -> dict:
        """Summary of the fitted null distributions."""
        return {
            'global_null_mean': self._global_null[0] if self._global_null else None,
            'global_null_std': self._global_null[1] if self._global_null else None,
            'n_centers_with_null': len(self._per_center_null),
        }

    # ------------------------------------------------------------------ scoring

    def _zscore(self, log_prob: float, center: str) -> float:
        """Z-score relative to per-center null (falls back to global null).

        Per-center std is floored at 10% of the global std to prevent
        blow-up when a center has very consistent (near-identical) training log-probs.
        """
        null = self._per_center_null.get(center) or self._global_null
        if null is None:
            return float('nan')
        mean, std = null
        # Floor std at a fraction of the global std to avoid division by near-zero
        if self._global_null is not None:
            std = max(std, 0.1 * self._global_null[1])
        return (log_prob - mean) / std

    def score_window(self, window: tuple) -> WindowScore:
        """Score a single window tuple."""
        prob, tier = self.lookup_prob(window)
        lp = max(math.log(prob), self.floor) if prob > 0 else self.floor
        center = window[len(window) // 2]
        return WindowScore(
            window=window,
            center=center,
            prob=prob,
            log_prob=lp,
            tier=tier,
            zscore=self._zscore(lp, center),
        )

    def score_genome(
        self,
        genome_id: str,
        windows: list[tuple],
    ) -> GenomeReport:
        """
        Score all windows for one genome.

        Parameters
        ----------
        genome_id : str
            Identifier for this genome (used in the report only).
        windows : list of tuples
            Ordered list of function-name windows (triples and/or quintuples).
            Produce these with extract_windows().

        Returns
        -------
        GenomeReport
        """
        scored = [self.score_window(w) for w in windows]

        if not scored:
            return GenomeReport(
                genome_id=genome_id,
                window_scores=[],
                hotspots=[],
                global_log_prob=float('nan'),
                global_zscore=float('nan'),
                fraction_unseen=float('nan'),
                fraction_r2=float('nan'),
                n_windows=0,
            )

        n = len(scored)
        log_probs = [s.log_prob for s in scored]
        valid_z = [s.zscore for s in scored if not math.isnan(s.zscore)]

        return GenomeReport(
            genome_id=genome_id,
            window_scores=scored,
            hotspots=self._detect_hotspots(scored),
            global_log_prob=sum(log_probs) / n,
            global_zscore=sum(valid_z) / len(valid_z) if valid_z else float('nan'),
            fraction_unseen=sum(1 for s in scored if s.is_unseen) / n,
            fraction_r2=sum(1 for s in scored if s.tier.startswith('r2')) / n,
            n_windows=n,
        )

    # ------------------------------------------------------------------ hotspots

    def _detect_hotspots(self, scored: list[WindowScore]) -> list[Hotspot]:
        """
        Identify runs of consecutive anomalous windows.

        A window is anomalous if:
          - its z-score < hotspot_zscore_threshold, OR
          - it is unseen (P = 0, always anomalous regardless of threshold)
        """
        hotspots: list[Hotspot] = []
        run_start: Optional[int] = None

        def _is_anomalous(ws: WindowScore) -> bool:
            if ws.is_unseen:
                return True
            return not math.isnan(ws.zscore) and ws.zscore < self.hotspot_threshold

        def _close_run(end_exclusive: int) -> None:
            length = end_exclusive - run_start  # type: ignore[operator]
            if length < self.hotspot_min_length:
                return
            run = scored[run_start:end_exclusive]  # type: ignore[index]
            valid_z = [s.zscore for s in run if not math.isnan(s.zscore)]
            hotspots.append(Hotspot(
                window_start=run_start,  # type: ignore[arg-type]
                window_end=end_exclusive - 1,
                length=length,
                mean_zscore=sum(valid_z) / len(valid_z) if valid_z else float('nan'),
                min_zscore=min(valid_z, default=float('nan')),
                centers=[s.center for s in run],
            ))

        for i, ws in enumerate(scored):
            if _is_anomalous(ws):
                if run_start is None:
                    run_start = i
            else:
                if run_start is not None:
                    _close_run(i)
                    run_start = None

        if run_start is not None:
            _close_run(len(scored))

        return hotspots

    # ------------------------------------------------------------------ batch

    def score_corpus(
        self,
        genome_windows: dict[str, list[tuple]],
    ) -> dict[str, GenomeReport]:
        """
        Score multiple genomes at once.

        Parameters
        ----------
        genome_windows : dict mapping genome_id → list of window tuples.

        Returns
        -------
        dict mapping genome_id → GenomeReport.
        """
        return {gid: self.score_genome(gid, wins) for gid, wins in genome_windows.items()}

    def corpus_summary(
        self,
        reports: dict[str, GenomeReport],
    ) -> list[dict]:
        """
        Tabular summary across multiple GenomeReports.

        Returns a list of dicts (one per genome), suitable for DataFrame construction.
        """
        rows = []
        for gid, r in reports.items():
            rows.append({
                'genome_id': gid,
                'n_windows': r.n_windows,
                'global_log_prob': r.global_log_prob,
                'global_zscore': r.global_zscore,
                'fraction_unseen': r.fraction_unseen,
                'fraction_r2': r.fraction_r2,
                'n_hotspots': len(r.hotspots),
                'hotspot_windows': sum(h.length for h in r.hotspots),
            })
        rows.sort(key=lambda r: r['global_zscore']
                  if not math.isnan(r['global_zscore']) else -1e9)
        return rows
