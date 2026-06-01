"""
imputation.py — Missing gene-function imputation from neighborhood context.

Given a genomic window where one position is unknown or uncertain, predict the
most likely function using the inverse of the NeighborhoodEstimator's probability.

Mathematical basis
------------------
The estimator stores P(X, Y | F).  For imputation we want P(F | X, Y).
By Bayes' theorem and the way the estimator's counts are built:

    P(F | X, Y)  ∝  P(X, Y | F) · P(F)
                 =  (count(X, F, Y) / count(F)) · (count(F) / N)
                 =  count(X, F, Y) / N

The prior P(F) cancels in the normalization.  The posterior is therefore simply:

    P(F | X, Y)  =  count(X, F, Y) / Σ_{F'} count(X, F', Y)

i.e. "among all training genomes where the X–?–Y pattern is observed, what
fraction have F as the unknown center?"

The same simplification holds for the one-sided and quintuple-level queries.

Query types
-----------
1. impute_center(up, down [, f1, f2])
        Predict F given immediate neighbors (and optionally outer context).
        Tiers: r2_conditional → r2_marginal → r1_joint → r1_marginal

2. impute_upstream(center [, downstream])
        Predict X given F (and optionally Y).
        Tiers: r1_joint (if Y given) → r1_marginal

3. impute_downstream(center [, upstream])
        Predict Y given F (and optionally X).
        Tiers: r1_joint (if X given) → r1_marginal
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Optional

from .genome_net import NeighborhoodEstimator


# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------

@dataclass
class ImputationResult:
    """
    Result of a single imputation query.

    predictions : list of (function, posterior_prob, support, tier)
        Sorted by posterior_prob descending.
        - posterior_prob : normalised probability in [0, 1]
        - support        : raw genome count backing this candidate
        - tier           : data source ('r2_conditional', 'r2_marginal',
                           'r1_joint', 'r1_marginal', or 'unseen')
    """
    query_type: str            # 'center' | 'upstream' | 'downstream'
    known: dict                # position label → function string
    predictions: list          # [(fn, posterior, support, tier)]
    tier: str                  # tier that resolved the query
    total_support: int         # Σ support across all candidates
    n_candidates: int          # number of distinct candidate functions

    # ---- convenience ----

    @property
    def top1(self) -> Optional[str]:
        return self.predictions[0][0] if self.predictions else None

    @property
    def top1_prob(self) -> float:
        return self.predictions[0][1] if self.predictions else 0.0

    @property
    def is_resolved(self) -> bool:
        return self.tier != 'unseen'

    def rank_of(self, fn: str) -> Optional[int]:
        """1-based rank of fn in predictions; None if fn is not a candidate."""
        for i, (f, *_) in enumerate(self.predictions):
            if f == fn:
                return i + 1
        return None

    def in_top_k(self, fn: str, k: int) -> bool:
        r = self.rank_of(fn)
        return r is not None and r <= k

    def prob_of(self, fn: str) -> float:
        for f, p, *_ in self.predictions:
            if f == fn:
                return p
        return 0.0


# ---------------------------------------------------------------------------
# FunctionImputer
# ---------------------------------------------------------------------------

class FunctionImputer:
    """
    Predict missing gene functions from genomic neighborhood context.

    Builds inverted lookup tables from a NeighborhoodEstimator so that all
    queries are O(1) or O(|candidates|) — no full-vocabulary scan needed.

    Usage
    -----
    >>> imputer = FunctionImputer(estimator)
    >>> result  = imputer.impute_center(upstream='GeneA', downstream='GeneB')
    >>> print(result.top1, result.top1_prob)
    """

    def __init__(self, estimator: NeighborhoodEstimator):
        self.est = estimator
        self._build_indices()

    # ----------------------------------------------------------------- build

    def _build_indices(self) -> None:
        """
        Pre-compute four inverted tables from the estimator's internals.

        All tables map a key → dict[center_fn, raw_genome_count].
        Stored counts are len(genome_set) (deduped).
        """
        est = self.est

        # ── radius-1 ────────────────────────────────────────────────────────

        # (up, down) → {center: count}   (exact joint)
        r1_pair: dict[tuple, dict[str, int]] = defaultdict(dict)
        for center, pair_dict in est._r1_joint.items():
            for (up, down), gset in pair_dict.items():
                r1_pair[(up, down)][center] = len(gset)
        self._r1_from_pair: dict[tuple, dict[str, int]] = dict(r1_pair)

        # up_fn → {center: count}   (upstream-only marginal)
        r1_up: dict[str, dict[str, int]] = defaultdict(dict)
        for center, up_dict in est._r1_up.items():
            for fn, gset in up_dict.items():
                r1_up[fn][center] = len(gset)
        self._r1_from_up: dict[str, dict[str, int]] = dict(r1_up)

        # down_fn → {center: count}  (downstream-only marginal)
        r1_dn: dict[str, dict[str, int]] = defaultdict(dict)
        for center, dn_dict in est._r1_down.items():
            for fn, gset in dn_dict.items():
                r1_dn[fn][center] = len(gset)
        self._r1_from_down: dict[str, dict[str, int]] = dict(r1_dn)

        # ── radius-2 marginal ────────────────────────────────────────────────

        # (x, y) → {center: count}   (immediate neighbors from quintuples)
        r2_pair: dict[tuple, dict[str, int]] = defaultdict(dict)
        for center, pair_dict in est._r2_marginal.items():
            for (x, y), gset in pair_dict.items():
                r2_pair[(x, y)][center] = len(gset)
        self._r2_from_pair: dict[tuple, dict[str, int]] = dict(r2_pair)

        # ── radius-2 conditional ─────────────────────────────────────────────
        # (f1, x, y, f2) → {center: count}   (full 5-gene context)
        r2_cond: dict[tuple, dict[str, int]] = defaultdict(dict)
        for (center, f1, f2), pair_dict in est._r2_joint.items():
            for (x, y), gset in pair_dict.items():
                r2_cond[(f1, x, y, f2)][center] = len(gset)
        self._r2_from_context: dict[tuple, dict[str, int]] = dict(r2_cond)

    # ----------------------------------------------------------------- queries

    def impute_center(
        self,
        upstream: Optional[str] = None,
        downstream: Optional[str] = None,
        f1: Optional[str] = None,
        f2: Optional[str] = None,
        top_k: Optional[int] = None,
    ) -> ImputationResult:
        """
        Predict the center function given its genomic neighbors.

        Tier priority (highest first):
            r2_conditional  — if all four of (up, down, f1, f2) provided and
                              sum-support ≥ min_support_r2
            r2_marginal     — if (up, down) provided and sum-support ≥ min_r2
            r1_joint        — (up, down) with exact joint counts
            r1_upstream     — up only (marginal)
            r1_downstream   — down only (marginal)

        Parameters
        ----------
        upstream   : immediate upstream function (position X in …X–?–Y…)
        downstream : immediate downstream function (position Y)
        f1         : outer-upstream  (position F1 in F1–X–?–Y–F2)
        f2         : outer-downstream
        top_k      : limit to top-k predictions
        """
        known: dict = {}
        if f1 is not None:         known['f1']         = f1
        if upstream is not None:   known['upstream']   = upstream
        if downstream is not None: known['downstream'] = downstream
        if f2 is not None:         known['f2']         = f2

        candidates: dict[str, int] = {}
        tier = 'unseen'

        # ── tier 1: r2_conditional ──────────────────────────────────────────
        if (f1 is not None and f2 is not None
                and upstream is not None and downstream is not None):
            cond = self._r2_from_context.get((f1, upstream, downstream, f2), {})
            if sum(cond.values()) >= self.est.min_support_r2:
                candidates, tier = cond, 'r2_conditional'

        # ── tier 2: r2_marginal ──────────────────────────────────────────────
        if not candidates and upstream is not None and downstream is not None:
            r2 = self._r2_from_pair.get((upstream, downstream), {})
            if sum(r2.values()) >= self.est.min_support_r2:
                candidates, tier = r2, 'r2_marginal'

        # ── tier 3: r1_joint ────────────────────────────────────────────────
        if not candidates and upstream is not None and downstream is not None:
            r1 = self._r1_from_pair.get((upstream, downstream), {})
            if r1:
                candidates, tier = r1, 'r1_joint'

        # ── tier 4: upstream marginal ────────────────────────────────────────
        if not candidates and upstream is not None:
            r1 = self._r1_from_up.get(upstream, {})
            if r1:
                candidates, tier = r1, 'r1_upstream_only'

        # ── tier 5: downstream marginal ──────────────────────────────────────
        if not candidates and downstream is not None:
            r1 = self._r1_from_down.get(downstream, {})
            if r1:
                candidates, tier = r1, 'r1_downstream_only'

        return self._to_result('center', known, candidates, tier, top_k)

    def impute_upstream(
        self,
        center: str,
        downstream: Optional[str] = None,
        top_k: Optional[int] = None,
    ) -> ImputationResult:
        """
        Predict the upstream neighbor given the center.

        If downstream is also provided, conditions on Y (uses joint table);
        otherwise returns P(X | F) from the marginal.
        """
        known: dict = {'center': center}
        if downstream is not None:
            known['downstream'] = downstream

        candidates: dict[str, int] = {}
        tier = 'unseen'

        if downstream is not None:
            # Scan r1_joint[center] for (X, downstream) pairs
            for (x, y), gset in self.est._r1_joint.get(center, {}).items():
                if y == downstream:
                    candidates[x] = len(gset)
            if candidates:
                tier = 'r1_joint'

        if not candidates:
            denom = len(self.est._r1_center_genomes.get(center, ()))
            if denom >= self.est.min_support_r1:
                candidates = {
                    x: len(g)
                    for x, g in self.est._r1_up.get(center, {}).items()
                }
                tier = 'r1_marginal' if candidates else 'unseen'

        return self._to_result('upstream', known, candidates, tier, top_k)

    def impute_downstream(
        self,
        center: str,
        upstream: Optional[str] = None,
        top_k: Optional[int] = None,
    ) -> ImputationResult:
        """
        Predict the downstream neighbor given the center.

        If upstream is also provided, conditions on X (uses joint table);
        otherwise returns P(Y | F) from the marginal.
        """
        known: dict = {'center': center}
        if upstream is not None:
            known['upstream'] = upstream

        candidates: dict[str, int] = {}
        tier = 'unseen'

        if upstream is not None:
            for (x, y), gset in self.est._r1_joint.get(center, {}).items():
                if x == upstream:
                    candidates[y] = len(gset)
            if candidates:
                tier = 'r1_joint'

        if not candidates:
            denom = len(self.est._r1_center_genomes.get(center, ()))
            if denom >= self.est.min_support_r1:
                candidates = {
                    y: len(g)
                    for y, g in self.est._r1_down.get(center, {}).items()
                }
                tier = 'r1_marginal' if candidates else 'unseen'

        return self._to_result('downstream', known, candidates, tier, top_k)

    # ----------------------------------------------------------------- scoring

    def score_annotation(
        self,
        current_fn: str,
        upstream: Optional[str] = None,
        downstream: Optional[str] = None,
        f1: Optional[str] = None,
        f2: Optional[str] = None,
    ) -> dict:
        """
        Score how well the current annotation fits its neighborhood.

        Runs impute_center() and reports where current_fn ranks in the
        predictions.  High rank (low number) means the annotation is
        consistent with the model; absent or low-ranked means the annotation
        may be wrong *or* the neighborhood is uninformative.

        Returns
        -------
        dict with keys:
            function        : the annotation being scored
            rank            : 1-based rank in predictions (None if not a candidate)
            posterior_prob  : P(current_fn | context); 0.0 if absent
            top1, top1_prob : best prediction
            in_top1/3/5     : bool membership flags
            n_candidates    : total candidates found
            tier            : data tier used
        """
        result = self.impute_center(
            upstream=upstream, downstream=downstream, f1=f1, f2=f2
        )
        rank = result.rank_of(current_fn)
        prob = result.prob_of(current_fn)
        return {
            'function':       current_fn,
            'rank':           rank,
            'posterior_prob': prob,
            'top1':           result.top1,
            'top1_prob':      result.top1_prob,
            'in_top1':        rank == 1,
            'in_top3':        rank is not None and rank <= 3,
            'in_top5':        rank is not None and rank <= 5,
            'n_candidates':   result.n_candidates,
            'tier':           result.tier,
            'total_support':  result.total_support,
        }

    def score_genome_windows(
        self,
        windows: list[tuple],
        target_functions: Optional[set] = None,
    ) -> list[dict]:
        """
        Score the center annotation of every window in a genome.

        For each triple (up, center, down), runs score_annotation() on the
        center and returns enriched rows including the neighborhood context.

        Parameters
        ----------
        windows : list of (up, center, down) triples.
            Quintuples are also accepted; only the inner triple is used.
        target_functions : set or None
            If given, only score windows whose center is in this set.
            Pass {'hypothetical protein', 'none'} to focus on unannotated genes.
        """
        rows = []
        for window in windows:
            if len(window) == 5:
                f1_ctx, up, center, down, f2_ctx = window
            elif len(window) == 3:
                f1_ctx = f2_ctx = None
                up, center, down = window
            else:
                continue

            if target_functions is not None and center not in target_functions:
                continue

            row = self.score_annotation(
                center, upstream=up, downstream=down,
                f1=f1_ctx, f2=f2_ctx
            )
            row['upstream']   = up
            row['downstream'] = down
            rows.append(row)
        return rows

    def suggest_reannotations(
        self,
        windows: list[tuple],
        target_functions: Optional[set] = None,
        min_top1_prob: float = 0.5,
        min_support: int = 5,
    ) -> list[dict]:
        """
        Return windows where the neighborhood strongly suggests a different
        (or more specific) annotation than the current one.

        Filters score_genome_windows() to cases where:
        - The model resolves to a prediction (tier != 'unseen')
        - top1_prob >= min_top1_prob
        - total_support >= min_support
        - current function is NOT ranked #1  (i.e., model disagrees)

        If target_functions is None, includes all unannotated functions
        ('hypothetical protein', 'none') by default.

        Returns rows sorted by top1_prob descending (most confident first).
        """
        if target_functions is None:
            target_functions = {'hypothetical protein', 'none',
                                'hypothetical protein, conserved'}

        rows = self.score_genome_windows(windows, target_functions=target_functions)
        filtered = [
            r for r in rows
            if r['tier'] != 'unseen'
            and r['top1_prob'] >= min_top1_prob
            and r['total_support'] >= min_support
            and not r['in_top1']   # model predicts something different
        ]
        filtered.sort(key=lambda r: -r['top1_prob'])
        return filtered

    # ----------------------------------------------------------------- validation

    def validate(
        self,
        test_windows: list[tuple],
        k_values: Optional[list] = None,
    ) -> dict:
        """
        Mask-and-predict validation on a held-out set of triples.

        For each (up, center, down):
        1. Mask the center (pretend it is unknown).
        2. Run impute_center(up, down).
        3. Check if the true center appears in the top-k predictions.

        This is optimistic when test_windows overlap with training data
        (in-distribution recall), and honest when using a separate test
        corpus or leave-one-genome-out splits.

        Returns
        -------
        dict with:
            total_windows   : total triples evaluated
            seen            : windows the model could resolve
            unseen          : windows with no model coverage
            unseen_rate     : fraction unseen
            recall@k        : fraction of seen windows where truth in top-k
            per_tier        : {tier: {recall@k, total}} breakdown
        """
        if k_values is None:
            k_values = [1, 3, 5, 10]

        hits: dict[int, int] = {k: 0 for k in k_values}
        tier_hits: dict[str, dict[int, int]] = defaultdict(
            lambda: {k: 0 for k in k_values}
        )
        tier_total: dict[str, int] = defaultdict(int)
        total = unseen = 0

        for window in test_windows:
            if len(window) not in (3, 5):
                continue
            if len(window) == 5:
                _, up, center, down, _ = window
            else:
                up, center, down = window

            result = self.impute_center(upstream=up, downstream=down)
            total += 1

            if not result.is_resolved:
                unseen += 1
                continue

            tier_total[result.tier] += 1
            rank = result.rank_of(center)
            for k in k_values:
                if rank is not None and rank <= k:
                    hits[k]             += 1
                    tier_hits[result.tier][k] += 1

        seen = total - unseen
        recalls = {
            f'recall@{k}': hits[k] / seen if seen > 0 else float('nan')
            for k in k_values
        }
        per_tier = {
            t: {
                **{f'recall@{k}': tier_hits[t][k] / n if n > 0 else float('nan')
                   for k in k_values},
                'total': n,
            }
            for t, n in tier_total.items()
        }

        return {
            'total_windows': total,
            'seen':          seen,
            'unseen':        unseen,
            'unseen_rate':   unseen / total if total > 0 else float('nan'),
            **recalls,
            'per_tier':      per_tier,
        }

    # ----------------------------------------------------------------- helpers

    def _to_result(
        self,
        query_type: str,
        known: dict,
        candidates: dict[str, int],
        tier: str,
        top_k: Optional[int],
    ) -> ImputationResult:
        """Normalise raw counts into a sorted ImputationResult."""
        total = sum(candidates.values())

        if not candidates or total == 0:
            return ImputationResult(
                query_type=query_type,
                known=known,
                predictions=[],
                tier='unseen',
                total_support=0,
                n_candidates=0,
            )

        rows = [
            (fn, count / total, count, tier)
            for fn, count in candidates.items()
        ]
        rows.sort(key=lambda r: -r[1])
        if top_k:
            rows = rows[:top_k]

        return ImputationResult(
            query_type=query_type,
            known=known,
            predictions=rows,
            tier=tier,
            total_support=total,
            n_candidates=len(candidates),
        )

    # ----------------------------------------------------------------- index stats

    def index_stats(self) -> dict:
        """Summary statistics about the built inverted indices."""
        return {
            'r1_pair_keys':    len(self._r1_from_pair),
            'r1_up_keys':      len(self._r1_from_up),
            'r1_down_keys':    len(self._r1_from_down),
            'r2_pair_keys':    len(self._r2_from_pair),
            'r2_context_keys': len(self._r2_from_context),
        }
