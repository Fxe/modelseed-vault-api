"""Quantile/percentile scoring of RAST annotations against a reference sample.

The reference parquet has one row per distinct protein sequence annotated with a
RAST `function`. `hit_count` (raw kmer hits) and `weighted_hit_count` are the
support scores. For a given function those scores form an empirical distribution;
a "foreign" annotation (from some other genome) can then be placed within that
distribution to judge whether its support is low (weak/suspect) or high (typical).
"""

from dataclasses import dataclass

import numpy as np
import pandas as pd

VALUE_COLS = ("hit_count", "weighted_hit_count")

# Tukey/quartile-based ranges. See the module docstring at the bottom and the
# README discussion for *why* these cut points were chosen.
DEFAULT_QS = (0.0, 0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99, 1.0)


def load(path: str, value_cols=VALUE_COLS) -> pd.DataFrame:
    df = pd.read_parquet(path)
    for c in value_cols:
        df[c] = pd.to_numeric(df[c])
    return df


# --------------------------------------------------------------------------- #
# 1. compute quantiles for hit_count AND weighted_hit_count
# --------------------------------------------------------------------------- #
def compute_quantiles(
    df: pd.DataFrame,
    function: str,
    value_cols=VALUE_COLS,
    qs=DEFAULT_QS,
) -> pd.DataFrame:
    """Quantile cut points for a function, one column per score type.

    Returns a DataFrame indexed by quantile (0..1) with columns
    ``hit_count`` and ``weighted_hit_count``. Also tucks the IQR-based
    Tukey outlier fences into ``.attrs`` for the classifier.
    """
    sub = df.loc[df["function"] == function, list(value_cols)]
    if sub.empty:
        raise KeyError(f"no rows for function {function!r}")

    out = pd.DataFrame(
        {c: np.quantile(sub[c].to_numpy(), qs) for c in value_cols},
        index=pd.Index(qs, name="q"),
    )

    fences = {}
    for c in value_cols:
        q1, q3 = np.quantile(sub[c].to_numpy(), [0.25, 0.75])
        iqr = q3 - q1
        fences[c] = {
            "q1": q1,
            "q3": q3,
            "iqr": iqr,
            "low_fence": q1 - 1.5 * iqr,   # Tukey lower outlier fence
            "high_fence": q3 + 1.5 * iqr,  # Tukey upper outlier fence
            "n": len(sub),
        }
    out.attrs["fences"] = fences
    out.attrs["function"] = function
    return out


# --------------------------------------------------------------------------- #
# 2. classify a foreign annotation against the quantiles
# --------------------------------------------------------------------------- #
@dataclass
class Classification:
    function: str
    value_col: str
    value: float
    percentile: float   # 0-100, midrank
    band: str           # coarse percentile band
    outlier: str        # 'low' | 'high' | 'normal' by Tukey fences
    n: int              # sample size backing the call

    def __str__(self) -> str:
        return (
            f"{self.value_col}={self.value:g}  P{self.percentile:.1f}  "
            f"[{self.band}]  outlier={self.outlier}  (n={self.n})"
        )


def _percentile(vals: np.ndarray, value: float) -> float:
    below = np.count_nonzero(vals < value)
    ties = np.count_nonzero(vals == value)
    return 100.0 * (below + 0.5 * ties) / vals.size


def _band(pct: float) -> str:
    # symmetric tail-focused bands (see ranges discussion)
    if pct < 1:
        return "extreme low"
    if pct < 5:
        return "very low"
    if pct < 25:
        return "low"
    if pct <= 75:
        return "typical"
    if pct <= 95:
        return "high"
    if pct <= 99:
        return "very high"
    return "extreme high"


def classify(
    df: pd.DataFrame,
    function: str,
    value: float,
    value_col: str = "weighted_hit_count",
) -> Classification:
    """Place a foreign annotation's score within the function's sample."""
    vals = df.loc[df["function"] == function, value_col].to_numpy()
    if vals.size == 0:
        raise KeyError(f"no rows for function {function!r}")

    pct = _percentile(vals, value)
    q1, q3 = np.quantile(vals, [0.25, 0.75])
    iqr = q3 - q1
    if value < q1 - 1.5 * iqr:
        out = "low"
    elif value > q3 + 1.5 * iqr:
        out = "high"
    else:
        out = "normal"

    return Classification(function, value_col, float(value), pct, _band(pct), out, vals.size)


if __name__ == "__main__":
    PATH = r"C:\Users\filip\Downloads\1.2.1.2.RAST_theseed_org.parquet"
    df = load(PATH)
    fn = "Formate dehydrogenase O beta subunit (EC 1.2.1.2)"

    q = compute_quantiles(df, fn)
    print("Quantiles for:", fn)
    print(q.round(2).to_string())
    print("\nTukey fences:")
    for c, f in q.attrs["fences"].items():
        print(f"  {c:>20}: normal range "
              f"[{max(0, f['low_fence']):.1f} .. {f['high_fence']:.1f}]  (n={f['n']})")

    print("\nClassifying foreign annotations (weighted_hit_count):")
    for v in (5, 20, 120, 350, 600):
        print("  ", classify(df, fn, v, "weighted_hit_count"))
