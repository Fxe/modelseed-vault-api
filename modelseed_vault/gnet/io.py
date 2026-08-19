import json
from collections import defaultdict
from pathlib import Path

from modelseed_vault.gnet.genome_net import NeighborhoodEstimator


def save(estimator: NeighborhoodEstimator, path: str | Path) -> None:
    tk = lambda t: '\t'.join(t)
    s2l = lambda v: list(v)

    data = {
        'total_genomes': estimator.total_genomes,
        'min_support_r2': estimator.min_support_r2,
        'min_support_r1': estimator.min_support_r1,
        '_r1_joint': {
            k: {tk(ik): s2l(v) for ik, v in iv.items()}
            for k, iv in estimator._r1_joint.items()
        },
        '_r1_up': {
            k: {ik: s2l(v) for ik, v in iv.items()}
            for k, iv in estimator._r1_up.items()
        },
        '_r1_down': {
            k: {ik: s2l(v) for ik, v in iv.items()}
            for k, iv in estimator._r1_down.items()
        },
        '_r1_center_genomes': {k: s2l(v) for k, v in estimator._r1_center_genomes.items()},
        '_r2_joint': {
            tk(k): {tk(ik): s2l(v) for ik, v in iv.items()}
            for k, iv in estimator._r2_joint.items()
        },
        '_r2_context_genomes': {tk(k): s2l(v) for k, v in estimator._r2_context_genomes.items()},
        '_r2_marginal': {
            k: {tk(ik): s2l(v) for ik, v in iv.items()}
            for k, iv in estimator._r2_marginal.items()
        },
        '_r2_center_genomes': {k: s2l(v) for k, v in estimator._r2_center_genomes.items()},
    }

    with open(path, 'w') as fh:
        json.dump(data, fh)


def load(path: str | Path) -> NeighborhoodEstimator:
    with open(path, 'rb') as fh:
        is_json = fh.read(1) == b'{'
    if not is_json:
        import pickle
        with open(path, 'rb') as fh:
            return pickle.load(fh)
    with open(path) as fh:
        data = json.load(fh)

    tk = lambda s: tuple(s.split('\t'))

    est = NeighborhoodEstimator(
        total_genomes=data['total_genomes'],
        min_support_r2=data['min_support_r2'],
        min_support_r1=data['min_support_r1'],
    )

    est._r1_joint = defaultdict(lambda: defaultdict(set), {
        k: defaultdict(set, {tk(ik): set(v) for ik, v in iv.items()})
        for k, iv in data['_r1_joint'].items()
    })
    est._r1_up = defaultdict(lambda: defaultdict(set), {
        k: defaultdict(set, {ik: set(v) for ik, v in iv.items()})
        for k, iv in data['_r1_up'].items()
    })
    est._r1_down = defaultdict(lambda: defaultdict(set), {
        k: defaultdict(set, {ik: set(v) for ik, v in iv.items()})
        for k, iv in data['_r1_down'].items()
    })
    est._r1_center_genomes = defaultdict(set, {
        k: set(v) for k, v in data['_r1_center_genomes'].items()
    })
    est._r2_joint = defaultdict(lambda: defaultdict(set), {
        tk(k): defaultdict(set, {tk(ik): set(v) for ik, v in iv.items()})
        for k, iv in data['_r2_joint'].items()
    })
    est._r2_context_genomes = defaultdict(set, {
        tk(k): set(v) for k, v in data['_r2_context_genomes'].items()
    })
    est._r2_marginal = defaultdict(lambda: defaultdict(set), {
        k: defaultdict(set, {tk(ik): set(v) for ik, v in iv.items()})
        for k, iv in data['_r2_marginal'].items()
    })
    est._r2_center_genomes = defaultdict(set, {
        k: set(v) for k, v in data['_r2_center_genomes'].items()
    })

    return est
