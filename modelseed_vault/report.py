from modelseed_vault.core.genome import ProteinSequence
from modelseed_vault.elt.transform.text import extract_ec


def read_bakta(d):
    res = {
        'annotation_bakta_function': None,
    }
    if d:
        if len(d) == 1:
            _e, _n = d[0]
            res['annotation_bakta_function'] = _n['properties']['_value']
        elif len(d) > 1:
            print('!')
    return res


def read_kofam(d):
    res = {
        'annotation_kofam_ko': None,
        'annotation_kofam_full_sequence_score': None,
        'annotation_kofam_full_sequence_evalue': None,
    }
    if d:
        if len(d) == 1:

            _e, _n = d[0]
            # print(_e)
            res['annotation_kofam_ko'] = _n['entry']
            res['annotation_kofam_full_sequence_score'] = _e['properties']['full_sequence_score']
            res['annotation_kofam_full_sequence_evalue'] = _e['properties']['full_sequence_evalue']
        elif len(d) > 1:
            print('!')
    return res


def read_rast(d, q):
    res = {
        'annotation_rast_function': None,
        'annotation_rast_weighted_hit_count_band': None,
        'annotation_rast_weighted_hit_count_n': None,
        'annotation_rast_weighted_hit_count_outlier': None,
        'annotation_rast_weighted_hit_count_value': None,
        'annotation_rast_weighted_hit_count_percentile': None,
        'annotation_rast_hit_count_band': None,
        'annotation_rast_hit_count_n': None,
        'annotation_rast_hit_count_outlier': None,
        'annotation_rast_hit_count_value': None,
        'annotation_rast_hit_count_percentile': None,
    }
    if d:
        if len(d) == 1:

            _e, _n = d[0]
            _ecs = extract_ec(_n['properties']['_value'])
            res['annotation_rast_function'] = _n['properties']['_value']
            res['annotation_rast_function_ec'] = '; '.join(_ecs) if _ecs else None
            res['annotation_rast_weighted_hit_count_value'] = float(_e['properties'].get('weighted_hit_count')) if _e[
                'properties'].get('weighted_hit_count') else None
            res['annotation_rast_hit_count_value'] = int(_e['properties'].get('hit_count')) if _e['properties'].get(
                'hit_count') else None
            if q is not None:
                _any_ec = list(q.get('q'))[0]
                qc_weighted_hit_count, qc_hit_count = q.get('q')[_any_ec]
                res['annotation_rast_weighted_hit_count_band'] = qc_weighted_hit_count.band
                res['annotation_rast_weighted_hit_count_n'] = qc_weighted_hit_count.n
                res['annotation_rast_weighted_hit_count_outlier'] = qc_weighted_hit_count.outlier
                res['annotation_rast_weighted_hit_count_percentile'] = qc_weighted_hit_count.percentile

                res['annotation_rast_hit_count_band'] = qc_hit_count.band
                res['annotation_rast_hit_count_n'] = qc_hit_count.n
                res['annotation_rast_hit_count_outlier'] = qc_hit_count.outlier
                res['annotation_rast_hit_count_percentile'] = qc_hit_count.percentile
        elif len(d) > 1:
            print('!')

    return res


def build_report(genome, feature_quantile, feature_rast, feature_ko, feature_bakta, feature_z, feature_to_reannotation):
    c_data = []
    for f in genome.features:
        h = ProteinSequence(f.seq).hash_value
        _a_rast_q = feature_quantile.get(f.id)
        _a_rast = feature_rast.get(f.id)
        _a_rast_data = read_rast(_a_rast, _a_rast_q)
        _curation = feature_to_reannotation.get(f.id, {})
        _curation_ec = extract_ec(_curation.get('annotation')) if _curation.get('annotation') else None
        record = {
            'feature_id': f.id,
            'sequence': f.seq,
            'annotation_curated': _curation.get('annotation'),
            'annotation_curated_ec': '; '.join(_curation_ec) if _curation_ec else None
        }
        record.update(_a_rast_data)

        _a_kofam = feature_ko.get(f.id)
        record.update(read_kofam(_a_kofam))
        z_score = feature_z.get(f.id.split(':')[1])
        record['neighborhood_conformity_z_score'] = z_score
        _a_bakta = feature_bakta.get(f.id)
        record.update(read_bakta(_a_bakta))
        c_data.append(record)

    return c_data
