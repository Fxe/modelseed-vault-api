import re


def extract_ec1(s: str) -> list[str]:
    match = re.search(r'\[EC:([^\]]+)\]', s)
    if not match:
        return []
    return re.findall(r'[\d-]+\.[\d-]+\.[\d-]+\.[\d-]+', match.group(1))


def extract_ec2(s: str) -> list[str]:
    return re.findall(r'\bEC[: ]([\d-]+\.[\d-]+\.[\d-]+\.[\d-]+)', s)


def extract_ec3(s: str) -> list[str]:
    # Handles: [EC:3.1.3.16 3.1.3.48], (EC 2.1.1.100), (EC 3.2.1.n1)
    bracket_match = re.search(r'\[EC:([^\]]+)\]', s)
    if bracket_match:
        return re.findall(r'[\dn-]+\.[\dn-]+\.[\dn-]+\.[\dn-]+', bracket_match.group(1))
    return re.findall(r'\bEC[: ]([\dn-]+\.[\dn-]+\.[\dn-]+\.[\dn-]+)', s)


def extract_ec(text: str) -> list[str] | None:
    v = text
    ecs1 = extract_ec1(v)
    ecs2 = extract_ec2(v)
    ecs3 = extract_ec3(v)
    if len(ecs1) > 0:
        return ecs1
    elif len(ecs2) > 0:
        return ecs2
    elif len(ecs3) > 0:
        return ecs3
    else:
        return None


def auto_ec(vault):
    annotation_values = vault.list_nodes('FunctionalAnnotation')

    maybe_has_ec = {}
    for o in annotation_values:
        _v = o['properties']['_value']
        if '.' in _v:
            maybe_has_ec[o['entry']] = _v

    found_ec = {}
    maybe_not_ec = {}
    for k in maybe_has_ec:
        text = maybe_has_ec[k]
        ecs = extract_ec(text)
        if ecs is not None and len(ecs) > 0:
            found_ec[k] = ecs
        else:
            maybe_not_ec[k] = text

    all_ecs = set()
    for ecs in found_ec.values():
        all_ecs |= set(ecs)

    ref_to_eid = {}
    refs = [['ECNumber', x] for x in all_ecs]
    eids = vault.query_eid(refs)
    for o in eids:
        ref_to_eid[(o['type'], o['key'])] = o['elementId']
    refs = [['FunctionalAnnotation', x] for x in found_ec]
    eids = vault.query_eid(refs)
    for o in eids:
        ref_to_eid[(o['type'], o['key'])] = o['elementId']

    payload = []
    for k, l_ec in found_ec.items():
        for ec in l_ec:
            edge = [ref_to_eid[('FunctionalAnnotation', k)], 'contains', ref_to_eid[('ECNumber', ec)], {}]
            payload.append(edge)

    return maybe_has_ec, found_ec, maybe_not_ec, payload
