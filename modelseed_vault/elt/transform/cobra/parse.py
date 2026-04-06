from __future__ import annotations
from typing import Any, Dict, List, Tuple, Union, BinaryIO
import hashlib
import re
from lxml import etree


def parse_model_tag(fh) -> dict:
    """
    Return the attributes of the top-level <model> element as a dict
    (Clark-notation keys converted to prefixed form, e.g. fbc:strict).
    """
    raw_bytes = fh.read()
    parser = etree.XMLParser(remove_blank_text=False, recover=False, huge_tree=True)
    root = etree.fromstring(raw_bytes, parser=parser)
    nodes = root.xpath("//*[local-name()='model']")
    if not nodes:
        return {}
    return get_node_attrib(nodes[0])


def parse_parameters(fh) -> list:
    """
    Return a list of dicts for every <parameter> element inside
    <listOfParameters>, each containing all XML attributes with
    prefixed keys (e.g. sboTerm, id, value, constant).
    """
    result = parse_elements_with_provenance(fh, "parameter",
                                            xpath="//*[local-name()='listOfParameters']"
                                                  "/*[local-name()='parameter']")
    return result["elements"]


def parse_fbc_objectives(fh) -> list:
    """
    Return a list of dicts for every <fbc:objective> element.
    Each dict contains the objective attributes plus a
    '_flux_objectives' key holding a list of fluxObjective attribute dicts.
    """
    raw_bytes = fh.read()
    parser = etree.XMLParser(remove_blank_text=False, recover=False, huge_tree=True)
    root = etree.fromstring(raw_bytes, parser=parser)

    objectives = root.xpath("//*[local-name()='objective']")
    results = []
    for obj in objectives:
        d = get_node_attrib(obj)
        flux_objs = obj.xpath("*[local-name()='listOfFluxObjectives']/*[local-name()='fluxObjective']")
        d["_flux_objectives"] = [get_node_attrib(fo) for fo in flux_objs]
        results.append(d)
    return results


def parse_fbc_gene_products(fh) -> list:
    """
    Return a list of dicts for every <fbc:geneProduct> element,
    each containing all XML attributes with prefixed keys
    (e.g. fbc:id, fbc:label, metaid, sboTerm).
    """
    raw_bytes = fh.read()
    parser = etree.XMLParser(remove_blank_text=False, recover=False, huge_tree=True)
    root = etree.fromstring(raw_bytes, parser=parser)

    gene_products = root.xpath("//*[local-name()='geneProduct']")
    return [get_node_attrib(gp) for gp in gene_products]


def parse_unit_definitions(fh) -> list:
    """
    Return a list of dicts for every <unitDefinition> element.
    Each dict contains the unitDefinition attributes plus a
    '_units' key holding a list of <unit> attribute dicts.
    """
    raw_bytes = fh.read()
    parser = etree.XMLParser(remove_blank_text=False, recover=False, huge_tree=True)
    root = etree.fromstring(raw_bytes, parser=parser)

    unit_defs = root.xpath("//*[local-name()='unitDefinition']")
    results = []
    for ud in unit_defs:
        d = get_node_attrib(ud)
        units = ud.xpath("*[local-name()='listOfUnits']/*[local-name()='unit']")
        d["_units"] = [get_node_attrib(u) for u in units]
        results.append(d)
    return results


def parse_groups(fh) -> list:
    """
    Return a list of dicts for every <groups:group> element.
    Each dict contains the group attributes plus a
    '_members' key holding a list of idRef strings from <groups:member>.
    """
    raw_bytes = fh.read()
    parser = etree.XMLParser(remove_blank_text=False, recover=False, huge_tree=True)
    root = etree.fromstring(raw_bytes, parser=parser)

    groups = root.xpath("//*[local-name()='group']")
    results = []
    for grp in groups:
        d = get_node_attrib(grp)
        members = grp.xpath("*[local-name()='listOfMembers']/*[local-name()='member']")
        d["_members"] = [get_node_attrib(m) for m in members]
        results.append(d)
    return results


def parse_species_list(xml: str, list_tag: str = "listOfReactants", wrapper: str = None) -> list:
    wrapped = wrapper.format(xml) if wrapper else xml
    root = etree.fromstring(wrapped.encode("utf-8"))
    nodes = root.xpath(f"//*[local-name()='{list_tag}']/*[local-name()='speciesReference']")
    return [dict(node.attrib) for node in nodes]


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _build_line_offsets(text: str) -> List[int]:
    offsets = [0]
    for i, ch in enumerate(text):
        if ch == "\n":
            offsets.append(i + 1)
    return offsets


def _find_tag_end(text: str, start: int) -> int:
    in_quote = None
    i = start
    while i < len(text):
        c = text[i]
        if in_quote:
            if c == in_quote:
                in_quote = None
        else:
            if c in ("'", '"'):
                in_quote = c
            elif c == ">":
                return i
        i += 1
    raise ValueError("Malformed XML: could not find end of start tag ('>').")


def _slice_element_raw_xml(text: str, line_offsets: List[int], start_line: int, tag: str) -> Tuple[str, int, int]:
    """
    Slice exact raw XML for the <tag ...> element that starts near start_line.
    Handles both normal and self-closing tags.
    Returns: (raw_xml, start_char, end_char)
    """
    pos = line_offsets[start_line - 1] if 1 <= start_line <= len(line_offsets) else 0

    start = text.find(f"<{tag}", pos)
    if start == -1:
        raise ValueError(f"Could not locate '<{tag}' start tag near line {start_line}.")

    tag_end = _find_tag_end(text, start)
    start_tag = text[start: tag_end + 1]

    # Self-closing: <tag .../>
    if start_tag.rstrip().endswith("/>"):
        end = tag_end + 1
        return text[start:end], start, end

    # Find matching closing tag — track nesting depth to handle nested same-name tags
    depth = 1
    search_pos = tag_end + 1
    while depth > 0:
        next_open = text.find(f"<{tag}", search_pos)
        next_close = text.find(f"</{tag}", search_pos)

        if next_close == -1:
            raise ValueError(f"Malformed XML: missing '</{tag}>' closing tag.")

        if next_open != -1 and next_open < next_close:
            # Another opening tag of the same type comes first — go deeper
            inner_tag_end = _find_tag_end(text, next_open)
            inner_start_tag = text[next_open: inner_tag_end + 1]
            if not inner_start_tag.rstrip().endswith("/>"):
                depth += 1
            search_pos = inner_tag_end + 1
        else:
            depth -= 1
            close_end = text.find(">", next_close)
            if close_end == -1:
                raise ValueError(f"Malformed XML: missing '>' after '</{tag}'.")
            search_pos = close_end + 1
            if depth == 0:
                return text[start: close_end + 1], start, close_end + 1


def _clark_to_prefixed(attrib: dict, nsmap: dict) -> dict:
    """
    Convert Clark-notation attribute keys like
    {http://www.sbml.org/sbml/level3/version1/fbc/version2}lowerFluxBound
    to prefixed form like fbc:lowerFluxBound.

    Falls back to the bare local name if no prefix is found.
    """
    # Build reverse map: URI -> prefix (skip default namespace where prefix is None)
    uri_to_prefix = {uri: pfx for pfx, uri in nsmap.items() if pfx is not None}

    result = {}
    for key, value in attrib.items():
        if key.startswith("{"):
            uri, localname = key[1:].split("}", 1)
            prefix = uri_to_prefix.get(uri)
            new_key = f"{prefix}:{localname}" if prefix else localname
        else:
            new_key = key
        result[new_key] = value
    return result


def get_node_attrib(node) -> Dict[str, Any]:
    return _clark_to_prefixed(node.attrib, node.nsmap)


def parse_elements_with_provenance(fh, tag: str, xpath: str | None = None) -> Dict[str, Any]:
    """
    Generic evidence-grade extraction for any XML element tag.

    Args:
        fh:    File-like object (binary) to read from.
        tag:   Local tag name to search for, e.g. 'species', 'compartment', 'reaction'.
        xpath: Optional lxml XPath to locate nodes. Defaults to searching the entire
               document for any element with the given local name:
               "//*[local-name()='<tag>']"

    Returns a dict with:
        - file_sha256: hex digest of the raw file bytes
        - <tag>s: list of dicts, each with all XML attributes plus:
            _row      : 1-based source line number
            _xpath    : structural XPath locator
            _raw_xml  : exact raw XML slice from the original file
            _start    : char offset (start) in decoded text
            _end      : char offset (end) in decoded text
    """
    raw_bytes = fh.read()
    raw_text = raw_bytes.decode("utf-8")
    line_offsets = _build_line_offsets(raw_text)

    parser = etree.XMLParser(remove_blank_text=False, recover=False, huge_tree=True)
    root = etree.fromstring(raw_bytes, parser=parser)

    if xpath is None:
        xpath = f"//*[local-name()='{tag}']"

    nodes = root.xpath(xpath)
    results: List[Dict[str, Any]] = []
    tree = root.getroottree()

    for node in nodes:
        #d: Dict[str, Any] = dict(node.attrib)
        d: Dict[str, Any] = get_node_attrib(node)
        d["_row"] = node.sourceline
        d["_xpath"] = tree.getpath(node)
        raw_xml, start, end = _slice_element_raw_xml(raw_text, line_offsets, node.sourceline or 1, tag)
        d["_raw_xml"] = raw_xml
        d["_start"] = start
        d["_end"] = end
        results.append(d)

    return {"file_sha256": _sha256(raw_bytes), "elements": results}

def _parse_species_references(node, ns_tag: str) -> List[List]:
    """
    Extract speciesReference entries from a listOfReactants or listOfProducts node.
    Returns list of [species_id, stoichiometry_or_None].
    """
    results = []
    for ref in node:
        local = etree.QName(ref.tag).localname if ref.tag else ""
        if local != "speciesReference":
            continue
        species = ref.attrib.get("species", None)
        stoich_raw = ref.attrib.get("stoichiometry", None)
        stoich = stoich_raw if stoich_raw is not None else None
        results.append([species, stoich])
    return results


def _raw_slice(text: str, node) -> str:
    """
    Return the raw XML text of a single lxml element by its sourceline.
    Uses the element's own tag local name for slicing.
    """
    tag = etree.QName(node.tag).localname
    line_offsets = _build_line_offsets(text)
    raw_xml, _, _ = _slice_element_raw_xml(text, line_offsets, node.sourceline or 1, tag)
    return raw_xml


def parse_reactions_with_provenance(fh, xpath: str | None = None) -> Dict[str, Any]:
    """
    Reaction-specific evidence-grade extraction.

    Each reaction record contains:
        - all XML attributes (metaid, id, name, reversible, ...)
        - _row          : 1-based source line of the <reaction> tag
        - _xpath        : structural lxml XPath locator
        - _raw_xml      : full raw XML of the <reaction> element
        - _start        : char offset (start) in decoded text
        - _end          : char offset (end) in decoded text
        - _raw_notes    : raw XML of the <notes> child, or None
        - _raw_annotation : raw XML of the <annotation> child, or None
        - _raw_stoichiometry : raw XML of <listOfReactants> + <listOfProducts>, or None
        - _reactants    : [ [species_id, stoichiometry_or_None], ... ]
        - _products     : [ [species_id, stoichiometry_or_None], ... ]
    """
    raw_bytes = fh.read()
    raw_text = raw_bytes.decode("utf-8")
    line_offsets = _build_line_offsets(raw_text)

    parser = etree.XMLParser(remove_blank_text=False, recover=False, huge_tree=True)
    root = etree.fromstring(raw_bytes, parser=parser)

    if xpath is None:
        xpath = "//*[local-name()='reaction']"

    nodes = root.xpath(xpath)
    tree = root.getroottree()
    results: List[Dict[str, Any]] = []

    for rxn in nodes:
        d: Dict[str, Any] = get_node_attrib(rxn)
        d["_row"] = rxn.sourceline
        d["_xpath"] = tree.getpath(rxn)

        # Full raw XML + char offsets
        raw_xml, start, end = _slice_element_raw_xml(raw_text, line_offsets, rxn.sourceline or 1, "reaction")
        d["_raw_xml"] = raw_xml
        d["_start"] = start
        d["_end"] = end

        # Per-child extraction — keyed by local-name
        children = {etree.QName(child.tag).localname: child for child in rxn}

        # <notes>
        #notes_node = children.get("notes")
        #d["_raw_notes"] = _raw_slice(raw_text, notes_node) if notes_node is not None else None

        # <annotation>
        #ann_node = children.get("annotation")
        #d["_raw_annotation"] = _raw_slice(raw_text, ann_node) if ann_node is not None else None

        # <listOfReactants> + <listOfProducts> — store combined raw XML block
        lor = children.get("listOfReactants")
        lop = children.get("listOfProducts")

        """
        if lor is not None or lop is not None:
            parts = []
            if lor is not None:
                parts.append(_raw_slice(raw_text, lor))
            if lop is not None:
                parts.append(_raw_slice(raw_text, lop))
            d["_raw_stoichiometry"] = "\n".join(parts)
        else:
            d["_raw_stoichiometry"] = None
        """

        d["_reactants"] = _parse_species_references(lor, "speciesReference") if lor is not None else []
        d["_products"] = _parse_species_references(lop, "speciesReference") if lop is not None else []

        results.append(d)

    return {"file_sha256": _sha256(raw_bytes), "elements": results}


def _parse_gpa_node(node) -> List[List[str]]:
    """
    Recursively parse an fbc gene product association node.
    Returns a list of complexes; each complex is a list of gene fbc:id strings.
    - geneProductRef  → [[gene_id]]
    - and             → [all children flattened into one complex]
    - or              → one complex per child
    """
    local = etree.QName(node.tag).localname
    if local == "geneProductRef":
        gene_id = next(
            (v for k, v in node.attrib.items() if k.endswith("}geneProduct") or k == "geneProduct"),
            None,
        )
        return [[gene_id]] if gene_id else []
    if local == "and":
        genes: List[str] = []
        for child in node:
            for grp in _parse_gpa_node(child):
                genes.extend(grp)
        return [genes] if genes else []
    if local == "or":
        result: List[List[str]] = []
        for child in node:
            result.extend(_parse_gpa_node(child))
        return result
    return []


def parse_gene_associations(fh) -> Dict[str, List[List[str]]]:
    """
    Extract fbc:geneProductAssociation for every reaction in an SBML file.

    Returns:
        dict mapping reaction id → list of complexes,
        where each complex is a list of gene fbc:id strings (AND group).
        OR relationships between complexes are represented as separate list entries.
    """
    raw_bytes = fh.read()
    parser = etree.XMLParser(remove_blank_text=False, recover=False, huge_tree=True)
    root = etree.fromstring(raw_bytes, parser=parser)

    result: Dict[str, List[List[str]]] = {}
    for rxn in root.xpath("//*[local-name()='reaction']"):
        rxn_id = rxn.attrib.get("id", "")
        gpas = rxn.xpath("*[local-name()='geneProductAssociation']")
        if not gpas:
            result[rxn_id] = []
            continue
        children = list(gpas[0])
        result[rxn_id] = _parse_gpa_node(children[0]) if children else []
    return result


# ── Notes-based GPR scanner ────────────────────────────────────────────────────

def _tokenize_gpr(expr: str) -> List[str]:
    return [t for t in re.findall(r'\(|\)|and|or|[^\s()]+', expr, re.IGNORECASE) if t]


def _gpr_parse(tokens: List[str], pos: int) -> Tuple[List[List[str]], int]:
    """OR level — returns list of complexes (isozyme alternatives)."""
    complexes, pos = _gpr_and(tokens, pos)
    while pos < len(tokens) and tokens[pos].lower() == 'or':
        pos += 1
        right, pos = _gpr_and(tokens, pos)
        complexes = complexes + right
    return complexes, pos


def _gpr_and(tokens: List[str], pos: int) -> Tuple[List[List[str]], int]:
    """AND level — distributes AND over OR sub-expressions."""
    result, pos = _gpr_atom(tokens, pos)
    while pos < len(tokens) and tokens[pos].lower() == 'and':
        pos += 1
        right, pos = _gpr_atom(tokens, pos)
        result = [l + r for l in result for r in right]
    return result, pos


def _gpr_atom(tokens: List[str], pos: int) -> Tuple[List[List[str]], int]:
    """Atom: parenthesized sub-expression or single gene id."""
    if pos >= len(tokens):
        return [], pos
    if tokens[pos] == '(':
        pos += 1
        result, pos = _gpr_parse(tokens, pos)
        if pos < len(tokens) and tokens[pos] == ')':
            pos += 1
        return result, pos
    return [[tokens[pos]]], pos + 1


def scan_gpr_nodes(reaction_xml: str) -> List[List[str]]:
    """
    Extract GPR complexes from a reaction's raw XML by scanning the legacy
    ``<notes><html:p>GENE ASSOCIATION: ...</html:p></notes>`` pattern.

    Parses the boolean expression using AND/OR/parentheses:
    - AND  → genes belong to the same protein complex
    - OR   → separate isozyme alternatives

    Returns a list of complexes in the same format as ``_complexes``:
        [ [gene_id, ...], ... ]
    Returns ``[]`` if no GENE ASSOCIATION note is found or the expression is empty.
    """
    m = re.search(r'GENE[_ ]ASSOCIATION\s*:\s*([^\n<]+)', reaction_xml, re.IGNORECASE)
    if not m:
        return []
    expr = m.group(1).strip()
    if not expr or expr.lower() == 'none':
        return []
    tokens = _tokenize_gpr(expr)
    if not tokens:
        return []
    complexes, _ = _gpr_parse(tokens, 0)
    return [c for c in complexes if c]
