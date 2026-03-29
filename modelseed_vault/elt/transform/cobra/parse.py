from __future__ import annotations
from typing import Any, Dict, List, Tuple, Union, BinaryIO
import hashlib
from lxml import etree


def parse_model_tag() -> dict:
    res = {}
    pass


def parse_parameters() -> list:
    res = {}
    pass

def parse_fbc_objectives() -> list:
    pass

def parse_fbc_gene_products() -> list:
    pass

def parse_unit_definitions() -> list:
    pass

def parse_groups() -> list:
    pass


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
