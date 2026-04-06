from dataclasses import dataclass, field
from typing import Optional
import re


@dataclass
class NodeRecord:
    """A node in the knowledge subgraph."""
    id: str
    type: str
    properties: dict = field(default_factory=dict)

    def to_prompt_str(self) -> str:
        props = ""
        if self.properties:
            props_str = ", ".join(f"{k}={v}" for k, v in self.properties.items())
            props = f" | properties: {props_str}"
        return f"  - [{self.type}] {self.id}{props}"


@dataclass
class EdgeRecord:
    """An edge in the knowledge subgraph."""
    source: str
    relation: str
    target: str
    properties: dict = field(default_factory=dict)

    def to_prompt_str(self) -> str:
        props = ""
        if self.properties:
            props_str = ", ".join(f"{k}={v}" for k, v in self.properties.items())
            props = f" | properties: {props_str}"
        return f"  - {self.source} -[{self.relation}]-> {self.target}{props}"


SYSTEM_PROMPT = """\
You are a computational biology expert evaluating knowledge subgraphs extracted from \
integrated biological databases. Your role is to assess the consistency, completeness, \
and biological plausibility of the graph structure and annotations.

You will receive a subgraph containing nodes (genomic features, genes, reactions, \
metabolites, proteins, etc.) and edges (relationships between them). The subgraph may \
integrate data from multiple sources (COBRA metabolic models, SBML models, genomic \
annotations, sequence databases).

You MUST respond ONLY with a valid XML document matching the schema below. \
Do not include any text outside the XML tags. Do not wrap in markdown code fences."""

EVAL_PROMPT_TEMPLATE = """\
<task>
Evaluate the following biological knowledge subgraph. Identify conflicts between \
data sources, assess biological plausibility, determine the gene function and \
associated reaction(s), and flag any structural or annotation issues in the graph.
</task>

<subgraph>
<nodes>
{nodes_block}
</nodes>

<edges>
{edges_block}
</edges>
</subgraph>

{context_block}

<instructions>
1. CONFLICTS: Identify any conflicts or inconsistencies between data sources, \
annotations, or model representations. For each conflict, specify the type, \
the nodes/edges involved, severity (critical/warning/info), and a description.

2. CONCLUSION: Based on the subgraph, determine:
   a. The biological function of the central gene/feature.
   b. The reaction(s) associated with this gene (there may be 1 or more). \
For each reaction, provide the reaction ID, a human-readable name, \
directionality, substrates, products, and your confidence level (0.0-1.0).

3. GRAPH ISSUES: Flag any structural problems with the knowledge graph itself, such as:
   - Missing expected edges (e.g., gene without protein, reaction without metabolites)
   - Orphan nodes (nodes with no connecting edges)
   - Redundant or duplicate representations
   - Inconsistent identifier schemes
   - Missing annotations or cross-references

4. EVIDENCE SUMMARY: Provide a brief summary of the evidence supporting your conclusions.

5. RECOMMENDATIONS: Suggest specific actions to resolve conflicts or improve the subgraph.
</instructions>

<response_format>
Respond ONLY with the following XML structure:

<evaluation>
  <conflicts>
    <conflict>
      <type>annotation_mismatch | stoichiometry_error | missing_link | \
cross_reference_conflict | model_disagreement | identifier_inconsistency</type>
      <severity>critical | warning | info</severity>
      <description>Human-readable description of the conflict</description>
      <nodes_involved>comma-separated list of node IDs</nodes_involved>
      <suggested_resolution>How to resolve this conflict</suggested_resolution>
    </conflict>
    <!-- repeat for each conflict found, or leave empty if none -->
  </conflicts>

  <conclusion>
    <gene_function>Description of the gene's biological function</gene_function>
    <reactions>
      <reaction>
        <id>Reaction identifier from the graph</id>
        <name>Human-readable reaction name (e.g., EC number or common name)</name>
        <direction>reversible | irreversible_forward | irreversible_reverse</direction>
        <substrates>comma-separated list of substrate names/IDs</substrates>
        <products>comma-separated list of product names/IDs</products>
        <confidence>0.0 to 1.0</confidence>
        <evidence>Brief justification for this reaction assignment</evidence>
      </reaction>
      <!-- repeat for each reaction -->
    </reactions>
  </conclusion>

  <graph_issues>
    <issue>
      <type>missing_edge | orphan_node | redundant_representation | \
identifier_inconsistency | missing_annotation | schema_violation</type>
      <severity>critical | warning | info</severity>
      <description>Description of the graph issue</description>
      <affected_elements>comma-separated list of affected node/edge IDs</affected_elements>
      <suggested_fix>How to fix this issue</suggested_fix>
    </issue>
    <!-- repeat for each issue found -->
  </graph_issues>

  <evidence_summary>
    Brief summary of the overall evidence quality and consistency
  </evidence_summary>

  <overall_confidence>0.0 to 1.0</overall_confidence>

  <recommendations>
    <recommendation>Specific action item</recommendation>
    <!-- repeat -->
  </recommendations>
</evaluation>
</response_format>"""


def build_prompt(
        nodes: list[NodeRecord],
        edges: list[EdgeRecord],
        context: Optional[str] = None
) -> tuple[str, str]:
    """
    Build the system prompt and user prompt for subgraph evaluation.

    Returns:
        (system_prompt, user_prompt)
    """
    nodes_block = "\n".join(n.to_prompt_str() for n in nodes)
    edges_block = "\n".join(e.to_prompt_str() for e in edges)

    context_block = ""
    if context:
        context_block = f"\n<additional_context>\n{context}\n</additional_context>\n"

    user_prompt = EVAL_PROMPT_TEMPLATE.format(
        nodes_block=nodes_block,
        edges_block=edges_block,
        context_block=context_block
    )

    return SYSTEM_PROMPT, user_prompt


# ─────────────────────────────────────────────
# Response Parser
# ─────────────────────────────────────────────

def _extract_tag(xml_str: str, tag: str) -> Optional[str]:
    """Extract content of a single XML tag."""
    pattern = rf"<{tag}>(.*?)</{tag}>"
    match = re.search(pattern, xml_str, re.DOTALL)
    return match.group(1).strip() if match else None


def _extract_all_tags(xml_str: str, tag: str) -> list[str]:
    """Extract content of all occurrences of a tag."""
    pattern = rf"<{tag}>(.*?)</{tag}>"
    return [m.strip() for m in re.findall(pattern, xml_str, re.DOTALL)]


def _parse_conflict(conflict_xml: str) -> dict:
    """Parse a single <conflict> block."""
    return {
        "type": _extract_tag(conflict_xml, "type") or "unknown",
        "severity": _extract_tag(conflict_xml, "severity") or "info",
        "description": _extract_tag(conflict_xml, "description") or "",
        "nodes_involved": [
            n.strip()
            for n in (_extract_tag(conflict_xml, "nodes_involved") or "").split(",")
            if n.strip()
        ],
        "suggested_resolution": _extract_tag(conflict_xml, "suggested_resolution") or "",
    }


def _parse_reaction(reaction_xml: str) -> dict:
    """Parse a single <reaction> block."""
    return {
        "id": _extract_tag(reaction_xml, "id") or "",
        "name": _extract_tag(reaction_xml, "name") or "",
        "direction": _extract_tag(reaction_xml, "direction") or "unknown",
        "substrates": [
            s.strip()
            for s in (_extract_tag(reaction_xml, "substrates") or "").split(",")
            if s.strip()
        ],
        "products": [
            p.strip()
            for p in (_extract_tag(reaction_xml, "products") or "").split(",")
            if p.strip()
        ],
        "confidence": float(_extract_tag(reaction_xml, "confidence") or 0.0),
        "evidence": _extract_tag(reaction_xml, "evidence") or "",
    }


def _parse_issue(issue_xml: str) -> dict:
    """Parse a single <issue> block."""
    return {
        "type": _extract_tag(issue_xml, "type") or "unknown",
        "severity": _extract_tag(issue_xml, "severity") or "info",
        "description": _extract_tag(issue_xml, "description") or "",
        "affected_elements": [
            e.strip()
            for e in (_extract_tag(issue_xml, "affected_elements") or "").split(",")
            if e.strip()
        ],
        "suggested_fix": _extract_tag(issue_xml, "suggested_fix") or "",
    }


@dataclass
class SubgraphEvalResult:
    """Parsed result from LLM evaluation."""
    gene_function: str
    reactions: list  # list of dicts with id, name, direction, confidence
    conflicts: list  # list of dicts with type, description, severity, nodes_involved
    issues: list     # list of dicts with type, description, severity, affected_elements
    confidence: float
    evidence_summary: str
    recommendations: list  # list of strings


def parse_evaluation_response(response_text: str) -> SubgraphEvalResult:
    """
    Parse the XML response from the LLM into a structured result.

    Args:
        response_text: Raw XML string from LLM response

    Returns:
        SubgraphEvalResult with all parsed fields

    Raises:
        ValueError: If the response cannot be parsed
    """
    # Strip any markdown fences if present
    text = response_text.strip()
    if text.startswith("```"):
        text = re.sub(r"^```\w*\n?", "", text)
        text = re.sub(r"\n?```$", "", text)
        text = text.strip()

    # Validate top-level structure
    if "<evaluation>" not in text:
        raise ValueError("Response missing <evaluation> root element")

    eval_block = _extract_tag(text, "evaluation")
    if not eval_block:
        raise ValueError("Could not extract <evaluation> block")

    # Parse conflicts
    conflicts = [_parse_conflict(c) for c in _extract_all_tags(eval_block, "conflict")]

    # Parse conclusion
    gene_function = _extract_tag(eval_block, "gene_function") or "Unknown"
    reactions = [_parse_reaction(r) for r in _extract_all_tags(eval_block, "reaction")]

    # Parse graph issues
    issues = [_parse_issue(i) for i in _extract_all_tags(eval_block, "issue")]

    # Parse metadata
    confidence_str = _extract_tag(eval_block, "overall_confidence") or "0.0"
    try:
        confidence = float(confidence_str)
    except ValueError:
        confidence = 0.0

    evidence_summary = _extract_tag(eval_block, "evidence_summary") or ""
    recommendations = _extract_all_tags(eval_block, "recommendation")

    return SubgraphEvalResult(
        gene_function=gene_function,
        reactions=reactions,
        conflicts=conflicts,
        issues=issues,
        confidence=confidence,
        evidence_summary=evidence_summary,
        recommendations=recommendations,
    )
