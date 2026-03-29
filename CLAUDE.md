# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Install in development mode
pip install -e .

# Run all tests
pytest

# Run a single test
pytest tests/test_client.py::test_client_initialization

# Format code
black modelseed_vault/
isort modelseed_vault/

# Type check
mypy modelseed_vault/
```

## Architecture

This is the **ModelSEED Vault Python API** — a bioinformatics library for interacting with ModelSEED graph-based annotation services. The package (both the importable module and the installable) is `modelseed_vault`.

### Dual-backend data layer

The library has two distinct persistence backends that are used at different levels:

1. **`Vault` (`vault.py`)** — HTTP REST client for a running Vault API server (default `http://192.168.1.22:12022/`). It manages graph nodes and edges via `/graph/node/` and `/graph/edge/` endpoints, and can fetch COBRA metabolic models via `/cobra/model/`.

2. **`Neo4jDAO` (`dao_neo4j.py`)** — Direct Neo4j driver connection for Cypher queries. Nodes use a `key` property as the domain identifier. `GraphNode` and its subclasses (`GraphNodeProtein`, `GraphNodeRastFunction`, `GraphNodeRastExecution`) wrap Neo4j nodes with fetch/create/delete lifecycle.

### Core domain objects

- **`Node` (`core/node.py`)** — Generic labeled graph node with `key`, `label`, and `data` dict. Keys have spaces replaced with underscores. Identity is `label/key`.
- **`AnnotationFunction` (`core/base.py`)** — Protein function annotation with synonyms, sub-functions, and function groups. Serializes to/from JSON.
- **`ModelSEEDAnnotationClient` (`client.py`)** — High-level client wrapping `Neo4jDAO`. Currently implements RAST annotation ingestion (`get_rast_annotation`) and genome set retrieval; most annotation query methods raise `NotImplementedError`.

### ETL pipeline

`elt/sbml/` contains an SBML metabolic model parser and transformer:
- `parse.py` — XML parsing of SBML elements (species, compartments, reactions) with provenance
- `etl_sbml.py` — `ETLSBML` transforms parsed SBML into `Node` objects with labels `SBMLModel`, `SBMLSpecies`, `SBMLCompartment`, `SBMLReaction`, and builds typed edges (`has_reactant`, `has_product`, `has_sbml_*`)

### Other modules

- `seq_store_hdf5.py` / `seq_store_mongo.py` — Sequence storage backends (HDF5 and MongoDB)
- `biodb/uniprot/` — UniProt database parsers
- `api_curation.py` — Curation API helpers
- `object_factory.py` — Stub factory for constructing domain objects from Neo4j node labels (incomplete)
