"""
Tests for the ModelSEEDAnnotation client.
"""

import pytest
from modelseed_annotation_api import ModelSEEDAnnotationClient


def test_client_initialization():
    """Test that the client can be initialized."""
    client = ModelSEEDAnnotationClient()
    assert client is not None


def test_annotate_sequence_not_implemented():
    """Test that annotate_sequence raises NotImplementedError."""
    client = ModelSEEDAnnotationClient()
    with pytest.raises(NotImplementedError):
        client.annotate_sequence("ATGC")


def test_get_annotation_not_implemented():
    """Test that get_annotation raises NotImplementedError."""
    client = ModelSEEDAnnotationClient()
    with pytest.raises(NotImplementedError):
        client.get_annotation("test_id") 