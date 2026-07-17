"""
Unit and regression test for the openawsem package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import openawsem


def test_openawsem_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "openawsem" in sys.modules
