"""Test PDB fetching."""
import logging
import pytest
from fetch import get_structure


_LOGGER = logging.getLogger(__name__)


@pytest.mark.parametrize("pdb_id", ["1FAS", "7KMK"], ids=str)
def test_fetch(pdb_id):
    """Test fetching of PDB- and CIF-format files."""
    for pdb_fmt in ["CIF", "PDB"]:
        _LOGGER.debug(f"Fetching {pdb_id} in {pdb_fmt} format.")
        get_structure(pdb_id, pdb_fmt)
