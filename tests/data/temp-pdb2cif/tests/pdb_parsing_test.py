"""Test old-format PDB parsing."""
import logging
import pytest
from pdb2cif.fetch import get_structure
from old_pdb.pdb_entry import Entry

_LOGGER = logging.getLogger(__name__)


@pytest.mark.parametrize("pdb_id", ["1FAS", "3U7T"], ids=str)
def test_pdb(pdb_id):
    """Test old-format PDB parsing."""
    pdb_fmt = "PDB"
    _LOGGER.info(f"Fetching {pdb_id} in {pdb_fmt} format.")
    pdb_file = get_structure(pdb_id, pdb_fmt)
    entry = Entry()
    entry.parse_file(pdb_file)
    entry.check_master()
