"""Test old-format PDB parsing."""
import logging
import pytest
from common import get_pdb, TEST_DATA
from old_pdb.pdb_entry import Entry


_LOGGER = logging.getLogger(__name__)


@pytest.mark.parametrize("pdb_id", TEST_DATA.index, ids=str)
def test_pdb(pdb_id):
    """Test old-format PDB parsing."""
    _LOGGER.debug(f"Parsing {pdb_id} in PDB format.")
    pdb_file = get_pdb(pdb_id)
    entry = Entry()
    entry.parse_file(pdb_file)
    entry.check_master()

    test_row = TEST_DATA.loc[pdb_id, :]
    assert entry.num_atoms() == test_row["Num_heavy_atoms"]
    assert entry.num_residues() == test_row["Num_residues"]
    assert entry.num_chains() == test_row["Num_chains"]
    num_model = test_row["Num_models"]
    if num_model > 1:
        assert entry.num_model.model_number == test_row["Num_models"]
    assert len(entry.model) == test_row["Num_models"]
