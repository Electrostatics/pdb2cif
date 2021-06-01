"""Test old-format PDB output formatting."""
import logging
import pytest
import difflib
from pprint import pprint
from common import get_cif, TEST_DATA, REF_LINE
from pdb2cif.pdb_entry import Entry


_LOGGER = logging.getLogger(__name__)


@pytest.mark.parametrize("pdb_id", TEST_DATA.index, ids=str)
def test_cif(pdb_id):
    """Test old-format PDB parsing."""
    _LOGGER.debug(f"Parsing {pdb_id} in PDB format.")
    cif_file = get_cif(pdb_id)
    entry = Entry()
    entry.parse_cif_file(cif_file)
    _LOGGER.debug(str(entry))
    if pdb_id not in ["1AFS", "1E7G"]:
        err = f"Haven't checked results for {pdb_id} yet."
        raise NotImplementedError(err)
