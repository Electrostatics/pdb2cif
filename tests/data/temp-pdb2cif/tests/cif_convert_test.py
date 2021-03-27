"""Test conversion between formats."""
import logging
import pytest
import pdbx
from pdb2cif.fetch import get_structure
from pdb2cif.convert_cif import convert

# from pdb2cif.cif2pdb import convert as cif_convert


_LOGGER = logging.getLogger(__name__)


@pytest.mark.parametrize("pdb_id", ["3WXV", "1FAS", "3U7T"], ids=str)
def test_pdb(pdb_id):
    """Test old-format PDB parsing."""
    pdb_fmt = "CIF"
    _LOGGER.info(f"Fetching {pdb_id} in {pdb_fmt} format.")
    cif_file = get_structure(pdb_id, pdb_fmt)
    containers = pdbx.load(cif_file)
    _LOGGER.info(f"Got {len(containers)} CIF containers.")
    entry = convert(containers)
