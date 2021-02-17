"""Functions used by all tests."""
import logging
from pathlib import Path
import pandas as pd
from fetch import get_structure


_LOGGER = logging.getLogger(__name__)
REF_LINE = (
    "0        1         2         3         4         5         6         7         8\n"
    "12345678901234567890123456789012345678901234567890123456789012345678901234567890"
)
DATA_DIR = Path("tests/data")
TEST_CSV_PATH = DATA_DIR / Path("tests.csv")
TEST_DATA = pd.read_csv(TEST_CSV_PATH, index_col="PDB_ID")


def get_pdb(pdb_id):
    """Get PDB file from local data or remotely.

    :param str pdb_id:  4-character PDB ID
    :returns:  file object ready for reading
    """
    pdb_path = DATA_DIR / f"{pdb_id.lower()}.pdb"
    if not pdb_path.exists():
        _LOGGER.warning(f"Fetching {pdb_id} remotely.")
        return get_structure(pdb_id, "PDB")
    else:
        return open(pdb_path, "rt", encoding="utf8")
