"""Functions used by all tests."""
import logging
from pathlib import Path
from fetch import get_structure


_LOGGER = logging.getLogger(__name__)
REF_LINE = (
    "0        1         2         3         4         5         6         7         8\n"
    "12345678901234567890123456789012345678901234567890123456789012345678901234567890"
)
PDB_DIR = Path("tests/data")
TEST_LIST = {
    "1FAS",
    "3U7T",
    "2HSY",
    "2HXH",
    "6PU4",
    "7KLH",
    "7KMK",
    "7JK8",
    "6XN9",
    "7LB7",
    "6ZHN",
    "6VC1",
    "3PDM",
    "6VAR",
    "4J7Z",
    "1RMN",
    "1SSZ",
    "1VYC",
    "4UN3", "1K1I", "1AFS", "1FAS", "5DV8", "5D8V", "1E7G"
}


def get_pdb(pdb_id):
    """Get PDB file from local data or remotely.

    :param str pdb_id:  4-character PDB ID
    :returns:  file object ready for reading
    """
    pdb_path = PDB_DIR / f"{pdb_id.lower()}.pdb"
    if not pdb_path.exists():
        _LOGGER.error(f"Fetching {pdb_id} remotely.")
        return get_structure(pdb_id, "PDB")
    else:
        return open(pdb_path, "rt", encoding="utf8")
