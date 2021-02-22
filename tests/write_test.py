"""Test old-format PDB output formatting."""
import logging
import pytest
import difflib
from pprint import pprint
from common import get_pdb, TEST_DATA, REF_LINE
from old_pdb.pdb_entry import Entry


_LOGGER = logging.getLogger(__name__)


@pytest.mark.parametrize("pdb_id", TEST_DATA.index, ids=str)
def test_pdb(pdb_id):
    """Test old-format PDB parsing."""
    _LOGGER.debug(f"Parsing {pdb_id} in PDB format.")
    pdb_file = get_pdb(pdb_id)
    entry = Entry()
    entry.parse_file(pdb_file)

    new_strings = [line.strip() for line in str(entry).splitlines()]
    old_strings = [line.strip() for line in get_pdb(pdb_id).readlines()]
    for iline in range(len(new_strings)):
        old_line = old_strings[iline]
        new_line = new_strings[iline]
        diff = list(
            difflib.unified_diff([old_line], [new_line])
        )
        if len(list(diff)) > 0:
            diff_string = "\n".join(diff)
            old_lines = "\n".join(old_strings[iline-1:iline+2])
            new_lines = "\n".join(new_strings[iline-1:iline+2])
            err = "\n".join(
                [
                    f"Difference in text output:\n{diff_string}",
                    f"Offending line:\n{REF_LINE}\n{new_lines}",
                    f"Line should be:\n{REF_LINE}\n{old_lines}"
                ]
            )
            _LOGGER.error(err)
            raise AssertionError()
