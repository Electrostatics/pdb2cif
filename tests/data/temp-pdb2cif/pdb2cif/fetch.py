"""Functions for getting PDB files from the archive."""
import logging
import io
import requests


_LOGGER = logging.getLogger(__name__)
PDB_URL_FMT = "http://files.rcsb.org/download/{pdb_id}.pdb"
CIF_URL_FMT = "http://files.rcsb.org/download/{pdb_id}.cif"


def get_structure(pdb_id, pdb_fmt) -> io.StringIO:
    """Obtain a structure in PDB format.

    Fetch the file from the PDB webserver at http://www.rcsb.org/pdb/

    :param str pdb_id:  PDB ID
    :param str pdb_fmt:  format of structure (CIF or PDB)
    :returns:  file-like object containing PDB file
    """
    if pdb_fmt.upper() == "CIF":
        url_path = CIF_URL_FMT.format(pdb_id=pdb_id)
    elif pdb_fmt.upper() == "PDB":
        url_path = PDB_URL_FMT.format(pdb_id=pdb_id)
    _LOGGER.debug(f"Attempting to fetch PDB from {url_path}")
    resp = requests.get(url_path)
    if resp.status_code != 200:
        errstr = f"Got code {resp.status_code} while retrieving {url_path}"
        raise IOError(errstr)
    return io.StringIO(resp.text)
