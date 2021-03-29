"""Functions and classes for all PDB record types.

.. codeauthor::  Todd Dolinsky
.. codeauthor::  Yong Huang
.. codeauthor::  Nathan Baker
"""
import logging
from datetime import datetime
from itertools import zip_longest
from pandas import DataFrame


_LOGGER = logging.getLogger(__name__)
DATE_FMT = r"%d-%b-%y"


def grouper(iterable, block_size, fillvalue=None) -> list:
    """Group an iterable into chunks of block_size.

    Adapted from
    https://docs.python.org/3/library/itertools.html#itertools-recipes.

    :param list iterable:  list to break into chunks
    :param int block_size:  chunk size
    :param fillvalue:  fill in chunks without values
    :returns:  list of chunks
    """
    args = [iter(iterable)] * block_size
    return zip_longest(*args, fillvalue=fillvalue)


def date_parse(date_string) -> datetime:
    """Parse date string in PDB format into date object.

    :param str date_string:  date string in PDB format
    :returns:  datetime object
    """
    return datetime.strptime(date_string, DATE_FMT)


def date_format(date) -> str:
    """Return date formatted according to PDB guidelines.

    :param datetime date:  input date
    :returns:  formatted date string
    """
    if date is not None:
        str_ = date.strftime(DATE_FMT).upper()
    else:
        str_ = "         "
    return str_


def atom_format(record) -> str:
    """Atom name logic is very complicated.

    https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html

    :param BaseRecord record:  atom record
    :returns formatted atom name:
    """
    if len(record.name) == 4:
        return record.name
    if record.name[:2] == record.element[:2]:
        return f"{record.name:>2}  "[:4]
    return f" {record.name:<3}"[:4]


def cif_df(cif_object) -> DataFrame:
    """Convert a CIF object to a DataFrame.

    :param :class:`pdbx.containers.DataCategory` cif_object:  object to convert
    :returns:  DataFrame with CIF object data
    """
    if cif_object is None:
        return DataFrame()
    row_list = cif_object.row_list
    attr_list = cif_object.attribute_list
    return DataFrame(data=row_list, columns=attr_list)


class BaseRecord:
    """Base class for all PDB records."""

    def __init__(self):
        self.original_text = []

    def parse_pdb(self, line):
        """Parse line of PDB file.

        :param str line:  PDB file line to parse
        """
        if line is not None:
            self.original_text.append(line.rstrip("\r\n"))

    def parse_cif(self, container) -> bool:
        """Parse CIF container for information about this record.

        :param :class:`pdbx.containers.DataContainer` container:  container to
            parse
        :returns:  True if useful information was extracted from container
        """
        raise NotImplementedError("BaseRecord does not implement parse_cif.")

    def to_cif(self, container):
        """Add information about this record to CIF container.

        :param :class:`pdbx.containers.DataContainer` container:  container to
            augment
        """
        raise NotImplementedError("BaseRecord does not implement to_cif.")

    def __str__(self):
        raise NotImplementedError("BaseRecord does not implement __str__.")
