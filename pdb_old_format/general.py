"""Functions and classes for all PDB record types.

.. codeauthor::  Todd Dolinsky
.. codeauthor::  Yong Huang
.. codeauthor::  Nathan Baker
"""
import logging
from datetime import datetime
from itertools import zip_longest


_LOGGER = logging.getLogger(__name__)
DATE_FMT = r"%d-%b-%y"


def grouper(iterable, block_size, fillvalue=None) -> list:
    """Group an iterable into chunks of block_size.

    Adapted from https://docs.python.org/3/library/itertools.html#itertools-recipes.

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


class BaseRecord:
    """Base class for all PDB records."""

    def __init__(self):
        self.original_text = []

    def parse_line(self, line):
        """Parse line of PDB file.

        :param str line:  PDB file line to parse
        """
        if line is not None:
            self.original_text.append(line.rstrip("\r\n"))

    def __str__(self):
        raise NotImplementedError("BaseRecord does not implement __str__.")
