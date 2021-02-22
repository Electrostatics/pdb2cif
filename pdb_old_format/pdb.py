""" PDB parsing class

This module parses PDBs in accordance to PDB Format Description Version 2.2
(1996); it is not very forgiving.  Each class in this module corresponds
to a record in the PDB Format Description.  Much of the documentation for
the classes is taken directly from the above PDB Format Description.

.. codeauthor::  Todd Dolinsky
.. codeauthor::  Yong Huang
.. codeauthor::  Nathan Baker
"""
import logging


_LOGGER = logging.getLogger(__name__)
LINE_PARSERS = {}


def register_line_parser(klass):
    """Register a line parser in the global dictionary.

    :param klass:  class for line parser
    """
    LINE_PARSERS[klass.__name__] = klass
    return klass


@register_line_parser
class END(BaseRecord):
    """END class

    The END records are paired with MODEL records to group individual
    structures found in a coordinate entry.
    """

    def __str__(self):
        return "END"


@register_line_parser
class ENDMDL(BaseRecord):
    """ENDMDL class

    The ENDMDL records are paired with MODEL records to group individual
    structures found in a coordinate entry.
    """

    def __str__(self):
        return "ENDMDL"


@register_line_parser
class SIGUIJ(BaseRecord):
    """SIGUIJ class

    The SIGUIJ records present the anisotropic temperature factors.
    """

    def __init__(self):
        self.serial = None
        self.name = None
        self.alt_loc = None
        self.res_name = None
        self.chain_id = None
        self.res_seq = None
        self.ins_code = None
        self.sig11 = None
        self.sig22 = None
        self.sig33 = None
        self.sig12 = None
        self.sig13 = None
        self.sig23 = None
        self.seg_id = None
        self.element = None
        self.charge = None

    def parse_line(self, line):
        """Initialize by parsing line:

        +---------+--------+----------+-------------------------------------+
        | COLUMNS | TYPE   | FIELD    | DEFINITION                          |
        +=========+========+==========+=====================================+
        | 7-11    | int    | serial   | Atom serial number.                 |
        +---------+--------+----------+-------------------------------------+
        | 13-16   | string | name     | Atom name.                          |
        +---------+--------+----------+-------------------------------------+
        | 17      | string | alt_loc  | Alternate location indicator.       |
        +---------+--------+----------+-------------------------------------+
        | 18-20   | string | res_name | Residue name.                       |
        +---------+--------+----------+-------------------------------------+
        | 22      | string | chain_id | Chain identifier.                   |
        +---------+--------+----------+-------------------------------------+
        | 23-26   | int    | res_seq  | Residue sequence number.            |
        +---------+--------+----------+-------------------------------------+
        | 27      | string | ins_code | Insertion code.                     |
        +---------+--------+----------+-------------------------------------+
        | 29-35   | int    | sig11    | Sigma U(1,1)                        |
        +---------+--------+----------+-------------------------------------+
        | 36-42   | int    | sig22    | Sigma U(2,2)                        |
        +---------+--------+----------+-------------------------------------+
        | 43-49   | int    | sig33    | Sigma U(3,3)                        |
        +---------+--------+----------+-------------------------------------+
        | 50-56   | int    | sig12    | Sigma U(1,2)                        |
        +---------+--------+----------+-------------------------------------+
        | 57-63   | int    | sig13    | Sigma U(1,3)                        |
        +---------+--------+----------+-------------------------------------+
        | 64-70   | int    | sig23    | Sigma U(2,3)                        |
        +---------+--------+----------+-------------------------------------+
        | 73-76   | string | seg_id   | Segment identifier, left-justified. |
        +---------+--------+----------+-------------------------------------+
        | 77-78   | string | el.ment  | Element symbol, right-justified.    |
        +---------+--------+----------+-------------------------------------+
        | 79-80   | string | charge   | Charge on the atom.                 |
        +---------+--------+----------+-------------------------------------+

        :param str line:  line with PDB class
        """
        super().parse_line(line)
        self.serial = int(line[6:11].strip())
        self.name = line[12:16].strip()
        self.alt_loc = line[16].strip()
        self.res_name = line[17:20].strip()
        self.chain_id = line[21].strip()
        self.res_seq = int(line[22:26].strip())
        self.ins_code = line[26].strip()
        self.sig11 = int(line[28:35].strip())
        self.sig22 = int(line[35:42].strip())
        self.sig33 = int(line[42:49].strip())
        self.sig12 = int(line[49:56].strip())
        self.sig13 = int(line[56:63].strip())
        self.sig23 = int(line[63:70].strip())
        self.seg_id = line[72:76].strip()
        self.element = line[76:78].strip()
        self.charge = line[78:80].strip()

    def __str__(self):
        raise NotImplementedError()


@register_line_parser
class SIGATM(BaseRecord):
    """SIGATM class

    The SIGATM records present the standard deviation of atomic parameters
    as they appear in ATOM and HETATM records.
    """

    def __init__(self):
        self.serial = None
        self.name = None
        self.alt_loc = None
        self.res_name = None
        self.chain_id = None
        self.res_seq = None
        self.ins_code = None
        self.sig_x = None
        self.sig_y = None
        self.sig_z = None
        self.sig_occ = None
        self.sig_temp = None
        self.seg_id = None
        self.element = None
        self.charge = None

    def parse_line(self, line):
        """Initialize by parsing line

        +---------+--------+----------+-------------------------------------+
        | COLUMNS | TYPE   | FIELD    | DEFINITION                          |
        +=========+========+==========+=====================================+
        | 7-11    | int    | serial   | Atom serial number.                 |
        +---------+--------+----------+-------------------------------------+
        | 13-16   | string | name     | Atom name.                          |
        +---------+--------+----------+-------------------------------------+
        | 17      | string | alt_loc  | Alternate location indicator.       |
        +---------+--------+----------+-------------------------------------+
        | 18-20   | string | res_name | Residue name.                       |
        +---------+--------+----------+-------------------------------------+
        | 22      | string | chain_id | Chain identifier.                   |
        +---------+--------+----------+-------------------------------------+
        | 23-26   | int    | res_seq  | Residue sequence number.            |
        +---------+--------+----------+-------------------------------------+
        | 27      | string | ins_code | Code for insertion of residues.     |
        +---------+--------+----------+-------------------------------------+
        | 31-38   | float  | sig_x    | Standard deviation of orthogonal    |
        |         |        |          | coordinates for X in Angstroms.     |
        +---------+--------+----------+-------------------------------------+
        | 39-46   | float  | sig_y    | Standard deviation of orthogonal    |
        |         |        |          | coordinates for Y in Angstroms.     |
        +---------+--------+----------+-------------------------------------+
        | 47-54   | float  | sig_z    | Standard deviation of orthogonal    |
        |         |        |          | coordinates for Z in Angstroms.     |
        +---------+--------+----------+-------------------------------------+
        | 55-60   | float  | sig_occ  | Standard deviation of occupancy.    |
        +---------+--------+----------+-------------------------------------+
        | 61-66   | float  | sig_temp | Standard deviation of temperature   |
        |         |        |          | factor.                             |
        +---------+--------+----------+-------------------------------------+
        | 73-76   | string | seg_id   | Segment identifier, left-justified. |
        +---------+--------+----------+-------------------------------------+
        | 77-78   | string | element  | Element symbol, right-justified.    |
        +---------+--------+----------+-------------------------------------+
        | 79-80   | string | charge   | Charge on the atom.                 |
        +---------+--------+----------+-------------------------------------+

        :param str line:  line with PDB class
        """
        super().parse_line(line)
        self.serial = int(line[6:11].strip())
        self.name = line[12:16].strip()
        self.alt_loc = line[16].strip()
        self.res_name = line[17:20].strip()
        self.chain_id = line[21].strip()
        self.res_seq = int(line[22:26].strip())
        self.ins_code = line[26].strip()
        self.sig_x = float(line[30:38].strip())
        self.sig_y = float(line[38:46].strip())
        self.sig_z = float(line[46:54].strip())
        self.sig_occ = float(line[54:60].strip())
        self.sig_temp = float(line[60:66].strip())
        self.seg_id = line[72:76].strip()
        self.element = line[76:78].strip()
        self.charge = line[78:80].strip()

    def __str__(self):
        raise NotImplementedError()


@register_line_parser
class TVECT(BaseRecord):
    """TVECT class

    The TVECT records present the translation vector for infinite
    covalently connected structures.
    """

    def __init__(self):
        self.serial = None
        self.trans1 = None
        self.trans2 = None
        self.trans3 = None
        self.text = None

    def parse_line(self, line):
        """Initialize by parsing line

        +---------+--------+--------+----------------------------------+
        | COLUMNS | TYPE   | FIELD  | DEFINITION                       |
        +=========+========+========+==================================+
        | 8-10    | int    | serial | Serial number                    |
        +---------+--------+--------+----------------------------------+
        | 11-20   | float  | t1     | Components of translation vector |
        +---------+--------+--------+----------------------------------+
        | 21-30   | float  | t2     | Components of translation vector |
        +---------+--------+--------+----------------------------------+
        | 31-40   | float  | t2     | Components of translation vector |
        +---------+--------+--------+----------------------------------+
        | 41-70   | string | text   | Comments                         |
        +---------+--------+--------+----------------------------------+

        :param str line:  line with PDB class
        """
        super().parse_line(line)
        self.serial = int(line[7:10].strip())
        self.trans1 = float(line[10:20].strip())
        self.trans2 = float(line[20:30].strip())
        self.trans3 = float(line[30:40].strip())
        self.text = line[40:70].strip()


@register_line_parser
class SLTBRG(BaseRecord):
    """SLTBRG field

    The SLTBRG records specify salt bridges in the entry.
    records and is provided here for convenience in searching.
    """

    def __init__(self):
        self.name1 = None
        self.alt_loc1 = None
        self.res_name1 = None
        self.chain_id1 = None
        self.res_seq1 = None
        self.ins_code1 = None
        self.name2 = None
        self.alt_loc2 = None
        self.res_name2 = None
        self.chain_id2 = None
        self.res_seq2 = None
        self.ins_code2 = None
        self.sym1 = None
        self.sym2 = None

    def parse_line(self, line):
        """Initialize by parsing line

        +---------+--------+-----------+---------------------------------+
        | COLUMNS | TYPE   | FIELD     | DEFINITION                      |
        +=========+========+===========+=================================+
        | 13-16   | string | name1     | Atom name.                      |
        +---------+--------+-----------+---------------------------------+
        | 17      | string | alt_loc1  | Alternate location indicator.   |
        +---------+--------+-----------+---------------------------------+
        | 18-20   | string | res_name1 | Residue name.                   |
        +---------+--------+-----------+---------------------------------+
        | 22      | string | chain_id1 | Chain identifier.               |
        +---------+--------+-----------+---------------------------------+
        | 23-26   | int    | res_seq1  | Residue sequence number.        |
        +---------+--------+-----------+---------------------------------+
        | 27      | string | ins_code1 | Insertion code.                 |
        +---------+--------+-----------+---------------------------------+
        | 43-46   | string | name2     | Atom name.                      |
        +---------+--------+-----------+---------------------------------+
        | 47      | string | alt_loc2  | Alternate location indicator.   |
        +---------+--------+-----------+---------------------------------+
        | 48-50   | string | res_name2 | Residue name.                   |
        +---------+--------+-----------+---------------------------------+
        | 52      | string | chain_id2 | Chain identifier.               |
        +---------+--------+-----------+---------------------------------+
        | 53-56   | int    | res_seq2  | Residue sequence number.        |
        +---------+--------+-----------+---------------------------------+
        | 57      | string | ins_code2 | Insertion code.                 |
        +---------+--------+-----------+---------------------------------+
        | 60-65   | string | sym1      | Symmetry operator for 1st atom. |
        +---------+--------+-----------+---------------------------------+
        | 67-72   | string | sym2      | Symmetry operator for 2nd atom. |
        +---------+--------+-----------+---------------------------------+

        :param str line:  line with PDB class
        """
        super().parse_line(line)
        self.name1 = line[12:16].strip()
        self.alt_loc1 = line[16].strip()
        self.res_name1 = line[17:20].strip()
        self.chain_id1 = line[21].strip()
        self.res_seq1 = int(line[22:26].strip())
        self.ins_code1 = line[26].strip()
        self.name2 = line[42:46].strip()
        self.alt_loc2 = line[46].strip()
        self.res_name2 = line[47:50].strip()
        self.chain_id2 = line[51].strip()
        self.res_seq2 = int(line[52:56].strip())
        self.ins_code2 = line[56].strip()
        self.sym1 = line[59:65].strip()
        self.sym2 = line[66:72].strip()


@register_line_parser
class HYDBND(BaseRecord):
    """HYDBND field

    The HYDBND records specify hydrogen bonds in the entry.
    """

    def __init__(self):
        self.name1 = None
        self.alt_loc1 = None
        self.res_name1 = None
        self.chain1 = None
        self.res_seq1 = None
        self.i_code1 = None
        self.name_h = None
        self.alt_loc_h = None
        self.chain_h = None
        self.res_seq_h = None
        self.i_code_h = None
        self.name2 = None
        self.alt_loc2 = None
        self.res_name2 = None
        self.chain2 = None
        self.res_seq2 = None
        self.i_code2 = None
        self.sym1 = None
        self.sym2 = None

    def parse_line(self, line):
        """Initialize by parsing line

        +---------+--------+-----------+-------------------------------------+
        | COLUMNS | TYPE   | FIELD     | DEFINITION                          |
        +=========+========+===========+=====================================+
        | 13-16   | string | name1     | Atom name.                          |
        +---------+--------+-----------+-------------------------------------+
        | 17      | string | alt_loc1  | Alternate location indicator.       |
        +---------+--------+-----------+-------------------------------------+
        | 18-20   | string | res_name1 | Residue name.                       |
        +---------+--------+-----------+-------------------------------------+
        | 22      | string | chain1    | Chain identifier.                   |
        +---------+--------+-----------+-------------------------------------+
        | 23-27   | int    | res_seq1  | Residue sequence number.            |
        +---------+--------+-----------+-------------------------------------+
        | 28      | string | i_code1   | Insertion code.                     |
        +---------+--------+-----------+-------------------------------------+
        | 30-33   | string | name_h    | Hydrogen atom name.                 |
        +---------+--------+-----------+-------------------------------------+
        | 34      | string | alt_loc_h | Alternate location indicator.       |
        +---------+--------+-----------+-------------------------------------+
        | 36      | string | chain_h   | Chain identifier.                   |
        +---------+--------+-----------+-------------------------------------+
        | 37-41   | int    | res_seq_h | Residue sequence number.            |
        +---------+--------+-----------+-------------------------------------+
        | 42      | string | ins_codeH | Insertion code.                     |
        +---------+--------+-----------+-------------------------------------+
        | 44-47   | string | name2     | Atom name.                          |
        +---------+--------+-----------+-------------------------------------+
        | 48      | string | alt_loc2  | Alternate location indicator.       |
        +---------+--------+-----------+-------------------------------------+
        | 49-51   | string | res_name2 | Residue name.                       |
        +---------+--------+-----------+-------------------------------------+
        | 53      | string | chain_id2 | Chain identifier.                   |
        +---------+--------+-----------+-------------------------------------+
        | 54-58   | int    | res_seq2  | Residue sequence number.            |
        +---------+--------+-----------+-------------------------------------+
        | 59      | string | ins_code2 | Insertion code.                     |
        +---------+--------+-----------+-------------------------------------+
        | 60-65   | string | sym1      | Symmetry operator for 1st           |
        |         |        |           | non-hydrogen atom.                  |
        +---------+--------+-----------+-------------------------------------+
        | 67-72   | string | sym2      | Symmetry operator for 2nd           |
        |         |        |           | non-hydrogen atom.                  |
        +---------+--------+-----------+-------------------------------------+

        :param str line:  line with PDB class
        """
        super().parse_line(line)
        self.name1 = line[12:16].strip()
        self.alt_loc1 = line[16].strip()
        self.res_name1 = line[17:20].strip()
        self.chain1 = line[21].strip()
        self.res_seq1 = line[22:27].strip()
        self.i_code1 = line[27].strip()
        self.name_h = line[29:33].strip()
        self.alt_loc_h = line[33].strip()
        self.chain_h = line[35].strip()
        self.res_seq_h = line[36:41].strip()
        self.i_code_h = line[41].strip()
        self.name2 = line[43:47].strip()
        self.alt_loc2 = line[47].strip()
        self.res_name2 = line[48:51].strip()
        self.chain2 = line[52].strip()
        self.res_seq2 = line[53:58].strip()
        self.i_code2 = line[58].strip()
        self.sym1 = line[59:65].strip()
        self.sym2 = line[66:72].strip()


@register_line_parser
class TURN(BaseRecord):
    """TURN field

    The TURN records identify turns and other short loop turns which normally
    connect other secondary structure segments.
    """

    def __init__(self):
        self.seq = None
        self.turn_id = None
        self.init_res_name = None
        self.init_chain_id = None
        self.init_seq_num = None
        self.init_i_code = None
        self.end_res_name = None
        self.end_chain_id = None
        self.end_seq_num = None
        self.end_i_code = None
        self.comment = None

    def parse_line(self, line):
        """Initialize by parsing line

        +---------+--------+---------------+---------------------------------+
        | COLUMNS | TYPE   | FIELD         | DEFINITION                      |
        +=========+========+===============+=================================+
        | 8-10    | int    | seq           | Turn number; starts with 1 and  |
        |         |        |               | increments by one.              |
        +---------+--------+---------------+---------------------------------+
        | 12-14   | string | turn_id       | Turn identifier.                |
        +---------+--------+---------------+---------------------------------+
        | 16-18   | string | init_res_name | Residue name of initial residue |
        |         |        |               | in turn.                        |
        +---------+--------+---------------+---------------------------------+
        | 20      | string | init_chain_id | Chain identifier for the chain  |
        |         |        |               | containing this turn.           |
        +---------+--------+---------------+---------------------------------+
        | 21-24   | int    | init_seq_num  | Sequence number of initial      |
        |         |        |               | residue in turn.                |
        +---------+--------+---------------+---------------------------------+
        | 25      | string | init_i_code   | Insertion code of initial       |
        |         |        |               | residue in turn.                |
        +---------+--------+---------------+---------------------------------+
        | 27-29   | string | end_res_name  | Residue name of terminal residue|
        |         |        |               | of turn.                        |
        +---------+--------+---------------+---------------------------------+
        | 31      | string | end_chain_id  | Chain identifier for the chain  |
        |         |        |               | containing this turn.           |
        +---------+--------+---------------+---------------------------------+
        | 32-35   | int    | end_seq_num   | Sequence number of terminal     |
        |         |        |               | residue of turn.                |
        +---------+--------+---------------+---------------------------------+
        | 36      | string | end_i_code    | Insertion code of terminal      |
        |         |        |               | residue of turn.                |
        +---------+--------+---------------+---------------------------------+
        | 41-70   | string | comment       | Associated comment.             |
        +---------+--------+---------------+---------------------------------+

        :param str line:  line with PDB class
        """
        super().parse_line(line)
        self.seq = int(line[7:10].strip())
        self.turn_id = line[11:14].strip()
        self.init_res_name = line[15:18].strip()
        self.init_chain_id = line[19].strip()
        self.init_seq_num = int(line[20:24].strip())
        self.init_i_code = line[24].strip()
        self.end_res_name = line[26:29].strip()
        self.end_chain_id = line[30].strip()
        self.end_seq_num = int(line[31:35].strip())
        self.end_i_code = line[35].strip()
        self.comment = line[40:70].strip()


def read_file(file_, tolerant=True) -> tuple:
    """Parse PDB-format data into array of Atom objects.

    :param file file_:  open File-like object
    :param bool tolerant:  whether to supress or raise parsing errors
    :rtype:  (list, list)
    """
    pdblist = []  # Array of parsed lines (as objects)
    errlist = []  # List of records we can't parse

    # We can come up with nothing if can't get our file off the web.
    if file_ is None:
        return pdblist, errlist

    while True:
        line = file_.readline().strip()
        if line == "":
            break

        # We assume we have a method for each PDB record and can therefore
        # parse them automatically
        record = ""
        try:
            record = line[0:6].strip()
            if record not in errlist:
                klass = LINE_PARSERS[record]
                obj = klass()
                obj.parse_line(line)
                pdblist.append(obj)
        except (KeyError, ValueError) as details:
            if (record not in ["HETATM", "ATOM"]) and tolerant:
                errlist.append(record)
                _LOGGER.error(f"Error parsing line: {details}")
                _LOGGER.error(f"<{line.strip()}>")
                _LOGGER.error(
                    f"Truncating remaining errors for record type:{record}"
                )
            else:
                raise details
    return pdblist, errlist
