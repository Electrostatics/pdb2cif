"""Classes for PDB records that provide primary structure information.

.. codeauthor::  Todd Dolinsky
.. codeauthor::  Yong Huang
.. codeauthor::  Nathan Baker
"""
from itertools import chain
import logging
from collections import OrderedDict
from .general import BaseRecord, grouper


_LOGGER = logging.getLogger(__name__)


class DatabaseReference(BaseRecord):
    """DBREF record.

    The DBREF record provides cross-reference links between PDB sequences
    (what appears in SEQRES record) and a corresponding database sequence.

    +---------+-------------+-------------+-----------------------------------+
    | COLUMNS | DATA TYPE   | FIELD       | DEFINITION                        |
    +=========+=============+=============+===================================+
    | 1-6     | Record name | "DBREF "    |                                   |
    +---------+-------------+-------------+-----------------------------------+
    | 8-11    | IDcode      | idCode      | ID code of this entry.            |
    +---------+-------------+-------------+-----------------------------------+
    | 13      | Character   | chainID     | Chain identifier.                 |
    +---------+-------------+-------------+-----------------------------------+
    | 15-18   | Integer     | seqBegin    | Initial sequence number of the    |
    |         |             |             | PDB sequence segment.             |
    +---------+-------------+-------------+-----------------------------------+
    | 19      | AChar       | insertBegin | Initial insertion code of the     |
    |         |             |             | PDB sequence segment.             |
    +---------+-------------+-------------+-----------------------------------+
    | 21-24   | Integer     | seqEnd      | Ending sequence number of the PDB |
    |         |             |             | sequence segment.                 |
    +---------+-------------+-------------+-----------------------------------+
    | 25      | AChar       | insertEnd   | Ending insertion code of the PDB  |
    |         |             |             | sequence segment.                 |
    +---------+-------------+-------------+-----------------------------------+
    | 27-32   | LString     | database    | Sequence database name.           |
    +---------+-------------+-------------+-----------------------------------+
    | 34-41   | LString     | dbAccession | Sequence database accession code. |
    +---------+-------------+-------------+-----------------------------------+
    | 43-54   | LString     | dbIdCode    | Sequence database identification  |
    |         |             |             | code.                             |
    +---------+-------------+-------------+-----------------------------------+
    | 56-60   | Integer     | dbseqBegin  | Initial sequence number of the    |
    |         |             |             | database seqment.                 |
    +---------+-------------+-------------+-----------------------------------+
    | 61      | AChar       | idbnsBeg    | Insertion code of initial residue |
    |         |             |             | the segment, if PDB is the        |
    |         |             |             | reference.                        |
    +---------+-------------+-------------+-----------------------------------+
    | 63-67   | Integer     | dbseqEnd    | Ending sequence number of the     |
    |         |             |             | segment.                          |
    +---------+-------------+-------------+-----------------------------------+
    | 68      | AChar       | dbinsEnd    | Insertion code of the ending of   |
    |         |             |             | the segment, if PDB is the        |
    |         |             |             | reference.                        |
    +---------+-------------+-------------+-----------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.id_code = None
        self.chain_id = None
        self.seq_begin = None
        self.ins_begin = None
        self.seq_end = None
        self.ins_end = None
        self.database = None
        self.database_accession = None
        self.database_id_code = None
        self.database_seq_begin = None
        self.database_ins_begin = None
        self.database_seq_end = None
        self.database_ins_end = None

    def parse_line(self, line):
        """Parse DBREF line.

        :param str line:  line to parse
        """
        super().parse_line(line)
        self.id_code = line[7:11].strip()
        self.chain_id = line[12].strip()
        self.seq_begin = int(line[14:18])
        self.ins_begin = line[18].strip()
        self.seq_end = int(line[20:24])
        self.ins_end = line[24].strip()
        self.database = line[26:32].strip()
        self.database_accession = line[33:41].strip()
        self.database_id_code = line[42:54].strip()
        self.database_seq_begin = int(line[55:60])
        self.database_ins_begin = line[60].strip()
        self.database_seq_end = int(line[62:67])
        self.database_ins_end = line[67].strip()

    def __str__(self):
        return (
            f"DBREF  {self.id_code:4} {self.chain_id:1} {self.seq_begin:>4}"
            f"{self.ins_begin:1} {self.seq_end:>4}{self.ins_end:1} "
            f"{self.database:6} {self.database_accession:8} "
            f"{self.database_id_code:12} {self.database_seq_begin:>5}"
            f"{self.database_ins_begin:1} {self.database_seq_end:>5}"
            f"{self.database_ins_end:1}"
        )


class DatabaseReference1(BaseRecord):
    """Provides cross-reference links between PDB sequences (what appears in
    SEQRES record) and a corresponding database sequence.

    This updated two-line format is used when the accession code or sequence
    numbering does not fit the space allotted in the standard DBREF format.
    This includes some GenBank sequence numbering (greater than 5 characters)
    and UNIMES accession numbers (greater than 12 characters).

    +---------+-------------+-------------+-----------------------------------+
    | COLUMNS | DATA TYPE   | FIELD       | DEFINITION                        |
    +=========+=============+=============+===================================+
    | 1-6     | Record name | "DBREF1"    |                                   |
    +---------+-------------+-------------+-----------------------------------+
    | 8-11    | IDcode      | idCode      | ID code of this entry.            |
    +---------+-------------+-------------+-----------------------------------+
    | 13      | Character   | chainID     | Chain identifier.                 |
    +---------+-------------+-------------+-----------------------------------+
    | 15-18   | Integer     | seqBegin    | Initial sequence number of the    |
    |         |             |             | PDB sequence segment, right       |
    |         |             |             | justified.                        |
    +---------+-------------+-------------+-----------------------------------+
    | 19      | AChar       | insertBegin | Initial insertion code of the PDB |
    |         |             |             | sequence segment.                 |
    +---------+-------------+-------------+-----------------------------------+
    | 21-24   | Integer     | seqEnd      | Ending sequence number of the PDB |
    |         |             |             | sequence segment, right           |
    |         |             |             | justified.                        |
    +---------+-------------+-------------+-----------------------------------+
    | 25      | AChar       | insertEnd   | Ending insertion code of the PDB  |
    |         |             |             | sequence segment.                 |
    +---------+-------------+-------------+-----------------------------------+
    | 27-32   | LString     | database    | Sequence database name.           |
    +---------+-------------+-------------+-----------------------------------+
    | 48-67   | LString     | dbIdCode    | Sequence database identification  |
    |         |             |             | code, left justified.             |
    +---------+-------------+-------------+-----------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.ins_code = None
        self.chain_id = None
        self.seq_begin = None
        self.ins_begin = None
        self.seq_end = None
        self.ins_end = None
        self.database = None
        self.db_ins_code = None

    def parse_line(self, line):
        """Parse PDB-format line.

        :param str line:  line with PDB class
        """
        super().parse_line(line)
        self.id_code = line[7:11].strip()
        self.chain_id = line[12].strip()
        self.seq_begin = int(line[14:18])
        self.ins_begin = line[18].strip()
        self.seq_end = int(line[20:24])
        self.ins_end = line[24].strip()
        self.database = line[26:32].strip()
        self.db_ins_code = line[47:67].strip()

    def __str__(self):
        return (
            f"DBREF1 {self.id_code:4} {self.chain_id:1} {self.seq_begin:4}"
            f"{self.ins_begin:1} {self.seq_end:4}{self.ins_end:1}"
            f" {self.database:6}               {self.db_ins_code:20}"
        )


class DatabaseReference2(BaseRecord):
    """Provides cross-reference links between PDB sequences (what appears in
    SEQRES record) and a corresponding database sequence.

    This updated two-line format is used when the accession code or sequence
    numbering does not fit the space allotted in the standard DBREF format.
    This includes some GenBank sequence numbering (greater than 5 characters)
    and UNIMES accession numbers (greater than 12 characters).

    +---------+-------------+-------------+-----------------------------------+
    | COLUMNS | DATA TYPE   | FIELD       | DEFINITION                        |
    +=========+=============+=============+===================================+
    | 1-6     | Record name | "DBREF2"    |                                   |
    +---------+-------------+-------------+-----------------------------------+
    | 8-11    | IDcode      | idCode      | ID code of this entry.            |
    +---------+-------------+-------------+-----------------------------------+
    | 13      | Character   | chainID     | Chain identifier.                 |
    +---------+-------------+-------------+-----------------------------------+
    | 19-40   | LString     | dbAccession | Sequence database accession code, |
    |         |             |             | left justified.                   |
    +---------+-------------+-------------+-----------------------------------+
    | 46-55   | Integer     | seqBegin    | Initial sequence number of the    |
    |         |             |             | Database segment, right           |
    |         |             |             | justified.                        |
    +---------+-------------+-------------+-----------------------------------+
    | 58-67   | Integer     | seqEnd      | Ending sequence number of the     |
    |         |             |             | Database segment, right           |
    |         |             |             | justified.                        |
    +---------+-------------+-------------+-----------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.id_code = None
        self.chain_id = None
        self.db_accession = None
        self.seq_begin = None
        self.seq_seq_end = None

    def parse_line(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_line(line)
        self.id_code = line[7:11].strip()
        self.chain_id = line[12].strip()
        self.db_accession = line[18:40].strip()
        self.seq_begin = int(line[45:55])
        self.seq_end = int(line[57:67])

    def __str__(self):
        return (
            f"DBREF2 {self.id_code:4} {self.chain_id:1}"
            f"     {self.db_accession:22}     {self.seq_begin:10}"
            f"  {self.seq_end:10}"
        )


class ModifiedResidue(BaseRecord):
    """MODRES field

    The MODRES record provides descriptions of modifications (e.g.,
    chemical or post-translational) to protein and nucleic acid residues.
    Included are a mapping between residue names given in a PDB entry and
    standard residues.

    +---------+--------------+----------+-------------------------------------+
    | COLUMNS | DATA TYPE    | FIELD    | DEFINITION                          |
    +=========+==============+==========+=====================================+
    | 1-6     | Record name  | "MODRES" |                                     |
    +---------+--------------+----------+-------------------------------------+
    | 8-11    | IDcode       | idCode   | ID code of this entry.              |
    +---------+--------------+----------+-------------------------------------+
    | 13-15   | Residue name | resName  | Residue name used in this entry.    |
    +---------+--------------+----------+-------------------------------------+
    | 17      | Character    | chainID  | Chain identifier.                   |
    +---------+--------------+----------+-------------------------------------+
    | 19-22   | Integer      | seqNum   | Sequence number.                    |
    +---------+--------------+----------+-------------------------------------+
    | 23      | AChar        | iCode    | Insertion code.                     |
    +---------+--------------+----------+-------------------------------------+
    | 25-27   | Residue name | stdRes   | Standard residue name.              |
    +---------+--------------+----------+-------------------------------------+
    | 30-70   | String       | comment  | Description of the residue          |
    |         |              |          | modification.                       |
    +---------+--------------+----------+-------------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.id_code = None
        self.residue_name = None
        self.chain_id = None
        self.sequence_num = None
        self.ins_code = None
        self.standard_res = None
        self.comment = None

    def parse_line(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_line(line)
        self.id_code = line[7:11].strip()
        self.residue_name = line[12:15].strip()
        self.chain_id = line[16].strip()
        self.sequence_num = int(line[18:22].strip())
        self.ins_code = line[22].strip()
        self.standard_res = line[24:27].strip()
        self.comment = line[29:70].strip()

    def __str__(self):
        return (
            f"MODRES {self.id_code:4} {self.residue_name:3} {self.chain_id:1}"
            f" {self.sequence_num:4}{self.ins_code:1} "
            f"{self.standard_res:3}  {self.comment:41}"
        )


class SequenceDifferences(BaseRecord):
    """SEQADV field

    The SEQADV record identifies conflicts between sequence information in
    the ATOM records of the PDB entry and the sequence database entry given
    on DBREF. Please note that these records were designed to identify
    differences and not errors. No assumption is made as to which database
    contains the correct data. PDB may include REMARK records in the entry
    that reflect the depositor's view of which database has the correct
    sequence.

    +---------+--------------+-------------+----------------------------------+
    | COLUMNS | DATA TYPE    | FIELD       | DEFINITION                       |
    +=========+==============+=============+==================================+
    | 1-6     | Record name  | "SEQADV"    |                                  |
    +---------+--------------+-------------+----------------------------------+
    | 8-11    | IDcode       | idCode      | ID code of this entry.           |
    +---------+--------------+-------------+----------------------------------+
    | 13-15   | Residue name | resName     | Name of the PDB residue in       |
    |         |              |             | conflict.                        |
    +---------+--------------+-------------+----------------------------------+
    | 17      | Character    | chainID     | PDB chain identifier.            |
    +---------+--------------+-------------+----------------------------------+
    | 19-22   | Integer      | seqNum      | PDB sequence number.             |
    +---------+--------------+-------------+----------------------------------+
    | 23      | AChar        | iCode       | PDB insertion code.              |
    +---------+--------------+-------------+----------------------------------+
    | 25-28   | LString      | database    |                                  |
    +---------+--------------+-------------+----------------------------------+
    | 30-38   | LString      | dbAccession | Sequence database accession      |
    |         |              |             | number.                          |
    +---------+--------------+-------------+----------------------------------+
    | 40-42   | Residue name | dbRes       | Sequence database residue name.  |
    +---------+--------------+-------------+----------------------------------+
    | 44-48   | Integer      | dbSeq       | Sequence database sequence       |
    |         |              |             | number.                          |
    +---------+--------------+-------------+----------------------------------+
    | 50-70   | LString      | conflict    | Conflict comment.                |
    +---------+--------------+-------------+----------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.id_code = None
        self.res_name = None
        self.chain_id = None
        self.seq_num = None
        self.ins_code = None
        self.database = None
        self.db_id_code = None
        self.db_res = None
        self.db_seq = ""
        self.conflict = None

    def parse_line(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_line(line)
        self.id_code = line[7:11].strip()
        self.res_name = line[12:15].strip()
        self.chain_id = line[16].strip()
        try:
            self.seq_num = int(line[18:22].strip())
        except ValueError:
            pass
        self.ins_code = line[22].strip()
        self.database = line[24:28].strip()
        self.db_id_code = line[29:38].strip()
        self.db_res = line[39:42].strip()
        try:
            self.db_seq = int(line[43:48].strip())
        except ValueError:
            pass
        self.conflict = line[49:70].strip()

    def __str__(self):
        return (
            f"SEQADV {self.id_code:5}{self.res_name:3} {self.chain_id:1}"
            f" {self.seq_num:4}{self.ins_code:1} {self.database:4}"
            f" {self.db_id_code:9} {self.db_res:3} {self.db_seq:5}"
            f" {self.conflict:21}"
        )


class SequenceResidues(BaseRecord):
    """SEQRES field

    SEQRES records contain the amino acid or nucleic acid sequence of
    residues in each chain of the macromolecule that was studied.

    +---------+--------------+----------+-------------------------------------+
    | COLUMNS | DATA TYPE    | FIELD    | DEFINITION                          |
    +=========+==============+==========+=====================================+
    | 1-6     | Record name  | "SEQRES" |                                     |
    +---------+--------------+----------+-------------------------------------+
    | 8-10    | Integer      | serNum   | Serial number of the SEQRES record  |
    |         |              |          | for the current chain. Starts at 1  |
    |         |              |          | and increments by one each line.    |
    |         |              |          | Reset to 1 for each chain.          |
    +---------+--------------+----------+-------------------------------------+
    | 12      | Character    | chainID  | Chain identifier. This may be any   |
    |         |              |          | single legal character, including a |
    |         |              |          | blank which is is used if there is  |
    |         |              |          | only one chain.                     |
    +---------+--------------+----------+-------------------------------------+
    | 14-17   | Integer      | numRes   | Number of residues in the chain.    |
    |         |              |          | This value is repeated on every     |
    |         |              |          | record.                             |
    +---------+--------------+----------+-------------------------------------+
    | 20-22   | Residue name | resName  | Residue name.                       |
    +---------+--------------+----------+-------------------------------------+
    | 24-26   | Residue name | resName  | Residue name.                       |
    +---------+--------------+----------+-------------------------------------+
    | 28-30   | Residue name | resName  | Residue name.                       |
    +---------+--------------+----------+-------------------------------------+
    | 32-34   | Residue name | resName  | Residue name.                       |
    +---------+--------------+----------+-------------------------------------+
    | 36-38   | Residue name | resName  | Residue name.                       |
    +---------+--------------+----------+-------------------------------------+
    | 40-42   | Residue name | resName  | Residue name.                       |
    +---------+--------------+----------+-------------------------------------+
    | 44-46   | Residue name | resName  | Residue name.                       |
    +---------+--------------+----------+-------------------------------------+
    | 48-50   | Residue name | resName  | Residue name.                       |
    +---------+--------------+----------+-------------------------------------+
    | 52-54   | Residue name | resName  | Residue name.                       |
    +---------+--------------+----------+-------------------------------------+
    | 56-58   | Residue name | resName  | Residue name.                       |
    +---------+--------------+----------+-------------------------------------+
    | 60-62   | Residue name | resName  | Residue name.                       |
    +---------+--------------+----------+-------------------------------------+
    | 64-66   | Residue name | resName  | Residue name.                       |
    +---------+--------------+----------+-------------------------------------+
    | 68-70   | Residue name | resName  | Residue name.                       |
    +---------+--------------+----------+-------------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.residues = OrderedDict()
        self.num_residues = {}

    def parse_line(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_line(line)
        chain_id = line[11].strip()
        if chain_id not in self.residues:
            self.residues[chain_id] = []
        self.num_residues [chain_id]= int(line[13:17])
        start = 19
        end = 22
        while True:
            self.residues[chain_id].append(line[start:end].strip())
            start += 4
            end += 4
            if start > 68:
                break

    def num_chains(self) -> int:
        """Number of chains in sequence."""
        return len(self.residues)

    def __str__(self):
        strings = []
        for chain_id, residues in self.residues.items():
            for ichunk, chunk in enumerate(grouper(residues, 13)):
                serial_num = ichunk + 1
                string = f"SEQRES {serial_num:>3} {chain_id:1} "
                string += f"{self.num_residues[chain_id]:>4} "
                for residue in chunk:
                    string += f" {residue:>3}"
                strings.append(string.strip())
        return "\n".join(strings)
