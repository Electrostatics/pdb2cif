"""Classes for PDB records that provide primary structure information.

.. codeauthor::  Todd Dolinsky
.. codeauthor::  Yong Huang
.. codeauthor::  Nathan Baker
"""
import logging
from collections import OrderedDict
import pandas as pd

from pandas.core.algorithms import diff
from .general import BaseRecord, grouper, cif_df


_LOGGER = logging.getLogger(__name__)


class DatabaseReference(BaseRecord):
    """DBREF record.

    The DBREF record provides cross-reference links between PDB sequences
    (what appears in SEQRES record) and a corresponding database sequence.

    +---------+-------------+--------------------+-----------------------------+
    | COLUMNS | DATA TYPE   | FIELD              | DEFINITION                  |
    +=========+=============+====================+=============================+
    | 1-6     | Record name | "DBREF "           |                             |
    +---------+-------------+--------------------+-----------------------------+
    | 8-11    | IDcode      | id_code            | ID code of this entry.      |
    +---------+-------------+--------------------+-----------------------------+
    | 13      | Character   | chain_id           | Chain identifier.           |
    +---------+-------------+--------------------+-----------------------------+
    | 15-18   | Integer     | seq_begin          | Initial sequence number of  |
    |         |             |                    | PDB sequence segment.       |
    +---------+-------------+--------------------+-----------------------------+
    | 19      | AChar       | ins_begin          | Initial insertion code of   |
    |         |             |                    | PDB sequence segment.       |
    +---------+-------------+--------------------+-----------------------------+
    | 21-24   | Integer     | seq_end            | Ending sequence number of   |
    |         |             |                    | PDB sequence segment.       |
    +---------+-------------+--------------------+-----------------------------+
    | 25      | AChar       | ins_end            | Ending insertion code of    |
    |         |             |                    | PDB sequence segment.       |
    +---------+-------------+--------------------+-----------------------------+
    | 27-32   | LString     | database           | Sequence database name.     |
    +---------+-------------+--------------------+-----------------------------+
    | 34-41   | LString     | database_accession | Sequence database accession |
    |         |             |                    | code.                       |
    +---------+-------------+--------------------+-----------------------------+
    | 43-54   | LString     | database_id_code   | Sequence database id        |
    |         |             |                    | code.                       |
    +---------+-------------+--------------------+-----------------------------+
    | 56-60   | Integer     | database_seq_begin | Initial sequence number of  |
    |         |             |                    | database seqment.           |
    +---------+-------------+--------------------+-----------------------------+
    | 61      | AChar       | database_ins_begin | Insertion code of initial   |
    |         |             |                    | residue segment, if PDB is  |
    |         |             |                    | reference.                  |
    +---------+-------------+--------------------+-----------------------------+
    | 63-67   | Integer     | database_seq_end   | Ending sequence number of   |
    |         |             |                    | segment.                    |
    +---------+-------------+--------------------+-----------------------------+
    | 68      | AChar       | database_ins_end   | Insertion code of the       |
    |         |             |                    | the segment end, if PDB is  |
    |         |             |                    | reference.                  |
    +---------+-------------+--------------------+-----------------------------+
    """

    def __init__(self):
        super().__init__()
        self.id_code = ""
        self.chain_id = ""
        self.seq_begin = ""
        self.ins_begin = ""
        self.seq_end = ""
        self.ins_end = ""
        self.database = ""
        self.database_accession = ""
        self.database_id_code = ""
        self.database_seq_begin = ""
        self.database_ins_begin = ""
        self.database_seq_end = ""
        self.database_ins_end = ""

    @staticmethod
    def parse_cif(container) -> list:
        """Parse CIF container for information about this record.

        :param :class:`pdbx.containers.DataContainer` container:  container to
            parse
        :returns:  list of :class:`DatabaseReference` or similar objects
        """
        db_refs = []
        struct_ref_seq_df = cif_df(container.get_object("struct_ref_seq"))
        struct_ref_df = cif_df(container.get_object("struct_ref"))
        df = struct_ref_df.merge(
            struct_ref_seq_df,
            how="left",
            left_on="id",
            right_on="align_id",
            suffixes=("", "_seq"),
        )
        for _, row in df.iterrows():
            row = row.dropna()
            id_code = row["pdbx_PDB_id_code"]
            chain_id = row["pdbx_strand_id"]
            seq_begin = row["seq_align_beg"]
            try:
                ins_begin = row["pdbx_seq_align_beg_ins_code"]
            except KeyError:
                ins_begin = ""
            seq_end = row["seq_align_end"]
            try:
                ins_end = row["pdbx_seq_align_end_ins_code"]
            except KeyError:
                ins_end = ""
            database = row["db_name"]
            database_accession = row["pdbx_db_accession"]
            database_id_code = row["db_code"]
            database_seq_begin = row["db_align_beg"]
            database_seq_end = row["db_align_end"]
            if len(database_accession) < 13:
                ref = DatabaseReference()
                ref.id_code = id_code
                ref.chain_id = chain_id
                ref.seq_begin = seq_begin
                ref.ins_begin = ins_begin
                ref.seq_end = seq_end
                ref.ins_end = ins_end
                ref.database = database
                ref.database_accession = database_accession
                ref.database_id_code = database_id_code
                ref.database_seq_begin = database_seq_begin
                ref.database_seq_end = database_seq_end
                db_refs.append(ref)
            else:
                ref1 = DatabaseReference1()
                ref1.id_code = id_code
                ref1.chain_id = chain_id
                ref1.seq_begin = seq_begin
                ref1.ins_begin = ins_begin
                ref1.seq_end = seq_end
                ref1.ins_end = ins_end
                ref1.database = database
                ref1.database_id_code = database_id_code
                db_refs.append(ref1)
                ref2 = DatabaseReference2()
                ref2.id_code = id_code
                ref2.chain_id = chain_id
                ref2.database_accession = database_accession
                ref2.database_seq_begin = database_seq_begin
                ref2.database_seq_end = database_seq_end
                db_refs.append(ref2)
        return db_refs

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
        str = f"DBREF  {self.id_code:4} {self.chain_id:1} {self.seq_begin:>4}"
        str += f"{self.ins_begin:1} {self.seq_end:>4}{self.ins_end:1} "
        str += f"{self.database:6} {self.database_accession:8} "
        str += f"{self.database_id_code:12} {self.database_seq_begin:>5}"
        str += f"{self.database_ins_begin:1} {self.database_seq_end:>5}"
        str += f"{self.database_ins_end:1}"
        return str


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
    | 8-11    | IDcode      | id_code     | ID code of this entry.            |
    +---------+-------------+-------------+-----------------------------------+
    | 13      | Character   | chain_id    | Chain identifier.                 |
    +---------+-------------+-------------+-----------------------------------+
    | 15-18   | Integer     | seq_begin   | Initial sequence number of the    |
    |         |             |             | PDB sequence segment, right       |
    |         |             |             | justified.                        |
    +---------+-------------+-------------+-----------------------------------+
    | 19      | AChar       | ins_begin   | Initial insertion code of the PDB |
    |         |             |             | sequence segment.                 |
    +---------+-------------+-------------+-----------------------------------+
    | 21-24   | Integer     | seq_end     | Ending sequence number of the PDB |
    |         |             |             | sequence segment, right           |
    |         |             |             | justified.                        |
    +---------+-------------+-------------+-----------------------------------+
    | 25      | AChar       | ins_end     | Ending insertion code of the PDB  |
    |         |             |             | sequence segment.                 |
    +---------+-------------+-------------+-----------------------------------+
    | 27-32   | LString     | database    | Sequence database name.           |
    +---------+-------------+-------------+-----------------------------------+
    | 48-67   | LString     | db_id_code  | Sequence database identification  |
    |         |             |             | code, left justified.             |
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
        self.database_id_code = None

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
        self.database_id_code = line[47:67].strip()

    def __str__(self):
        return (
            f"DBREF1 {self.id_code:4} {self.chain_id:1} {self.seq_begin:4}"
            f"{self.ins_begin:1} {self.seq_end:4}{self.ins_end:1}"
            f" {self.database:6}               {self.database_id_code:20}"
        )


class DatabaseReference2(BaseRecord):
    """Provides cross-reference links between PDB sequences (what appears in
    SEQRES record) and a corresponding database sequence.

    This updated two-line format is used when the accession code or sequence
    numbering does not fit the space allotted in the standard DBREF format.
    This includes some GenBank sequence numbering (greater than 5 characters)
    and UNIMES accession numbers (greater than 12 characters).

    +---------+-------------+--------------+----------------------------------+
    | COLUMNS | DATA TYPE   | FIELD        | DEFINITION                       |
    +=========+=============+==============+==================================+
    | 1-6     | Record name | "DBREF2"     |                                  |
    +---------+-------------+--------------+----------------------------------+
    | 8-11    | IDcode      | id_code      | ID code of this entry.           |
    +---------+-------------+--------------+----------------------------------+
    | 13      | Character   | chain_id     | Chain identifier.                |
    +---------+-------------+--------------+----------------------------------+
    | 19-40   | LString     | db_accession | Sequence database accession code |
    |         |             |              | left justified.                  |
    +---------+-------------+--------------+----------------------------------+
    | 46-55   | Integer     | seq_begin    | Initial sequence number of the   |
    |         |             |              | Database segment, right          |
    |         |             |              | justified.                       |
    +---------+-------------+--------------+----------------------------------+
    | 58-67   | Integer     | seq_end      | Ending sequence number of the    |
    |         |             |              | Database segment, right          |
    |         |             |              | justified.                       |
    +---------+-------------+--------------+----------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.id_code = None
        self.chain_id = None
        self.database_accession = None
        self.database_seq_begin = None
        self.database_seq_end = None

    def parse_line(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_line(line)
        self.id_code = line[7:11].strip()
        self.chain_id = line[12].strip()
        self.database_accession = line[18:40].strip()
        self.database_seq_begin = int(line[45:55])
        self.seq_end = int(line[57:67])

    def __str__(self):
        return (
            f"DBREF2 {self.id_code:4} {self.chain_id:1}"
            f"     {self.database_accession:22}     {self.database_seq_begin:10}"
            f"  {self.seq_end:10}"
        )


class ModifiedResidue(BaseRecord):
    """MODRES field

    The MODRES record provides descriptions of modifications (e.g.,
    chemical or post-translational) to protein and nucleic acid residues.
    Included are a mapping between residue names given in a PDB entry and
    standard residues.

    +---------+--------------+--------------+---------------------------------+
    | COLUMNS | DATA TYPE    | FIELD        | DEFINITION                      |
    +=========+==============+==============+=================================+
    | 1-6     | Record name  | "MODRES"     |                                 |
    +---------+--------------+--------------+---------------------------------+
    | 8-11    | IDcode       | id_code      | ID code of this entry.          |
    +---------+--------------+--------------+---------------------------------+
    | 13-15   | Residue name | res_name     | Residue name used in this entry |
    +---------+--------------+--------------+---------------------------------+
    | 17      | Character    | chain_id     | Chain identifier.               |
    +---------+--------------+--------------+---------------------------------+
    | 19-22   | Integer      | seq_num      | Sequence number.                |
    +---------+--------------+--------------+---------------------------------+
    | 23      | AChar        | ins_code     | Insertion code.                 |
    +---------+--------------+--------------+---------------------------------+
    | 25-27   | Residue name | standard_res | Standard residue name.          |
    +---------+--------------+--------------+---------------------------------+
    | 30-70   | String       | comment      | Description of the residue      |
    |         |              |              | modification.                   |
    +---------+--------------+--------------+---------------------------------+
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

    @staticmethod
    def parse_cif(container) -> list:
        """Parse CIF container for information about this record.

        :param :class:`pdbx.containers.DataContainer` container:  container to
            parse
        :returns:  list of :class:`SequenceDifferences` objects
        """
        mod_res = []
        df = cif_df(container.get_object("pdbx_struct_mod_residue"))
        if len(df) > 0:
            raise NotImplementedError("MODRES parsing not implemented.")
        return mod_res

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
    | 8-11    | IDcode       | id_code     | ID code of this entry.           |
    +---------+--------------+-------------+----------------------------------+
    | 13-15   | Residue name | res_name    | Name of the PDB residue in       |
    |         |              |             | conflict.                        |
    +---------+--------------+-------------+----------------------------------+
    | 17      | Character    | chain_id    | PDB chain identifier.            |
    +---------+--------------+-------------+----------------------------------+
    | 19-22   | Integer      | seq_num     | PDB sequence number.             |
    +---------+--------------+-------------+----------------------------------+
    | 23      | AChar        | ins_code    | PDB insertion code.              |
    +---------+--------------+-------------+----------------------------------+
    | 25-28   | LString      | database    |                                  |
    +---------+--------------+-------------+----------------------------------+
    | 30-38   | LString      | db_id_code  | Sequence database accession      |
    |         |              |             | number.                          |
    +---------+--------------+-------------+----------------------------------+
    | 40-42   | Residue name | db_res      | Sequence database residue name.  |
    +---------+--------------+-------------+----------------------------------+
    | 44-48   | Integer      | db_seq      | Sequence database sequence       |
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

    @staticmethod
    def parse_cif(container) -> list:
        """Parse CIF container for information about this record.

        :param :class:`pdbx.containers.DataContainer` container:  container to
            parse
        :returns:  list of :class:`SequenceDifferences` objects
        """
        diffs = []
        cif_obj = container.get_object("struct_ref_seq_dif")
        if cif_obj is None:
            return diffs
        df = cif_df(cif_obj)
        print(df)
        raise NotImplementedError()

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
        self._residues = OrderedDict()
        self.num_residues = {}

    @property
    def residues(self) -> dict:
        """Dictionary of residues indexed by chain id.

        :returns: dictionary with chain IDs as keys and lists of residue names
            as values.
        """
        return self._residues

    def parse_cif(self, container) -> bool:
        """Parse CIF container for information about this record.

        :param :class:`pdbx.containers.DataContainer` container:  container to
            parse
        :returns:  True if useful information was extracted from container
        """
        value_added = False
        pdbx_poly_seq_df = cif_df(
            container.get_object("pdbx_poly_seq_scheme")
        )
        groups = pdbx_poly_seq_df.groupby("pdb_strand_id")
        for strand_id, df in groups:
            if len(df["mon_id"].values) > 0:
                value_added = True
                self._residues[strand_id] = df["mon_id"].values
        entity_poly_seq_df = cif_df(
            container.get_object("entity_poly_seq_df")
        )
        if len(entity_poly_seq_df) > 0:
            raise NotImplementedError(
                "Parsing entity_poly_seq not implemented."
            )
        return value_added

    @residues.setter
    def residues(self, value):
        self._residues = value

    def parse_line(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_line(line)
        chain_id = line[11].strip()
        if chain_id not in self._residues:
            self._residues[chain_id] = []
        self.num_residues[chain_id] = int(line[13:17])
        start = 19
        end = 22
        while True:
            self._residues[chain_id].append(line[start:end].strip())
            start += 4
            end += 4
            if start > 68:
                break

    def num_chains(self) -> int:
        """Number of chains in sequence."""
        return len(self._residues)

    def __str__(self):
        strings = []
        for chain_id, residues in self._residues.items():
            for ichunk, chunk in enumerate(grouper(residues, 13)):
                serial_num = ichunk + 1
                string = f"SEQRES {serial_num:>3} {chain_id:1} "
                string += f"{self.num_residues[chain_id]:>4} "
                for residue in chunk:
                    string += f" {residue:>3}"
                strings.append(string.strip())
        return "\n".join(strings)
