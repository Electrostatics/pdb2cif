"""Classes for records with secondary structure and connectivity information.

.. codeauthor::  Todd Dolinsky
.. codeauthor::  Yong Huang
.. codeauthor::  Nathan Baker
"""
import logging
from .general import BaseRecord, atom_format, cif_df


_LOGGER = logging.getLogger(__name__)


class CisPeptide(BaseRecord):
    """CISPEP field

    CISPEP records specify the prolines and other peptides found to be in
    the cis conformation. This record replaces the use of footnote records
    to list cis peptides.

    +---------+-------------+-----------+-------------------------------------+
    | COLUMNS | DATA TYPE   | FIELD     | DEFINITION                          |
    +=========+=============+===========+=====================================+
    | 1-6     | Record name | "CISPEP"  |                                     |
    +---------+-------------+-----------+-------------------------------------+
    | 8-10    | Integer     | ser_num   | Record serial number.               |
    +---------+-------------+-----------+-------------------------------------+
    | 12-14   | LString(3)  | pep1      | Residue name.                       |
    +---------+-------------+-----------+-------------------------------------+
    | 16      | Character   | chain_id1 | Chain identifier.                   |
    +---------+-------------+-----------+-------------------------------------+
    | 18-21   | Integer     | seq_num1  | Residue sequence number.            |
    +---------+-------------+-----------+-------------------------------------+
    | 22      | AChar       | icode1    | Insertion code.                     |
    +---------+-------------+-----------+-------------------------------------+
    | 26-28   | LString(3)  | pep2      | Residue name.                       |
    +---------+-------------+-----------+-------------------------------------+
    | 30      | Character   | chain_id2 | Chain identifier.                   |
    +---------+-------------+-----------+-------------------------------------+
    | 32-35   | Integer     | seq_num2  | Residue sequence number.            |
    +---------+-------------+-----------+-------------------------------------+
    | 36      | AChar       | icode2    | Insertion code.                     |
    +---------+-------------+-----------+-------------------------------------+
    | 44-46   | Integer     | mod_num   | Identifies the specific model.      |
    +---------+-------------+-----------+-------------------------------------+
    | 54-59   | Real(6.2)   | measure   | Angle measurement in degrees.       |
    +---------+-------------+-----------+-------------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.ser_num = None
        self.pep1 = None
        self.chain_id1 = None
        self.seq_num1 = None
        self.icode1 = None
        self.pep2 = None
        self.chain_id2 = None
        self.seq_num2 = None
        self.icode2 = None
        self.mod_num = None
        self.measure = None

    @staticmethod
    def parse_cif(container) -> list:
        """Parse CIF container for information about this record.

        :param :class:`pdbx.containers.DataContainer` container:  container to
            parse
        :returns:  list of objects of this class
        """
        cis_peps = []
        df = cif_df(container.get_object("struct_mon_prot_cis")).fillna("")
        for _, row in df.iterrows():
            pep = CisPeptide()
            pep.ser_num = int(row["pdbx_id"])
            pep.pep1 = row["auth_comp_id"]
            pep.chain_id1 = row["auth_asym_id"]
            pep.seq_num1 = int(row["auth_seq_id"])
            pep.icode1 = row["pdbx_PDB_ins_code"]
            pep.pep2 = row["pdbx_auth_comp_id_2"]
            pep.chain_id2 = row["pdbx_auth_asym_id_2"]
            pep.seq_num2 = int(row["pdbx_auth_seq_id_2"])
            pep.icode2 = row["pdbx_PDB_ins_code_2"]
            pep.mod_num = int(row["pdbx_PDB_model_num"])
            pep.measure = float(row["pdbx_omega_angle"])
            cis_peps.append(pep)
        return cis_peps

    def parse_pdb(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_pdb(line)
        self.ser_num = int(line[7:10].strip())
        self.pep1 = line[11:14].strip()
        self.chain_id1 = line[15].strip()
        self.seq_num1 = int(line[17:21].strip())
        self.icode1 = line[21].strip()
        self.pep2 = line[25:28].strip()
        self.chain_id2 = line[29].strip()
        self.seq_num2 = int(line[31:35].strip())
        self.icode2 = line[35].strip()
        self.mod_num = int(line[43:46].strip())
        self.measure = float(line[53:59].strip())

    def __str__(self):
        str_ = f"CISPEP {self.ser_num:3} {self.pep1:3} {self.chain_id1:1}"
        str_ += f" {self.seq_num1:4}{self.icode1:1}   {self.pep2:3}"
        str_ += f" {self.chain_id2:1} {self.seq_num2:4}{self.icode2:1}"
        str_ += f"       {self.mod_num:3}       {self.measure:6.2f}"
        return str_


class DisulfideBond(BaseRecord):
    """SSBOND field

    The SSBOND record identifies each disulfide bond in protein and
    polypeptide structures by identifying the two residues involved in the
    bond.

    +---------+-------------+-----------+-------------------------------------+
    | COLUMNS | DATA TYPE   | FIELD     | DEFINITION                          |
    +=========+=============+===========+=====================================+
    | 1-6     | Record name | "SSBOND"  |                                     |
    +---------+-------------+-----------+-------------------------------------+
    | 8-10    | Integer     | ser_num   | Serial number.                      |
    +---------+-------------+-----------+-------------------------------------+
    | 12-14   | LString(3)  | "CYS"     | Residue name.                       |
    +---------+-------------+-----------+-------------------------------------+
    | 16      | Character   | chain_id1 | Chain identifier.                   |
    +---------+-------------+-----------+-------------------------------------+
    | 18-21   | Integer     | seq_num1  | Residue sequence number.            |
    +---------+-------------+-----------+-------------------------------------+
    | 22      | AChar       | icode1    | Insertion code.                     |
    +---------+-------------+-----------+-------------------------------------+
    | 26-28   | LString(3)  | "CYS"     | Residue name.                       |
    +---------+-------------+-----------+-------------------------------------+
    | 30      | Character   | chain_id2 | Chain identifier.                   |
    +---------+-------------+-----------+-------------------------------------+
    | 32-35   | Integer     | seq_num2  | Residue sequence number.            |
    +---------+-------------+-----------+-------------------------------------+
    | 36      | AChar       | icode2    | Insertion code.                     |
    +---------+-------------+-----------+-------------------------------------+
    | 60-65   | SymOP       | sym1      | Symmetry operator for residue 1.    |
    +---------+-------------+-----------+-------------------------------------+
    | 67-72   | SymOP       | sym2      | Symmetry operator for residue 2.    |
    +---------+-------------+-----------+-------------------------------------+
    | 74-78   | Real(5.2)   | length    | Disulfide bond distance             |
    +---------+-------------+-----------+-------------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.ser_num = None
        self.chain_id1 = None
        self.seq_num1 = None
        self.icode1 = None
        self.chain_id2 = None
        self.seq_num2 = None
        self.icode2 = None
        self.sym1 = None
        self.sym2 = None
        self.length = None

    @staticmethod
    def parse_cif(container) -> list:
        """Parse CIF container for information about this record.

        :param :class:`pdbx.containers.DataContainer` container:  container to
            parse
        :returns:  list of objects of this class
        """
        bonds = []
        df = cif_df(container.get_object("struct_conn")).fillna("")
        for id_, row in df.iterrows():
            bond = DisulfideBond()
            bond.ser_num = id_ + 1
            bond.chain_id1 = row["ptnr1_auth_asym_id"]
            bond.seq_num1 = row["ptnr1_auth_seq_id"]
            bond.icode1 = row["pdbx_ptnr1_PDB_ins_code"]
            bond.chain_id2 = row["ptnr2_auth_asym_id"]
            bond.seq_num2 = row["ptnr2_auth_seq_id"]
            bond.icode2 = row["pdbx_ptnr2_PDB_ins_code"]
            bond.sym1 = row["ptnr1_symmetry"]
            bond.sym2 = row["ptnr2_symmetry"]
            distance = row["pdbx_dist_value"]
            if distance:
                bond.length = float(distance)
            bonds.append(bond)
        return bonds

    def parse_pdb(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_pdb(line)
        self.ser_num = int(line[7:10].strip())
        self.chain_id1 = line[15].strip()
        self.seq_num1 = int(line[17:21].strip())
        self.icode1 = line[21].strip()
        self.chain_id2 = line[29].strip()
        self.seq_num2 = int(line[31:35].strip())
        self.icode2 = line[35].strip()
        self.sym1 = line[59:65].strip()
        self.sym2 = line[66:72].strip()
        self.length = float(line[73:78])

    def __str__(self):
        str_ = (
            f"SSBOND {self.ser_num:3} CYS {self.chain_id1:1} {self.seq_num1:4}"
            f"{self.icode1:1}   CYS {self.chain_id2:1} {self.seq_num2:4}"
            f"{self.icode2:1}                         {self.sym1:6}"
            f" {self.sym2:6}"
        )
        if self.length:
            str_ += f"{self.length:4.2f}"
        return str_


class Helix(BaseRecord):
    """HELIX field

    HELIX records are used to identify the position of helices in the
    molecule. Helices are both named and numbered. The residues where the
    helix begins and ends are noted, as well as the total length.

    +---------+--------------+---------------+--------------------------------+
    | COLUMNS | DATA TYPE    | FIELD         | DEFINITION                     |
    +=========+==============+===============+================================+
    | 1-6     | Record name  | "HELIX "      |                                |
    +---------+--------------+---------------+--------------------------------+
    | 8-10    | Integer      | serNum        | Serial number of the helix.    |
    |         |              |               | starts at 1 and increases      |
    |         |              |               | incrementally.                 |
    +---------+--------------+---------------+--------------------------------+
    | 12-14   | LString(3)   | helix_id      | Helix identifier. In addition  |
    |         |              |               | to a serial number, each helix |
    |         |              |               | is given an alphanumeric       |
    |         |              |               | helix identifier.              | 
    +---------+--------------+---------------+--------------------------------+
    | 16-18   | Residue name | init_res_name | Name of the initial residue.   |
    +---------+--------------+---------------+--------------------------------+
    | 20      | Character    | init_chain_id | Chain identifier for the chain |
    |         |              |               | containing this helix.         |
    +---------+--------------+---------------+--------------------------------+
    | 22-25   | Integer      | init_seq_num  | Sequence number of the initial |
    |         |              |               | residue.                       |
    +---------+--------------+---------------+--------------------------------+
    | 26      | AChar        | init_i_code   | Insertion code of the initial  |
    |         |              |               | residue.                       |
    +---------+--------------+---------------+--------------------------------+
    | 28-30   | Residue name | end_res_name  | Name of the terminal residue   |
    |         |              |               | of the helix.                  |
    +---------+--------------+---------------+--------------------------------+
    | 32      | Character    | end_chain_id  | Chain identifier for the chain |
    |         |              |               | containing this helix.         |
    +---------+--------------+---------------+--------------------------------+
    | 34-37   | Integer      | end_seq_num   | Sequence number of the         |
    |         |              |               | terminal residue.              |
    +---------+--------------+---------------+--------------------------------+
    | 38      | AChar        | end_i_code    | Insertion code of the terminal |
    |         |              |               | residue.                       |
    +---------+--------------+---------------+--------------------------------+
    | 39-40   | Integer      | helix_class   | Helix class (see below).       |
    +---------+--------------+---------------+--------------------------------+
    | 41-70   | String       | comment       | Comment about this helix.      |
    +---------+--------------+---------------+--------------------------------+
    | 72-76   | Integer      | length        | Length of this helix.          |
    +---------+--------------+---------------+--------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.ser_num = None
        self.helix_id = None
        self.init_res_name = None
        self.init_chain_id = None
        self.init_seq_num = None
        self.init_i_code = None
        self.end_res_name = None
        self.end_chain_id = None
        self.end_seq_num = None
        self.end_i_code = None
        self.helix_class = None
        self.comment = None
        self.length = None

    @staticmethod
    def parse_cif(container) -> list:
        """Parse CIF container for information about this record.

        :param :class:`pdbx.containers.DataContainer` container:  container to
            parse
        :returns:  list of objects of this class
        """
        helices = []
        df = cif_df(container.get_object("struct_conf"))
        df = df.fillna("")
        for _, row in df.iterrows():
            helix = Helix()
            helix.ser_num = row["pdbx_PDB_helix_id"]
            helix.helix_id = row["pdbx_PDB_helix_id"]
            helix.init_res_name = row["beg_auth_comp_id"]
            helix.init_chain_id = row["beg_auth_asym_id"]
            helix.init_seq_num = row["beg_auth_seq_id"]
            helix.init_i_code = row["pdbx_beg_PDB_ins_code"]
            helix.end_res_name = row["end_auth_comp_id"]
            helix.end_chain_id = row["end_auth_asym_id"]
            helix.end_seq_num = row["end_auth_seq_id"]
            helix.end_i_code = row["pdbx_end_PDB_ins_code"]
            helix.helix_class = row["pdbx_PDB_helix_class"]
            helix.comment = row["details"]
            helix.length = row["pdbx_PDB_helix_length"]
            helices.append(helix)
        return helices

    def parse_pdb(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_pdb(line)
        self.ser_num = int(line[7:10].strip())
        self.helix_id = line[11:14].strip()
        self.init_res_name = line[15:18].strip()
        self.init_chain_id = line[19].strip()
        self.init_seq_num = int(line[21:25].strip())
        self.init_i_code = line[25].strip()
        self.end_res_name = line[27:30].strip()
        self.end_chain_id = line[31].strip()
        self.end_seq_num = int(line[33:37].strip())
        self.end_i_code = line[37].strip()
        try:
            self.helix_class = int(line[38:40].strip())
        except ValueError:
            pass
        self.comment = line[40:70].strip()
        try:
            self.length = int(line[71:76].strip())
        except ValueError:
            pass

    def __str__(self):
        return (
            f"HELIX  {self.ser_num:3} {self.helix_id:>3}"
            f" {self.init_res_name:3} {self.init_chain_id:1}"
            f" {self.init_seq_num:4}{self.init_i_code:1}"
            f" {self.end_res_name:3} {self.end_chain_id:1}"
            f" {self.end_seq_num:4}{self.end_i_code:1}{self.helix_class:2}"
            f"{self.comment:30} {self.length:5}"
        )


class Link(BaseRecord):
    """LINK field

    The LINK records specify connectivity between residues that is not
    implied by the primary structure. Connectivity is expressed in terms of
    the atom names. This record supplements information given in CONECT
    records and is provided here for convenience in searching.

    .. todo::  Clean up output using element/atom information and formatting

    +---------+--------------+-----------+------------------------------------+
    | COLUMNS | DATA TYPE    | FIELD     | DEFINITION                         |
    +=========+==============+===========+====================================+
    | 1-6     | Record name  | "LINK  "  |                                    |
    +---------+--------------+-----------+------------------------------------+
    | 13-16   | Atom         | name1     | Atom name.                         |
    +---------+--------------+-----------+------------------------------------+
    | 17      | Character    | alt_loc1  | Alternate location indicator.      |
    +---------+--------------+-----------+------------------------------------+
    | 18-20   | Residue name | res_name1 | Residue  name.                     |
    +---------+--------------+-----------+------------------------------------+
    | 22      | Character    | chain_id  | Chain identifier.                  |
    +---------+--------------+-----------+------------------------------------+
    | 23-26   | Integer      | res_seq1  | Residue sequence number.           |
    +---------+--------------+-----------+------------------------------------+
    | 27      | AChar        | ins_code1 | Insertion code.                    |
    +---------+--------------+-----------+------------------------------------+
    | 43-46   | Atom         | name2     | Atom name.                         |
    +---------+--------------+-----------+------------------------------------+
    | 47      | Character    | alt_loc2  | Alternate location indicator.      |
    +---------+--------------+-----------+------------------------------------+
    | 48-50   | Residue name | res_name2 | Residue name.                      |
    +---------+--------------+-----------+------------------------------------+
    | 52      | Character    | chain_id  | Chain identifier.                  |
    +---------+--------------+-----------+------------------------------------+
    | 53-56   | Integer      | res_seq2  | Residue sequence number.           |
    +---------+--------------+-----------+------------------------------------+
    | 57      | AChar        | ins_code2 | Insertion code.                    |
    +---------+--------------+-----------+------------------------------------+
    | 60-65   | SymOP        | sym1      | Symmetry operator atom 1.          |
    +---------+--------------+-----------+------------------------------------+
    | 67-72   | SymOP        | sym2      | Symmetry operator atom 2.          |
    +---------+--------------+-----------+------------------------------------+
    | 74-78   | Real(5.2)    | Length    | Link distance                      |
    +---------+--------------+-----------+------------------------------------+
    """

    def __init__(self):
        super().__init__()
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
        self.length = None

    @staticmethod
    def parse_cif(container) -> list:
        """Parse CIF container for information about this record.

        :param :class:`pdbx.containers.DataContainer` container:  container to
            parse
        :returns:  list of objects of this class
        """
        links = []
        df = cif_df(container.get_object("struct_conn")).fillna("")
        for _, row in df.iterrows():
            link = Link()
            if row["conn_type_id"] not in ["disulf"]:
                link.name1 = row["ptnr1_label_atom_id"]
                link.alt_loc1 = row["pdbx_ptnr1_label_alt_id"]
                link.res_name1 = row["ptnr1_auth_comp_id"]
                link.chain_id1 = row["ptnr1_auth_asym_id"]
                link.res_seq1 = row["ptnr1_auth_seq_id"]
                link.ins_code1 = row["pdbx_ptnr1_PDB_ins_code"]
                link.name2 = row["ptnr2_label_atom_id"]
                link.alt_loc2 = row["pdbx_ptnr2_label_alt_id"]
                link.res_name2 = row["ptnr2_auth_comp_id"]
                link.chain_id2 = row["ptnr2_auth_asym_id"]
                link.res_seq2 = row["ptnr2_auth_seq_id"]
                link.ins_code2 = row["pdbx_ptnr2_PDB_ins_code"]
                link.sym1 = row["ptnr1_symmetry"]
                link.sym2 = row["ptnr2_symmetry"]
                distance = row["pdbx_dist_value"]
                if distance:
                    link.length = float(distance)
                links.append(link)
        return links

    def parse_pdb(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_pdb(line)
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
        self.length = line[73:78]
        self.is_element1 = None
        self.is_element2 = None

    def __str__(self):
        # See atom-name formatting rules at
        # https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
        if None in [self.is_element1, self.is_element2]:
            err = (
                "Must run annotate_link first before correctly formatted "
                "strings can be produced."
            )
            raise ValueError(err)
        if self.is_element1:
            name1 = f"{self.name1:>2}  "[:4]
        elif len(self.name1) == 2:
            name1 = f" {self.name1} "
        else:
            name1 = f"{self.name1:>4}"
        if self.is_element2:
            name2 = f"{self.name2:>2}  "[:4]
        elif len(self.name2) == 2:
            name2 = f" {self.name2} "
        else:
            name2 = f"{self.name2:>4}"
        string = f"LINK        "
        string += (
            f"{name1}{self.alt_loc1:1}{self.res_name1:>3} {self.chain_id1:1}"
            f"{self.res_seq1:4}{self.ins_code1:1}               "
        )
        string += (
            f"{name2}{self.alt_loc2:1}{self.res_name2:>3} {self.chain_id2}"
            f"{self.res_seq2:4}{self.ins_code2:1}  {self.sym1:>6} "
            f"{self.sym2:>6}"
        )
        if self.length:
            string += f" {self.length:5}"
        return string


class Sheet(BaseRecord):
    """SHEET field

    SHEET records are used to identify the position of sheets in the
    molecule. Sheets are both named and numbered. The residues where the
    sheet begins and ends are noted.

    +---------+--------------+---------------+--------------------------------+
    | COLUMNS | DATA TYPE    | FIELD         | DEFINITION                     |
    +=========+==============+===============+================================+
    | 1-6     | Record name  | "SHEET "      |                                |
    +---------+--------------+---------------+--------------------------------+
    | 8-10    | Integer      | strand        | Strand number which starts at  |
    |         |              |               | 1 for each strand within a     |
    |         |              |               | sheet and increases by one.    |
    +---------+--------------+---------------+--------------------------------+
    | 12-14   | LString(3)   | sheet_id      | Sheet identifier.              |
    +---------+--------------+---------------+--------------------------------+
    | 15-16   | Integer      | num_strands   | Number of strands in sheet.    |
    +---------+--------------+---------------+--------------------------------+
    | 18-20   | Residue name | init_res_name | Name of initial residue.       |
    +---------+--------------+---------------+--------------------------------+
    | 22      | Character    | init_chain_id | Chain identifier of initial    |
    |         |              |               | residue in strand.             |
    +---------+--------------+---------------+--------------------------------+
    | 23-26   | Integer      | init_seq_num  | Sequence number of initial     |
    |         |              |               | residue in strand.             |
    +---------+--------------+---------------+--------------------------------+
    | 27      | AChar        | init_ins_code | Insertion code of initial      |
    |         |              |               | residue in strand.             |
    +---------+--------------+---------------+--------------------------------+
    | 29-31   | Residue name | end_res_name  | Name of terminal residue       |
    +---------+--------------+---------------+--------------------------------+
    | 33      | Character    | end_chain_id  | Chain identifier of terminal   |
    |         |              |               | residue                        |
    +---------+--------------+---------------+--------------------------------+
    | 34-37   | Integer      | end_seq_num   | Sequence number of terminal    |
    |         |              |               | residue.                       |
    +---------+--------------+---------------+--------------------------------+
    | 38      | AChar        | end_ins_code  | Insertion code of terminal     |
    |         |              |               | residue.                       |
    +---------+--------------+---------------+--------------------------------+
    | 39-40   | Integer      | sense         | Sense of strand with respect   |
    |         |              |               | to previous strand in the      |
    |         |              |               | sheet. 0 if first strand, 1 if |
    |         |              |               | parallel, and -1 if            |
    |         |              |               | anti-parallel.                 |
    +---------+--------------+---------------+--------------------------------+
    | 42-45   | Atom         | cur_atom      | Registration. Atom name in     |
    |         |              |               | current strand.                |
    +---------+--------------+---------------+--------------------------------+
    | 46-48   | Residue name | cur_res_name  | Registration. Residue name in  |
    |         |              |               | current strand                 |
    +---------+--------------+---------------+--------------------------------+
    | 50      | Character    | cur_chain_id  | Registration. Chain identifier |
    |         |              |               | in current strand.             |
    +---------+--------------+---------------+--------------------------------+
    | 51-54   | Integer      | cur_res_seq   | Registration. Residue sequence |
    |         |              |               | number in current strand.      |
    +---------+--------------+---------------+--------------------------------+
    | 55      | AChar        | cur_ins_code  | Registration. Insertion code   |
    |         |              |               | in current strand.             |
    +---------+--------------+---------------+--------------------------------+
    | 57-60   | Atom         | prev_atom     | Registration. Atom name in     |
    |         |              |               | previous strand.               |
    +---------+--------------+---------------+--------------------------------+
    | 61-63   | Residue name | prev_res_name | Registration. Residue name in  |
    |         |              |               | previous strand.               |
    +---------+--------------+---------------+--------------------------------+
    | 65      | Character    | prev_chain_id | Registration. Chain identifier |
    |         |              |               | in previous strand.            |
    +---------+--------------+---------------+--------------------------------+
    | 66-69   | Integer      | prev_res_seq  | Registration. Residue sequence |
    |         |              |               | number in previous strand.     |
    +---------+--------------+---------------+--------------------------------+
    | 70      | AChar        | prev_ins_code | Registration. Insertion code   |
    |         |              |               | in previous strand.            |
    +---------+--------------+---------------+--------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.range_id = None
        self.sheet_id = None
        self.num_strands = None
        self.init_res_name = None
        self.init_chain_id = None
        self.init_seq_num = None
        self.init_ins_code = None
        self.end_res_name = None
        self.end_chain_id = None
        self.end_seq_num = None
        self.end_ins_code = None
        self.sense = None
        self.curr_atom = ""
        self.curr_res_name = ""
        self.curr_chain_id = ""
        self.curr_res_seq = ""
        self.curr_ins_code = ""
        self.prev_atom = ""
        self.prev_res_name = ""
        self.prev_chain_id = ""
        self.prev_res_seq = ""
        self.prev_ins_code = ""

    # @staticmethod
    # def parse_cif(container) -> list:
    #     """Parse CIF container for information about this record.

    #     :param :class:`pdbx.containers.DataContainer` container:  container to
    #         parse
    #     :returns:  list of objects of this class
    #     """
    #     sheets = []
    #     hbond_df = cif_df(
    #         container.get_object("pdbx_struct_sheet_hbond")
    #     ).fillna("")
    #     range_df = cif_df(container.get_object("struct_sheet_range")).fillna(
    #         ""
    #     )
    #     order_df = cif_df(container.get_object("struct_sheet_order")).fillna(
    #         ""
    #     )
    #     sheet_df = cif_df(container.get_object("struct_sheet")).fillna("")
    #     sheet_ids = set()
    #     range_ids = set()
    #     if len(hbond_df) > 0:
    #         sheet_ids |= set(hbond_df["sheet_id"])
    #         range_ids |= set(hbond_df["range_id_2"])
    #     if len(sheet_df) > 0:
    #         sheet_ids |= set(sheet_df["id"])
    #     if len(order_df) > 0:
    #         sheet_ids |= set(order_df["sheet_id"])
    #         range_ids |= set(order_df["range_id_1"])
    #         range_ids |= set(order_df["range_id_2"])
    #     if len(range_df) > 0:
    #         sheet_ids |= set(range_df["sheet_id"])
    #         range_ids |= set(range_df["id"])
    #     sheet_ids = sorted(list(sheet_ids))
    #     range_ids = sorted(list(range_ids))
    #     for sheet_id in sheet_ids:
    #         num_strands = sheet_df[sheet_df["id"] == sheet_id][
    #             "number_strands"
    #         ].values[0]
    #         for range_id in range_ids:
    #             sheet = Sheet()
    #             sheet.range_id = range_id
    #             sheet.sheet_id = sheet_id
    #             sheet.num_strands = num_strands
    #             range_row = range_df[
    #                 (range_df["sheet_id"] == sheet_id)
    #                 & (range_df["id"] == range_id)
    #             ]
    #             if len(range_row) == 0:
    #                 continue
    #             sheet.init_res_name = range_row["beg_auth_comp_id"].values[0]
    #             sheet.init_chain_id = range_row["beg_auth_asym_id"].values[0]
    #             sheet.init_seq_num = range_row["beg_auth_seq_id"].values[0]
    #             sheet.init_ins_code = range_row[
    #                 "pdbx_beg_PDB_ins_code"
    #             ].values[0]
    #             sheet.end_res_name = range_row["end_auth_comp_id"].values[0]
    #             sheet.end_chain_id = range_row["end_auth_asym_id"].values[0]
    #             sheet.end_seq_num = range_row["end_auth_seq_id"].values[0]
    #             sheet.end_ins_code = range_row[
    #                 "pdbx_end_PDB_ins_code"
    #             ].values[0]
    #             order_row = order_df[
    #                 (order_df["sheet_id"] == sheet_id)
    #                 & (order_df["range_id_1"] == range_id)
    #             ]
    #             if len(order_row) == 0:
    #                 order_row = order_df[
    #                     (order_df["sheet_id"] == sheet_id)
    #                     & (order_df["range_id_2"] == range_id)
    #                 ]
    #             sense = order_row["sense"].values[0]
    #             if sense == "parallel":
    #                 sheet.sense = 1
    #             elif sense == "anti-parallel":
    #                 sheet.sense = -1
    #             else:
    #                 raise NotImplementedError(sense)
    #             hbond_row = hbond_df[
    #                 (hbond_df["sheet_id"] == sheet_id)
    #                 & (hbond_df["range_id_1"] == range_id)
    #             ]
    #             if len(hbond_row) == 0:
    #                 hbond_row = hbond_df[
    #                     (hbond_df["sheet_id"] == sheet_id)
    #                     & (hbond_df["range_id_2"] == range_id)
    #                 ]
    #             sheet.curr_atom = hbond_row["range_2_auth_atom_id"].values[0]
    #             sheet.curr_res_name = hbond_row[
    #                 "range_2_auth_comp_id"
    #             ].values[0]
    #             sheet.curr_chain_id = hbond_row[
    #                 "range_2_auth_asym_id"
    #             ].values[0]
    #             sheet.curr_res_seq = hbond_row["range_2_auth_seq_id"].values[
    #                 0
    #             ]
    #             sheet.curr_ins_code = hbond_row[
    #                 "range_2_PDB_ins_code"
    #             ].values[0]
    #             sheet.prev_atom = hbond_row["range_1_auth_atom_id"].values[0]
    #             sheet.prev_res_name = hbond_row[
    #                 "range_1_auth_comp_id"
    #             ].values[0]
    #             sheet.prev_chain_id = hbond_row[
    #                 "range_1_auth_asym_id"
    #             ].values[0]
    #             sheet.prev_res_name = hbond_row["range_1_auth_seq_id"].values[
    #                 0
    #             ]
    #             sheet.prev_ins_code = hbond_row[
    #                 "range_1_PDB_ins_code"
    #             ].values[0]
    #             sheets.append(sheet)
    #     return sheets

    def parse_pdb(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_pdb(line)
        self.range_id = int(line[7:10].strip())
        self.sheet_id = line[11:14].strip()
        self.num_strands = int(line[14:16].strip())
        self.init_res_name = line[17:20].strip()
        self.init_chain_id = line[21].strip()
        self.init_seq_num = int(line[22:26].strip())
        self.init_ins_code = line[26].strip()
        self.end_res_name = line[28:31].strip()
        self.end_chain_id = line[32].strip()
        self.end_seq_num = int(line[33:37].strip())
        self.end_ins_code = line[37].strip()
        self.sense = int(line[38:40].strip())
        self.curr_atom = line[41:45].strip()
        self.curr_res_name = line[45:48].strip()
        try:
            self.curr_chain_id = line[49].strip()
            try:
                self.curr_res_seq = int(line[50:54].strip())
            except ValueError:
                self.curr_res_seq = ""
            self.curr_ins_code = line[54].strip()
            self.prev_atom = line[56:60].strip()
            self.prev_res_name = line[60:63].strip()
            self.prev_chain_id = line[64].strip()
            try:
                self.prev_res_seq = int(line[65:69].strip())
            except ValueError:
                self.prev_res_seq = ""
            self.prev_ins_code = line[69].strip()
        except IndexError:
            pass

    def __str__(self):
        string = (
            f"SHEET  {self.range_id:3} {self.sheet_id:>3}{self.num_strands:>2}"
            f" {self.init_res_name:3} {self.init_chain_id:1}"
            f"{self.init_seq_num:>4}{self.init_ins_code:1} {self.end_res_name:3}"
            f" {self.end_chain_id:1}{self.end_seq_num:>4}{self.end_ins_code:1}"
            f"{self.sense:2}"
        )
        if len(self.curr_atom) == 1:
            string += f"  {self.curr_atom:3}"
        else:
            string += f" {self.curr_atom:4}"
        string += f"{self.curr_res_name:3} "
        string += f"{self.curr_chain_id:1}"
        string += f"{self.curr_res_seq:>4}"
        string += f"{self.curr_ins_code:1}"
        if len(self.prev_atom) == 1:
            string += f"  {self.prev_atom:3}"
        else:
            string += f" {self.prev_atom:4}"
        string += (
            f"{self.prev_res_name:3} {self.prev_chain_id:1}"
            f"{self.prev_res_seq:4}{self.prev_ins_code:1}"
        )
        return string
