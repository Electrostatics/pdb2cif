"""Classes for records with coordinate information.

.. codeauthor::  Todd Dolinsky
.. codeauthor::  Yong Huang
.. codeauthor::  Nathan Baker
"""
from itertools import count
import logging
from typing import OrderedDict
from .general import BaseRecord, atom_format

_LOGGER = logging.getLogger(__name__)


class Model(BaseRecord):
    """MODEL class.

    The MODEL record specifies the model serial number when multiple
    structures are presented in a single coordinate entry, as is often the
    case with structures determined by NMR.

    +---------+-------------+----------+--------------------------------------+
    | COLUMNS | DATA TYPE   | FIELD    | DEFINITION                           |
    +=========+=============+==========+======================================+
    | 1-6     | Record name | "MODEL " |                                      |
    +---------+-------------+----------+--------------------------------------+
    | 11-14   | Integer     | serial   | Model serial number.                 |
    +---------+-------------+----------+--------------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.serial = None
        self.records = []

    def parse_line(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_line(line)
        name = line[0:6].strip()
        if name == "MODEL":
            self.serial = int(line[10:14].strip())
            return
        if name == "ENDMDL":
            return
        if name == "ATOM":
            record = Atom()
        elif name == "ANISOU":
            record = TemperatureFactor()
        elif name == "TER":
            record = ChainTerminus()
        elif name == "HETATM":
            record = HeterogenAtom()
        else:
            err = f"Unexpected line: {line}"
            raise ValueError(err)
        record.parse_line(line)
        self.records.append(record)

    @property
    def all_atoms(self) -> list:
        """Get all atoms in model.

        :returns:  list of :class:`Atom`-like objects
        """
        return [
            rec for rec in self.records if isinstance(
                rec, (Atom, HeterogenAtom)
            )
        ]

    @property
    def het_atoms(self) -> list:
        """Get HETATM atoms in model.

        :returns:  list of :class:`Atom`-like objects
        """
        return [
            rec for rec in self.records if isinstance(rec, HeterogenAtom)
        ]

    @property
    def atoms(self) -> list:
        """Get ATOM atoms in model.

        :returns:  list of :class:`Atom`-like objects
        """
        return [
            rec for rec in self.records if isinstance(rec, Atom)
        ]

    def num_atoms(self, heavy_only) -> int:
        """Number of ATOM and HETATM entries in all chains in model.

        :param bool heavy_only:  exclude hydrogen atoms from count
        """
        num_atom = 0
        for atom in self.all_atoms:
            if (atom.element not in ["H", "D"]) or (not heavy_only):
                num_atom += 1
        return num_atom

    def num_chains(self) -> int:
        """Count number of chains in model."""
        chains = set()
        for atom in self.all_atoms:
            chains.add(atom.chain_id)
        return len(chains)

    def num_residues(self, count_hetatm) -> int:
        """Number of residues in entry.

        :param bool count_hetatm:  include heterogen residues in count
        """
        residues = set()
        if count_hetatm:
            atom_list = self.all_atoms
        else:
            atom_list = self.atoms
        for atom in atom_list:
            key = f"{atom.chain_id}{atom.res_name}{atom.res_seq}"
            residues.add(key)
        return len(residues)

    def num_ter(self) -> int:
        """Count number of termini in entry."""
        return len(
            [
                rec for rec in self.records if isinstance(rec, ChainTerminus)
            ]
        )

    def __str__(self):
        strings = []
        if self.serial:
            strings.append(f"MODEL     {self.serial:4}".strip())
        for record in self.records:
            strings.append(str(record))
        return "\n".join(strings)


class Atom(BaseRecord):
    """ATOM class

    The ATOM records present the atomic coordinates for standard residues.
    They also present the occupancy and temperature factor for each atom.
    Heterogen coordinates use the HETATM record type. The element symbol is
    always present on each ATOM record; segment identifier and charge are
    optional.

    +---------+--------------+------------+-----------------------------------+
    | COLUMNS | DATA TYPE    | FIELD      | DEFINITION                        |
    +=========+==============+============+===================================+
    | 1-6     | Record name  |  "ATOM  "  |                                   |
    +---------+--------------+------------+-----------------------------------+
    | 7-11    | Integer      | serial     | Atom serial number.               |
    +---------+--------------+------------+-----------------------------------+
    | 13-16   | Atom         | name       | Atom name.                        |
    +---------+--------------+------------+-----------------------------------+
    | 17      | Character    | altLoc     | Alternate location indicator.     |
    +---------+--------------+------------+-----------------------------------+
    | 18-20   | Residue name | resName    | Residue name.                     |
    +---------+--------------+------------+-----------------------------------+
    | 22      | Character    | chainID    | Chain identifier.                 |
    +---------+--------------+------------+-----------------------------------+
    | 23-26   | Integer      | resSeq     | Residue sequence number.          |
    +---------+--------------+------------+-----------------------------------+
    | 27      | AChar        | iCode      | Code for insertion of residues.   |
    +---------+--------------+------------+-----------------------------------+
    | 31-38   | Real(8.3)    | x          | Orthogonal coordinates for X in   |
    |         |              |            | Angstroms.                        |
    +---------+--------------+------------+-----------------------------------+
    | 39-46   | Real(8.3)    | y          | Orthogonal coordinates for Y in   |
    |         |              |            | Angstroms.                        |
    +---------+--------------+------------+-----------------------------------+
    | 47-54   | Real(8.3)    | z          | Orthogonal coordinates for Z in   |
    |         |              |            | Angstroms.                        |
    +---------+--------------+------------+-----------------------------------+
    | 55-60   | Real(6.2)    | occupancy  | Occupancy.                        |
    +---------+--------------+------------+-----------------------------------+
    | 61-66   | Real(6.2)    | tempFactor | Temperature factor.               |
    +---------+--------------+------------+-----------------------------------+
    | 77-78   | LString(2)   | element    | Element symbol, right-justified.  |
    +---------+--------------+------------+-----------------------------------+
    | 79-80   | LString(2)   | charge     | Charge on the atom.               |
    +---------+--------------+------------+-----------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.serial = None
        self.name = None
        self.alt_loc = None
        self.res_name = None
        self.chain_id = None
        self.res_seq = None
        self.ins_code = None
        self.x = None
        self.y = None
        self.z = None
        self.occupancy = 0.00
        self.temp_factor = 0.00
        self.seg_id = ""
        self.element = ""
        self.charge = ""

    def parse_line(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_line(line)
        self.serial = int(line[6:11].strip())
        self.name = line[12:16].strip()
        self.alt_loc = line[16].strip()
        self.res_name = line[17:20].strip()
        self.chain_id = line[21].strip()
        self.res_seq = int(line[22:26].strip())
        self.ins_code = line[26].strip()
        self.x = float(line[30:38].strip())
        self.y = float(line[38:46].strip())
        self.z = float(line[46:54].strip())
        try:
            self.occupancy = float(line[54:60].strip())
            self.temp_factor = float(line[60:66].strip())
            self.seg_id = line[72:76].strip()
            self.element = line[76:78].strip()
            self.charge = line[78:80].strip()
        except (ValueError, IndexError):
            pass

    def __str__(self):
        return (
            f"ATOM  {self.serial:5} {atom_format(self)}{self.alt_loc:1}"
            f"{self.res_name:>3} {self.chain_id:1}{self.res_seq:4}"
            f"{self.ins_code:1}   {self.x:8.3f}{self.y:8.3f}"
            f"{self.z:8.3f}{self.occupancy:6.2f}{self.temp_factor:6.2f}"
            f"           {self.element}{self.charge:2}"
        )


class TemperatureFactor(BaseRecord):
    """ANISOU class

    The ANISOU records present the anisotropic temperature factors.

    +---------+--------------+----------+-------------------------------------+
    | COLUMNS | DATA TYPE    | FIELD    | DEFINITION                          |
    +=========+==============+==========+=====================================+
    | 1-6     | Record name  | "ANISOU" |                                     |
    +---------+--------------+----------+-------------------------------------+
    | 7-11    | Integer      | serial   | Atom serial number.                 |
    +---------+--------------+----------+-------------------------------------+
    | 13-16   | Atom         | name     | Atom name.                          |
    +---------+--------------+----------+-------------------------------------+
    | 17      | Character    | altLoc   | Alternate location indicator        |
    +---------+--------------+----------+-------------------------------------+
    | 18-20   | Residue name | resName  | Residue name.                       |
    +---------+--------------+----------+-------------------------------------+
    | 22      | Character    | chainID  | Chain identifier.                   |
    +---------+--------------+----------+-------------------------------------+
    | 23-26   | Integer      | resSeq   | Residue sequence number.            |
    +---------+--------------+----------+-------------------------------------+
    | 27      | AChar        | iCode    | Insertion code.                     |
    +---------+--------------+----------+-------------------------------------+
    | 29-35   | Integer      | u[0][0]  | U(1,1)                              |
    +---------+--------------+----------+-------------------------------------+
    | 36-42   | Integer      | u[1][1]  | U(2,2)                              |
    +---------+--------------+----------+-------------------------------------+
    | 43-49   | Integer      | u[2][2]  | U(3,3)                              |
    +---------+--------------+----------+-------------------------------------+
    | 50-56   | Integer      | u[0][1]  | U(1,2)                              |
    +---------+--------------+----------+-------------------------------------+
    | 57-63   | Integer      | u[0][2]  | U(1,3)                              |
    +---------+--------------+----------+-------------------------------------+
    | 64-70   | Integer      | u[1][2]  | U(2,3)                              |
    +---------+--------------+----------+-------------------------------------+
    | 77-78   | LString(2)   | element  | Element symbol, right-justified.    |
    +---------+--------------+----------+-------------------------------------+
    | 79-80   | LString(2)   | charge   | Charge on the atom.                 |
    +---------+--------------+----------+-------------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.serial = None
        self.name = None
        self.alt_loc = None
        self.res_name = None
        self.chain_id = None
        self.res_seq = None
        self.ins_code = None
        self.u00 = None
        self.u11 = None
        self.u22 = None
        self.u01 = None
        self.u02 = None
        self.u12 = None
        self.seg_id = None
        self.element = None
        self.charge = None

    def parse_line(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_line(line)
        self.serial = int(line[6:11].strip())
        self.name = line[12:16].strip()
        self.alt_loc = line[16].strip()
        self.res_name = line[17:20].strip()
        self.chain_id = line[21].strip()
        self.res_seq = int(line[22:26].strip())
        self.ins_code = line[26].strip()
        self.u00 = int(line[28:35].strip())
        self.u11 = int(line[35:42].strip())
        self.u22 = int(line[42:49].strip())
        self.u01 = int(line[49:56].strip())
        self.u02 = int(line[56:63].strip())
        self.u12 = int(line[63:70].strip())
        self.seg_id = line[72:76].strip()
        self.element = line[76:78].strip()
        self.charge = line[78:80].strip()

    def __str__(self):
        return(
            f"ANISOU{self.serial:5} {atom_format(self)}{self.alt_loc:1}"
            f"{self.res_name:>3} {self.chain_id:1}{self.res_seq:4}"
            f"{self.ins_code:1} {self.u00:7}{self.u11:7}{self.u22:7}"
            f"{self.u01:7}{self.u02:7}{self.u12:7}      {self.element:>2}"
            f"{self.charge:2}"
        )


class ChainTerminus(BaseRecord):
    """TER class

    The TER record indicates the end of a list of ATOM/HETATM records for a
    chain.

    +---------+--------------+----------+-------------------------------------+
    | COLUMNS | DATA TYPE    | FIELD    | DEFINITION                          |
    +=========+==============+==========+=====================================+
    | 1-6     | Record name  | "TER   " |                                     |
    +---------+--------------+----------+-------------------------------------+
    | 7-11    | Integer      | serial   | Serial number.                      |
    +---------+--------------+----------+-------------------------------------+
    | 18-20   | Residue name | resName  | Residue name.                       |
    +---------+--------------+----------+-------------------------------------+
    | 22      | Character    | chainID  | Chain identifier.                   |
    +---------+--------------+----------+-------------------------------------+
    | 23-26   | Integer      | resSeq   | Residue sequence number.            |
    +---------+--------------+----------+-------------------------------------+
    | 27      | AChar        | iCode    | Insertion code.                     |
    +---------+--------------+----------+-------------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.serial = None
        self.res_name = None
        self.chain_id = None
        self.res_seq = None
        self.ins_code = ""

    def parse_line(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_line(line)
        if line is None:
            line = ""
        try:
            self.serial = int(line[6:11].strip())
            self.res_name = line[17:20].strip()
            self.chain_id = line[21].strip()
            self.res_seq = int(line[22:26].strip())
            self.ins_code = line[26].strip()
        except (IndexError, ValueError):
            pass

    def __str__(self):
        return (
            f"TER   {self.serial:5}      {self.res_name:>3}"
            f" {self.chain_id:1}{self.res_seq:4}{self.ins_code:1}"
        )


class HeterogenAtom(BaseRecord):
    """HETATM class

    The HETATM records present the atomic coordinate records for atoms
    within "non-standard" groups. These records are used for water
    molecules and atoms presented in HET groups.

    +---------+--------------+------------+-----------------------------------+
    | COLUMNS | DATA TYPE    | FIELD      | DEFINITION                        |
    +=========+==============+============+===================================+
    | 1-6     | Record name  | "HETATM"   |                                   |
    +---------+--------------+------------+-----------------------------------+
    | 7-11    | Integer      | serial     | Atom serial number.               |
    +---------+--------------+------------+-----------------------------------+
    | 13-16   | Atom         | name       | Atom name.                        |
    +---------+--------------+------------+-----------------------------------+
    | 17      | Character    | altLoc     | Alternate location indicator.     |
    +---------+--------------+------------+-----------------------------------+
    | 18-20   | Residue name | resName    | Residue name.                     |
    +---------+--------------+------------+-----------------------------------+
    | 22      | Character    | chainID    | Chain identifier.                 |
    +---------+--------------+------------+-----------------------------------+
    | 23-26   | Integer      | resSeq     | Residue sequence number.          |
    +---------+--------------+------------+-----------------------------------+
    | 27      | AChar        | iCode      | Code for insertion of residues.   |
    +---------+--------------+------------+-----------------------------------+
    | 31-38   | Real(8.3)    | x          | Orthogonal coordinates for X.     |
    +---------+--------------+------------+-----------------------------------+
    | 39-46   | Real(8.3)    | y          | Orthogonal coordinates for Y.     |
    +---------+--------------+------------+-----------------------------------+
    | 47-54   | Real(8.3)    | z          | Orthogonal coordinates for Z.     |
    +---------+--------------+------------+-----------------------------------+
    | 55-60   | Real(6.2)    | occupancy  | Occupancy.                        |
    +---------+--------------+------------+-----------------------------------+
    | 61-66   | Real(6.2)    | tempFactor | Temperature factor.               |
    +---------+--------------+------------+-----------------------------------+
    | 77-78   | LString(2)   | element    | Element symbol; right-justified.  |
    +---------+--------------+------------+-----------------------------------+
    | 79-80   | LString(2)   | charge     | Charge on the atom.               |
    +---------+--------------+------------+-----------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.serial = None
        self.name = None
        self.alt_loc = None
        self.res_name = None
        self.chain_id = None
        self.res_seq = None
        self.ins_code = None
        self.x = None
        self.y = None
        self.z = None
        self.occupancy = 0.00
        self.temp_factor = 0.00
        self.seg_id = ""
        self.element = ""
        self.charge = ""

    def parse_line(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_line(line)
        self.serial = int(line[6:11].strip())
        self.name = line[12:16].strip()
        self.alt_loc = line[16].strip()
        try:
            self.res_name = line[17:20].strip()
            self.chain_id = line[21].strip()
            self.res_seq = int(line[22:26].strip())
            self.ins_code = line[26].strip()
        except IndexError:
            raise ValueError("Residue name must be less than 4 characters!")
        self.x = float(line[30:38].strip())
        self.y = float(line[38:46].strip())
        self.z = float(line[46:54].strip())
        try:
            self.occupancy = float(line[54:60].strip())
            self.temp_factor = float(line[60:66].strip())
            self.seg_id = line[72:76].strip()
            self.element = line[76:78].strip()
            self.charge = line[78:80].strip()
        except (ValueError, IndexError):
            pass

    def __str__(self):
        return (
            f"HETATM{self.serial:5} {atom_format(self)}{self.alt_loc:1}"
            f"{self.res_name:>3} {self.chain_id:1}{self.res_seq:4}"
            f"{self.ins_code:1}   {self.x:8.3f}{self.y:8.3f}{self.z:8.3f}"
            f"{self.occupancy:6.2f}{self.temp_factor:6.2f}"
            f"          {self.element:>2}{self.charge:2}"
        )
