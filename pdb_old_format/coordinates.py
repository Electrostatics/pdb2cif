"""Classes for records with coordinate information.

.. codeauthor::  Todd Dolinsky
.. codeauthor::  Yong Huang
.. codeauthor::  Nathan Baker
"""
import logging
from typing import OrderedDict
from .general import BaseRecord, atom_format

_LOGGER = logging.getLogger(__name__)


class Chain(BaseRecord):
    """Chain class.

    Represents PDB chain, containing ATOM, HETATM, ANISOU entries.
    """

    def __init__(self):
        super().__init__()
        self.chain_id = None
        self.atom = []
        self.temp_factor = []
        self.het_atom = []
        self.terminus = None
        self.has_atom = False
        self.has_het_atom = True

    def parse_line(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_line(line)
        name = line[0:6].strip()
        if name == "ATOM":
            self.has_atom = True
            atom = Atom()
            atom.parse_line(line)
            self.chain_id = atom.chain_id
            self.atom.append(atom)
        elif name == "HETATM":
            self.has_het_atom = True
            atom = Atom()
            atom.parse_line(line)
            self.chain_id = atom.chain_id
            self.het_atom.append(atom)
        elif name == "ANISOU":
            temp_factor = TemperatureFactor()
            temp_factor.parse_line(line)
            self.temp_factor.append(temp_factor)
        elif name == "TER":
            if not self.terminus:
                self.terminus = ChainTerminus()
                self.terminus.parse_line(line)
            else:
                err = (
                    f"Got new TER record: {line}\n"
                    f"Already have TER record: {self.terminus}"
                )
                raise ValueError(err)
        else:
            err = f"Unexpected line: {line}"
            raise ValueError(err)

    def num_residues(self, count_hetatm) -> int:
        """Number of residues in entry.

        :param bool count_hetam:  include heterogen residues in count
        """
        atom_res = {atom.res_seq for atom in self.atom}
        if count_hetatm:
            het_res = {atom.res_seq for atom in self.het_atom}
        else:
            het_res = set()
        return len(atom_res | het_res)

    def num_atoms(self, heavy_only) -> int:
        """Number of ATOM and HETATM entries in all chains.

        .. note::  PDB seems to only count heavy atoms, although doesn't
                   document this well.

        :param bool heavy_only:  exclude hydrogen atoms from count
        """
        num_atom = 0
        for atom in self.het_atom + self.atom:
            if not heavy_only:
                num_atom += 1
            elif atom.element not in ["H", "D"]:
                num_atom += 1
        return num_atom

    def __str__(self):
        raise NotImplementedError()


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
        self.chain = OrderedDict()

    def parse_line(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_line(line)
        name = line[0:6].strip()
        if name == "MODEL":
            self.serial = int(line[10:14].strip())
        elif name in ["ATOM", "ANISOU", "TER", "HETATM"]:
            chain_id = line[21].strip()
            chain = self.chain.get(chain_id, Chain())
            chain.parse_line(line)
            self.chain[chain_id] = chain
        elif name == "ENDMDL":
            pass
        else:
            err = f"Unexpected line: {line}"
            raise ValueError(err)

    def num_atoms(self, heavy_only) -> int:
        """Number of ATOM and HETATM entries in all chains in model.

        :param bool heavy_only:  exclude hydrogen atoms from count
        """
        num_atom = 0
        for chain in self.chain.values():
            num_atom += chain.num_atoms(heavy_only)
        return num_atom

    def num_chains(self) -> int:
        """Count number of chains in model."""
        return len(self.chain)

    def num_residues(self, count_hetatm) -> int:
        """Number of residues in entry.

        :param bool count_hetam:  include heterogen residues in count
        """
        num_residue = 0
        for chain in self.chain.values():
            num_residue += chain.num_residues(count_hetatm)
        return num_residue

    def __str__(self):
        strings = []
        if self.serial:
            strings.append(f"MODEL     {self.serial:4}".strip())
        for chain in self.chain.values():
            strings.append(chain.__str__())
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
        ).strip()


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
        string = f"ANISOU{self.serial:5}"
        if len(self.name) < 4:
            string += f"  {self.name:3}"
        else:
            string += f" {self.name:<4}"
        string += (
            f"{self.alt_loc:1}{self.res_name:3} {self.chain_id:1}"
            f"{self.res_seq:4}{self.ins_code:1} {self.u00:7}{self.u11:7}"
            f"{self.u22:7}{self.u01:7}{self.u02:7}{self.u12:7}"
            f"       {self.element:2}{self.charge:2}"
        )
        return string.strip()


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
        ).strip()


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
        ).strip()

