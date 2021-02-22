"""Classes for PDB records that provide heterogen information.

.. codeauthor::  Todd Dolinsky
.. codeauthor::  Yong Huang
.. codeauthor::  Nathan Baker
"""
import logging
from collections import OrderedDict
from .general import BaseRecord


_LOGGER = logging.getLogger(__name__)


class Heterogen(BaseRecord):
    """HET field

    HET records are used to describe non-standard residues, such as
    prosthetic groups, inhibitors, solvent molecules, and ions for which
    coordinates are supplied. Groups are considered HET if they are:

    * not one of the standard amino acids, and
    * not one of the nucleic acids (C, G, A, T, U, and I), and
    * not one of the modified versions of nucleic acids (+C, +G, +A, +T, +U,
      and +I), and
    * not an unknown amino acid or nucleic acid where UNK is used to indicate
      the unknown residue name.

    Het records also describe heterogens for which the chemical identity is
    unknown, in which case the group is assigned the hetatm_id UNK.

    +---------+-------------+-------------+-----------------------------------+
    | COLUMNS | DATA TYPE   | FIELD       | DEFINITION                        |
    +=========+=============+=============+===================================+
    | 1-6     | Record name | "HET   "    |                                   |
    +---------+-------------+-------------+-----------------------------------+
    | 8-10    | LString(3)  | hetID       | Het identifier, right-justified.  |
    +---------+-------------+-------------+-----------------------------------+
    | 13      | Character   | ChainID     | Chain identifier.                 |
    +---------+-------------+-------------+-----------------------------------+
    | 14-17   | Integer     | seqNum      | Sequence number.                  |
    +---------+-------------+-------------+-----------------------------------+
    | 18      | AChar       | iCode       | Insertion code.                   |
    +---------+-------------+-------------+-----------------------------------+
    | 21-25   | Integer     | numHetAtoms | Number of HETATM records for the  |
    |         |             |             | group present in the entry.       |
    +---------+-------------+-------------+-----------------------------------+
    | 31-70   | String      | text        | Text describing Het group.        |
    +---------+-------------+-------------+-----------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.hetatm_id = None
        self.chain_id = None
        self.seq_num = None
        self.ins_code = None
        self.num_het_atoms = None
        self.text = None

    def parse_line(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_line(line)
        self.hetatm_id = line[7:10].strip()
        self.chain_id = line[12].strip()
        self.seq_num = int(line[13:17].strip())
        self.ins_code = line[17].strip()
        self.num_het_atoms = int(line[20:25].strip())
        self.text = line[30:70].strip()

    def __str__(self):
        return (
            f"HET    {self.hetatm_id:>3}  {self.chain_id:1}{self.seq_num:4}"
            f"{self.ins_code:1}  {self.num_het_atoms:5}     {self.text:40}"
        )


class HeterogenName(BaseRecord):
    """HETNAM field

    This record gives the chemical name of the compound with the
    given hetatm_id.

    +---------+--------------+--------------+---------------------------------+
    | COLUMNS | DATA TYPE    | FIELD        | DEFINITION                      |
    +=========+==============+==============+=================================+
    | 1-6     | Record name  | "HETNAM"     |                                 |
    +---------+--------------+--------------+---------------------------------+
    | 9-10    | Continuation | continuation | Allows concatenation of         |
    |         |              |              | multiple records.               |
    +---------+--------------+--------------+---------------------------------+
    | 12-14   | LString(3)   | hetID        | Het identifier, right-          |
    |         |              |              | justified.                      |
    +---------+--------------+--------------+---------------------------------+
    | 16-70   | String       | text         | Chemical name.                  |
    +---------+--------------+--------------+---------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.heterogens = OrderedDict()

    def parse_line(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_line(line)
        hetatm_id = line[11:14].strip()
        string = line[15:70].strip()
        strings = self.heterogens.get(hetatm_id, [])
        strings.append(string)
        self.heterogens[hetatm_id] = strings

    def __str__(self):
        strings = []
        for hetatm, lines in self.heterogens.items():
            for iline, line in enumerate(lines):
                continuation = iline + 1
                if continuation > 1:
                    string = (
                        f"HETNAM  {continuation:>2} {hetatm:>3}  {line:54}"
                    )
                else:
                    string = f"HETNAM     {hetatm:>3} {line:55}"
                strings.append(string)
        return "\n".join(strings)


class HeterogenSynonym(BaseRecord):
    """HETSYN field

    This record provides synonyms, if any, for the compound in the
    corresponding (i.e., same hetatm_id) HETNAM record. This is to allow
    greater flexibility in searching for HET groups.

    +----------+--------------+--------------+--------------------------------+
    | COLUMNS  | DATA TYPE    | FIELD        | DEFINITION                     |
    +==========+==============+==============+================================+
    | 1-6      | Record name  | "HETSYN"     |                                |
    +----------+--------------+--------------+--------------------------------+
    | 9-10     | Continuation | continuation | Allows concatenation of        |
    |          |              |              | multiple records.              |
    +----------+--------------+--------------+--------------------------------+
    | 12-14    | LString(3)   | hetID        | Het identifier, right-         |
    |          |              |              | justified.                     |
    +----------+--------------+--------------+--------------------------------+
    | 16-70    | SList        | hetSynonyms  | List of synonyms.              |
    +----------+--------------+--------------+--------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.synonyms = OrderedDict()

    def parse_line(self, line):
        super().parse_line(line)
        het_id = line[11:14].strip()
        synonyms = self.synonyms.get(het_id, [])
        synonyms.append(line[15:70].strip())
        self.synonyms[het_id] = synonyms

    def __str__(self):
        lines = []
        for het_id, synonyms in self.synonyms.items():
            for isyn, syn in enumerate(synonyms):
                continuation = isyn + 1
                if continuation > 1:
                    line = f"HETSYN  {continuation:>2} {het_id:3}  {syn:54}"
                else:
                    line = f"HETSYN     {het_id:3} {syn:55}"
                lines.append(line)
        return "\n".join(lines)


class Formula(BaseRecord):
    """FORMUL field

    The FORMUL record presents the chemical formula and charge of a
    non-standard group.

    +---------+-------------+--------------+----------------------------------+
    | COLUMNS | DATA TYPE   | FIELD        | DEFINITION                       |
    +=========+=============+==============+==================================+
    | 1-6     | Record name | "FORMUL"     |                                  |
    +---------+-------------+--------------+----------------------------------+
    | 9-10    | Integer     | compNum      | Component number.                |
    +---------+-------------+--------------+----------------------------------+
    | 13-15   | LString(3)  | hetID        | Het identifier.                  |
    +---------+-------------+--------------+----------------------------------+
    | 17-18   | Integer     | continuation | Continuation number.             |
    +---------+-------------+--------------+----------------------------------+
    | 19      | Character   | asterisk     | "*" for water.                   |
    +---------+-------------+--------------+----------------------------------+
    | 20-70   | String      | text         | Chemical formula.                |
    +---------+-------------+--------------+----------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.components = OrderedDict()

    def parse_line(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_line(line)
        component_num = int(line[8:10].strip())
        if component_num not in self.components:
            self.components[component_num] = []
        hetatm_id = line[12:15].strip()
        text = line[18:70].rstrip()
        self.components[component_num].append((hetatm_id, text))

    def __str__(self):
        strings = []
        for component_num, component_list in self.components.items():
            for hetatm_id, text in component_list:
                string = (
                    f"FORMUL  {component_num:>2}  {hetatm_id:>3}   {text:52}"
                ).strip()
                strings.append(string)
        return "\n".join(strings)
