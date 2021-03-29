"""Classes for records with crystallographic information.

.. codeauthor::  Todd Dolinsky
.. codeauthor::  Yong Huang
.. codeauthor::  Nathan Baker
"""
import logging
from .general import BaseRecord, cif_df


_LOGGER = logging.getLogger(__name__)


class FractionalTransform(BaseRecord):
    """SCALEn baseclass

    The SCALEn (n = 1, 2, or 3) records present the transformation from the
    orthogonal coordinates as contained in the entry to fractional
    crystallographic coordinates. Non-standard coordinate systems should be
    explained in the remarks.

    +---------+-------------+------------------------+------------------------+
    | COLUMNS | DATA TYPE   | FIELD                  | DEFINITION             |
    +=========+=============+========================+========================+
    | 1 -  6  | Record name | "SCALEn" n=1,  2, or 3 |                        |
    +---------+-------------+------------------------+------------------------+
    | 11 - 20 | Real(10.6)  | sn1                    | Sn1                    |
    +---------+-------------+------------------------+------------------------+
    | 21 - 30 | Real(10.6)  | sn2                    | Sn2                    |
    +---------+-------------+------------------------+------------------------+
    | 31 - 40 | Real(10.6)  | sn3                    | Sn3                    |
    +---------+-------------+------------------------+------------------------+
    | 46 - 55 | Real(10.5)  | unif                   | Un                     |
    +---------+-------------+------------------------+------------------------+
    """

    def __init__(self, n):
        """Initialize with the n in SCALEn.

        :param int n:  n in SCALEn
        """
        super().__init__()
        self.n = n
        self.sn1 = None
        self.sn2 = None
        self.sn3 = None
        self.unif = None

    @staticmethod
    def parse_cif(container) -> list:
        """Parse CIF container for information about this record.

        :param :class:`pdbx.containers.DataContainer` container:  container to
            parse
        :returns:  list of objects of this class
        """
        transforms = []
        df = cif_df(container.get_object("atom_sites"))
        for n in [1, 2, 3]:
            transform = FractionalTransform(n)
            transform.sn1 = float(df[f"fract_transf_matrix[{n}][1]"])
            transform.sn2 = float(df[f"fract_transf_matrix[{n}][2]"])
            transform.sn3 = float(df[f"fract_transf_matrix[{n}][3]"])
            transform.unif = float(df[f"fract_transf_vector[{n}]"])
            transforms.append(transform)
        return transforms

    def parse_pdb(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_pdb(line)
        self.sn1 = float(line[10:20].strip())
        self.sn2 = float(line[20:30].strip())
        self.sn3 = float(line[30:40].strip())
        self.unif = float(line[45:55].strip())

    def __str__(self):
        return (
            f"SCALE{self.n}    {self.sn1:10.6f}{self.sn2:10.6f}"
            f"{self.sn3:10.6f}     {self.unif:10.5f}"
        )


class OriginalTransform(BaseRecord):
    """ORIGXn class

    The ORIGXn (n = 1, 2, or 3) records present the transformation from the
    orthogonal coordinates contained in the entry to the submitted
    coordinates.

    +---------+-------------+----------+--------------------------------------+
    | COLUMNS | DATA TYPE   | FIELD    | DEFINITION                           |
    +=========+=============+==========+======================================+
    | 1-6     | Record name | "ORIGXn" | n=1, 2, or 3                         |
    +---------+-------------+----------+--------------------------------------+
    | 11-20   | Real(10.6)  | on1      | On1                                  |
    +---------+-------------+----------+--------------------------------------+
    | 21-30   | Real(10.6)  | on2      | On2                                  |
    +---------+-------------+----------+--------------------------------------+
    | 31-40   | Real(10.6)  | on3      | On3                                  |
    +---------+-------------+----------+--------------------------------------+
    | 46-55   | Real(10.5)  | tn       | Tn                                   |
    +---------+-------------+----------+--------------------------------------+
    """

    def __init__(self, n):
        """Initialize with the n in ORIGXn.

        :param int n:  n in ORGIXn
        """
        super().__init__()
        self.n = n
        self.on1 = None
        self.on2 = None
        self.on3 = None
        self.tn = None

    @staticmethod
    def parse_cif(container) -> list:
        """Parse CIF container for information about this record.

        :param :class:`pdbx.containers.DataContainer` container:  container to
            parse
        :returns:  list of objects of this class
        """
        transforms = []
        df = cif_df(container.get_object("database_PDB_matrix"))
        for n in [1, 2, 3]:
            transform = OriginalTransform(n)
            transform.on1 = float(df[f"origx[{n}][1]"])
            transform.on2 = float(df[f"origx[{n}][2]"])
            transform.on3 = float(df[f"origx[{n}][3]"])
            transform.tn = float(df[f"origx_vector[{n}]"])
            transforms.append(transform)
        return transforms

    def parse_pdb(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_pdb(line)
        self.on1 = float(line[10:20].strip())
        self.on2 = float(line[20:30].strip())
        self.on3 = float(line[30:40].strip())
        self.tn = float(line[45:55].strip())

    def __str__(self):
        return (
            f"ORIGX{self.n:1}    {self.on1:10.6f}{self.on2:10.6f}"
            f"{self.on3:10.6f}     {self.tn:10.5f}"
        )


class NoncrystalTransform(BaseRecord):
    """MTRIXn baseclass

    The MTRIXn (n = 1, 2, or 3) records present transformations expressing
    non-crystallographic symmetry.

    +---------+-------------+----------+--------------------------------------+
    | COLUMNS | DATA TYPE   | FIELD    | DEFINITION                           |
    +=========+=============+==========+======================================+
    | 1-6     | Record name | "MTRIXn" | n=1, 2, or 3                         |
    +---------+-------------+----------+--------------------------------------+
    | 8-10    | Integer     | serial   | Serial number.                       |
    +---------+-------------+----------+--------------------------------------+
    | 11-20   | Real(10.6)  | mn1      | Mn1                                  |
    +---------+-------------+----------+--------------------------------------+
    | 21-30   | Real(10.6)  | mn2      | Mn2                                  |
    +---------+-------------+----------+--------------------------------------+
    | 31-40   | Real(10.6)  | mn3      | Mn3                                  |
    +---------+-------------+----------+--------------------------------------+
    | 46-55   | Real(10.5)  | vn       | Vn                                   |
    +---------+-------------+----------+--------------------------------------+
    | 60      | Integer     | i_given  | 1 if coordinates for the             |
    |         |             |          | representations which are            |
    |         |             |          | approximately related by the         |
    |         |             |          | transformations of the molecule are  |
    |         |             |          | contained in the entry. Otherwise,   |
    |         |             |          | blank.                               |
    +---------+-------------+----------+--------------------------------------+
    """

    def __init__(self, n):
        """Initialize with n in MTRIXn

        :param int n:  n in MTRIXn
        """
        super().__init__()
        self.n = n
        self.serial = None
        self.mn1 = None
        self.mn2 = None
        self.mn3 = None
        self.vecn = None
        self.i_given = None

    @staticmethod
    def parse_cif(container) -> list:
        """Parse CIF container for information about this record.

        :param :class:`pdbx.containers.DataContainer` container:  container to
            parse
        :returns:  list of objects of this class
        """
        transforms = []
        df = cif_df(container.get_object("struct_ncs_oper"))
        for n in [1, 2, 3]:
            transform = NoncrystalTransform(n)
            transform.serial = int(df["id"])
            transform.mn1 = float(df[f"matrix[{n}][1]"])
            transform.mn2 = float(df[f"matrix[{n}][2]"])
            transform.mn3 = float(df[f"matrix[{n}][3]"])
            transform.vecn = float(df[f"vector[{n}]"])
            if df["code"].values[0] == "given":
                transform.i_given = 1
            else:
                raise NotImplementedError(df["code"])
            transforms.append(transform)
        return transforms

    def parse_pdb(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_pdb(line)
        self.serial = int(line[7:10].strip())
        self.mn1 = float(line[10:20].strip())
        self.mn2 = float(line[20:30].strip())
        self.mn3 = float(line[30:40].strip())
        self.vecn = float(line[45:55].strip())
        try:
            self.i_given = int(line[59].strip())
        except (ValueError, IndexError):
            pass

    def __str__(self):
        return (
            f"MTRIX{self.n:1} {self.serial:3}{self.mn1:10.6f}{self.mn2:10.6f}"
            f"{self.mn3:10.6f}     {self.vecn:10.5f}    {self.i_given:1}"
        )


class UnitCell(BaseRecord):
    """CRYST1 class

    The CRYST1 record presents the unit cell parameters, space group, and Z
    value. If the structure was not determined by crystallographic means,
    CRYST1 simply defines a unit cube.

    +---------+-------------+----------+--------------------------------------+
    | COLUMNS | DATA TYPE   | FIELD    | DEFINITION                           |
    +=========+=============+==========+======================================+
    | 1-6     | Record name | "CRYST1" |                                      |
    +---------+-------------+----------+--------------------------------------+
    | 7-15    | Real(9.3)   | a        | a (Angstroms).                       |
    +---------+-------------+----------+--------------------------------------+
    | 16-24   | Real(9.3)   | b        | b (Angstroms).                       |
    +---------+-------------+----------+--------------------------------------+
    | 25-33   | Real(9.3)   | c        | c (Angstroms).                       |
    +---------+-------------+----------+--------------------------------------+
    | 34-40   | Real(7.2)   | alpha    | alpha (degrees).                     |
    +---------+-------------+----------+--------------------------------------+
    | 41-47   | Real(7.2)   | beta     | beta (degrees).                      |
    +---------+-------------+----------+--------------------------------------+
    | 48-54   | Real(7.2)   | gamma    | gamma (degrees).                     |
    +---------+-------------+----------+--------------------------------------+
    | 56-66   | LString     | sGroup   | Space group.                         |
    +---------+-------------+----------+--------------------------------------+
    | 67-70   | Integer     | z        | Z value.                             |
    +---------+-------------+----------+--------------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.a = None
        self.b = None
        self.c = None
        self.alpha = None
        self.beta = None
        self.gamma = None
        self.space_group = None
        self.z = ""

    def parse_cif(self, container) -> bool:
        """Parse CIF container for information about this record.

        :param :class:`pdbx.containers.DataContainer` container:  container to
            parse
        :returns:  True if useful information was extracted from container
        """
        cell_df = cif_df(container.get_object("cell"))
        sym_df = cif_df(container.get_object("symmetry"))
        self.a = float(cell_df["length_a"].values[0])
        self.b = float(cell_df["length_b"].values[0])
        self.c = float(cell_df["length_c"].values[0])
        self.alpha = float(cell_df["angle_alpha"].values[0])
        self.beta = float(cell_df["angle_beta"].values[0])
        self.gamma = float(cell_df["angle_gamma"].values[0])
        self.space_group = sym_df["space_group_name_H-M"].values[0]
        self.z = cell_df["Z_PDB"].values[0]
        return True

    def parse_pdb(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_pdb(line)
        self.a = float(line[6:15].strip())
        self.b = float(line[15:24].strip())
        self.c = float(line[24:33].strip())
        self.alpha = float(line[33:40].strip())
        self.beta = float(line[40:47].strip())
        self.gamma = float(line[47:54].strip())
        self.space_group = line[55:65].strip()
        try:
            self.z = int(line[66:70].strip())
        except ValueError:
            pass

    def __str__(self):
        return (
            f"CRYST1{self.a:9.3f}{self.b:9.3f}{self.c:9.3f}{self.alpha:7.2f}"
            f"{self.beta:7.2f}{self.gamma:7.2f} {self.space_group:11}"
            f"{self.z:4}"
        )
