"""Top-level module for PDB structure entries.

The specifications used in this class are derived from the `Protein Data Bank
Contents Guide: Atomic Coordinate Entry Format Description, Version 3.3
<https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html>`_.

.. codeauthor::  Nathan Baker
"""
import logging
from . import annotation, primary, heterogen, secondary, coordinates
from . import crystallography, bookkeeping


_LOGGER = logging.getLogger(__name__)
REF_LINE = (
    "0        1         2         3         4         5         6         7         8\n"
    "12345678901234567890123456789012345678901234567890123456789012345678901234567890"
)


class Entry:
    """Top-level class for PDB structure entry."""

    def __init__(self):
        # Title section
        self.header = None
        self.obsolete = None
        self.title = None
        self.split = None
        self.caveat = None
        self.compound = None
        self.source = None
        self.keywords = None
        self.experimental_data = None
        self.num_models = 1
        self.model_type = None
        self.author = None
        self.revision_data = []
        self.supersedes = None
        self.journal = []
        self.remark = []
        # Primary structure section
        self.database_reference = []
        self.sequence_differences = []
        self.sequence_residues = None
        self.modified_residues = []
        # Heterogen section
        self.heterogen = []
        self.heterogen_name = None
        self.heterogen_synonyms = None
        self.heterogen_formula = None
        # Secondary sequence section
        self.helix = []
        self.sheet = []
        # Connectivity annotation section
        self.disulfide_bond = []
        self.link = []
        self.cis_peptide = []
        # Miscellaneous section
        self.site = []
        # Crystallographic and coordinate transformation section
        self.unit_cell = None
        self.orig_transform = []
        self.frac_transform = []
        self.noncrystal_transform = []
        # Coordinate section
        self.model = []
        # Connectivity section
        self.connect = []
        # Bookeeping section
        self.master = None
        self.end = None

    def parse_file(self, file_):
        """Parse a PDB file.

        :param file file_:  file open for reading.
        """
        for line in file_:
            try:
                self.parse_line(line)
            except Exception as exc:
                err = f"Offending line:\n{REF_LINE}\n{line}"
                raise ValueError(err) from exc

    def parse_line(self, line):
        """Parse a line of a PDB file.

        :param str line:  line of PDB file
        """
        name = line[0:6].strip()
        if name == "HEADER":
            if self.header:
                err = f"HEADER already exists:\n{self.header}"
                raise ValueError(err)
            self.header = annotation.Header()
            self.header.parse_line(line)
        elif name == "OBSLTE":
            if not self.obsolete:
                self.obsolete = annotation.Obsolete()
            self.obsolete.parse_line(line)
        elif name == "TITLE":
            if not self.title:
                self.title = annotation.Title()
            self.title.parse_line(line)
        elif name == "SPLIT":
            if not self.split:
                self.split = annotation.Split()
            self.split.parse_line(line)
        elif name == "CAVEAT":
            if not self.caveat:
                self.caveat = annotation.Caveat()
            self.caveat.parse_line(line)
        elif name == "COMPND":
            if not self.compound:
                self.compound = annotation.Compound()
            self.compound.parse_line(line)
        elif name == "SOURCE":
            if not self.source:
                self.source = annotation.Source()
            self.source.parse_line(line)
        elif name == "KEYWDS":
            if not self.keywords:
                self.keywords = annotation.Keywords()
            self.keywords.parse_line(line)
        elif name == "EXPDTA":
            if not self.experimental_data:
                self.experimental_data = annotation.ExperimentalData()
            self.experimental_data.parse_line(line)
        elif name == "NUMMDL":
            if self.num_models > 1:
                err = f"NUMMDL already exists:\n{self.num_models}"
                raise ValueError(err)
            self.num_models = int(line[10:14])
        elif name == "MDLTYP":
            if not self.model_type:
                self.model_type = annotation.ModelType()
            self.model_type.parse_line(line)
        elif name == "AUTHOR":
            if not self.author:
                self.author = annotation.Author()
            self.author.parse_line(line)
        elif name == "REVDAT":
            revdat = annotation.RevisionData()
            revdat.parse_line(line)
            self.revision_data.append(revdat)
        elif name == "SPRSDE":
            if not self.supersedes:
                self.supersedes = annotation.Supersedes()
            self.supersedes.parse_line(line)
        elif name == "JRNL":
            journal = annotation.Journal()
            journal.parse_line(line)
            self.journal.append(journal)
        elif name == "REMARK":
            remark = annotation.Remark()
            remark.parse_line(line)
            self.remark.append(remark)
        elif name == "DBREF":
            database = primary.DatabaseReference()
            database.parse_line(line)
            self.database_reference.append(database)
        elif name == "DBREF1":
            database = primary.DatabaseReference1()
            database.parse_line(line)
            self.database_reference.append(database)
        elif name == "DBREF2":
            database = primary.DatabaseReference2()
            database.parse_line(line)
            self.database_reference.append(database)
        elif name == "SEQADV":
            seqadv = primary.SequenceDifferences()
            seqadv.parse_line(line)
            self.sequence_differences.append(seqadv)
        elif name == "SEQRES":
            if not self.sequence_residues:
                self.sequence_residues = primary.SequenceResidues()
            self.sequence_residues.parse_line(line)
        elif name == "HET":
            het = heterogen.Heterogen()
            het.parse_line(line)
            self.heterogen.append(het)
        elif name == "HETNAM":
            if not self.heterogen_name:
                self.heterogen_name = heterogen.HeterogenName()
            self.heterogen_name.parse_line(line)
        elif name == "HETSYN":
            if not self.heterogen_synonyms:
                self.heterogen_synonyms = heterogen.HeterogenSynonym()
            self.heterogen_synonyms.parse_line(line)
        elif name == "FORMUL":
            if not self.heterogen_formula:
                self.heterogen_formula = heterogen.Formula()
            self.heterogen_formula.parse_line(line)
        elif name == "HELIX":
            helix = secondary.Helix()
            helix.parse_line(line)
            self.helix.append(helix)
        elif name == "SHEET":
            sheet = secondary.Sheet()
            sheet.parse_line(line)
            self.sheet.append(sheet)
        elif name == "SSBOND":
            bond = secondary.DisulfideBond()
            bond.parse_line(line)
            self.disulfide_bond.append(bond)
        elif name == "LINK":
            link = secondary.Link()
            link.parse_line(line)
            self.link.append(link)
        elif name == "CISPEP":
            pep = secondary.CisPeptide()
            pep.parse_line(line)
            self.cis_peptide.append(pep)
        elif name == "SITE":
            site = annotation.Site()
            site.parse_line(line)
            self.site.append(site)
        elif name == "CRYST1":
            if self.unit_cell is not None:
                err = f"CRYST1 already exists:\n{self.unit_cell}"
                raise ValueError(err)
            self.unit_cell = crystallography.UnitCell()
            self.unit_cell.parse_line(line)
        elif name in ["ORIGX1", "ORIGX2", "ORIGX3"]:
            n = int(name[5])
            orig = crystallography.OriginalTransform(n)
            orig.parse_line(line)
            self.orig_transform.append(orig)
            if len(self.orig_transform) > 3:
                err = f"Too many ({len(self.orig_transform)}) transforms."
                raise ValueError(err)
        elif name in ["SCALE1", "SCALE2", "SCALE3"]:
            n = int(name[5])
            scale = crystallography.FractionalTransform(n)
            scale.parse_line(line)
            self.frac_transform.append(scale)
            if len(self.frac_transform) > 3:
                err = f"Too many ({len(self.frac_transform)}) transforms."
                raise ValueError(err)
        elif name in ["MTRIX1", "MTRIX2", "MTRIX3"]:
            n = int(name[5])
            matrix = crystallography.NoncrystalTransform(n)
            matrix.parse_line(line)
            self.noncrystal_transform.append(matrix)
            if len(self.noncrystal_transform) > 3:
                err = (
                    f"Too many ({len(self.noncrystal_transform)}) transforms."
                )
                raise ValueError(err)
        elif name == "MODEL":
            model = coordinates.Model()
            model.parse_line(line)
            self.model.append(model)
        elif name in ["ATOM", "ANISOU", "TER", "HETATM"]:
            if len(self.model) == 0:
                self.model = [coordinates.Model()]
            self.model[-1].parse_line(line)
        elif name == "CONECT":
            connect = bookkeeping.Connection()
            connect.parse_line(line)
            self.connect.append(connect)
        elif name == "MASTER":
            if self.master:
                err = f"MASTER record already exists. Got: {line}."
                raise ValueError(err)
            self.master = bookkeeping.Master()
            self.master.parse_line(line)
        elif name in ["ENDMDL", "END"]:
            pass
        else:
            err = f"Unexpected entry:\n{line}"
            raise ValueError(err)

    def num_transforms(self) -> int:
        """Return the number of optional transform records in entry.

        :returns:  number of ORGIXn + SCALEn + MTRIXn
        """
        return (
            len(self.orig_transform) + len(self.frac_transform)
            + len(self.noncrystal_transform)
        )

    def num_atoms(self, heavy_only=True) -> int:
        """Number of ATOM and HETATM entries in all chains in entry.

        :param bool heavy_only:  exclude hydrogen atoms from count
        """
        return self.model[0].num_atoms(heavy_only)

    def num_chains(self) -> int:
        """Number of chains in entry."""
        return self.model[0].num_chains()

    def num_residues(self, count_hetatm=False) -> int:
        """Number of residues in entry.

        :param bool count_hetam:  include heterogen residues in count
        """
        return self.model[0].num_residues(count_hetatm)

    def num_ter(self) -> int:
        """Number of TER records in entry."""
        num = 0
        for model in self.model:
            for chain in model.chain.values():
                if chain.terminus is not None:
                    num += 1
        return num

    def check_master(self):
        """Check the contents against internal bookkeeping records.

        :raises AssertionError:  if checks fail
        """
        master = self.master
        for field, expected, test in [
            ("model", self.num_models, len(self.model)),
            ("REMARK", master.num_remark, len(self.remark)),
            ("HETATM", master.num_het, len(self.heterogen)),
            ("HELIX", master.num_helix, len(self.helix)),
            ("SHEET", master.num_sheet, len(self.sheet)),
            ("SITE", master.num_site, len(self.site)),
            ("transform", master.num_xform, self.num_transforms()),
            ("coordinate", master.num_coord, self.num_atoms()),
            ("TER", master.num_ter, self.num_ter()),
            ("CONECT", master.num_conect, len(self.connect)),
        ]:
            test_str = (
                f"MASTER indicates {expected} {field} records; found {test}."
            )
            _LOGGER.debug(test_str)
            try:
                assert(expected == test)
            except AssertionError:
                err = (
                    f"{test_str}"
                    " However, the MASTER record is hard to interpret."
                )
                _LOGGER.warning(err)
