"""Top-level module for PDB structure entries.

The specifications used in this class are derived from the `Protein Data Bank
Contents Guide: Atomic Coordinate Entry Format Description, Version 3.3
<https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html>`_.

.. codeauthor::  Nathan Baker
"""
from old_pdb.primary import DatabaseReference
from old_pdb.secondary import Helix
from old_pdb.heterogen import HeterogenSynonym
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
        self._header = None
        self._obsolete = None
        self._title = None
        self._split = None
        self._caveat = None
        self._compound = None
        self._source = None
        self._keyword = None
        self._experimental_data = None
        self._num_model = None
        self._model_type = None
        self._author = None
        self._revision_data = None
        self._supersedes = None
        self._journal = []
        self._remark = []
        # Primary structure section
        self._database_reference = []
        self._sequence_difference = []
        self._sequence_residue = None
        self._modified_residue = []
        # Heterogen section
        self._heterogen = []
        self._heterogen_name = None
        self._heterogen_synonym = None
        self._heterogen_formula = None
        # Secondary structure section
        self._helix = []
        self._sheet = []
        # Connectivity annotation section
        self._disulfide_bond = []
        self._link = []
        self._cis_peptide = []
        # Miscellaneous section
        self.site = []
        # Crystallographic and coordinate transformation section
        self._unit_cell = None
        self._orig_transform = []
        self._frac_transform = []
        self._noncrystal_transform = []
        # Coordinate section
        self._model = []
        # Connectivity section
        self._connect = []
        # Bookeeping section
        self._master = None

    @property
    def header(self) -> annotation.Header:
        """:class:`.annotation.Header` HEADER record."""
        if self._header is None:
            self._header = annotation.Header()
        return self._header

    @header.setter
    def header(self, value):
        self._header = value

    @property
    def obsolete(self) -> annotation.Obsolete:
        """:class:`.annotation.Obsolete` OBSLTE record."""
        if self._obsolete is None:
            self._obsolete = annotation.Obsolete()
        return self._obsolete

    @obsolete.setter
    def obsolete(self, value):
        self._obsolete = value

    @property
    def title(self) -> annotation.Title:
        """:class:`.annotation.Title` TITLE record."""
        if self._title is None:
            self._title = annotation.Title()
        return self._title

    @title.setter
    def title(self, value):
        self._title = value

    @property
    def split(self) -> annotation.Split:
        """:class:`.annotation.Split` SPLIT record."""
        if self._split is None:
            self._split = annotation.Split()
        return self._split

    @split.setter
    def split(self, value):
        self._split = value

    @property
    def caveat(self) -> annotation.Caveat:
        """:class:`.annotation.Caveat` CAVEAT record."""
        if self._caveat is None:
            self._caveat = annotation.Caveat()
        return self._caveat

    @caveat.setter
    def caveat(self, value):
        self._caveat = value

    @property
    def compound(self) -> annotation.Compound:
        """:class:`.annotation.Compound` COMPND record."""
        if self._compound is None:
            self._compound = annotation.Compound()
        return self._compound

    @compound.setter
    def compound(self, value):
        self._compound = value

    @property
    def source(self) -> annotation.Source:
        """:class:`.annotation.Source` SOURCE record."""
        if self._source is None:
            self._source = annotation.Source()
        return self._source

    @source.setter
    def source(self, value):
        self._source = value

    @property
    def keyword(self) -> annotation.Keywords:
        """:class:`.annotation.Keywords` KEYWDS record."""
        if self._keyword is None:
            self._keyword = annotation.Keywords()
        return self._keyword

    @keyword.setter
    def keyword(self, value):
        self._keyword = value

    @property
    def experimental_data(self) -> annotation.ExperimentalData:
        """:class:`.annotation.ExperimentalData` EXPDTA record."""
        if self._experimental_data is None:
            self._experimental_data = annotation.ExperimentalData()
        return self._experimental_data

    @experimental_data.setter
    def experimental_data(self, value):
        self._experimental_data = value

    @property
    def num_model(self) -> annotation.NumModels:
        """:class:`.annotation.NumModels` NUMMDL record."""
        if self._num_model is None:
            _LOGGER.error("Creating new NumModel")
            self._num_model = annotation.NumModels()
        return self._num_model

    @num_model.setter
    def num_model(self, value):
        self._num_model = value

    @property
    def model_type(self) -> annotation.ModelType:
        """:class:`.annotation.ModelType` MDLTYP record."""
        if self._model_type is None:
            self._model_type = annotation.ModelType()
        return self._model_type

    @model_type.setter
    def model_type(self, value):
        self._model_type = value

    @property
    def author(self) -> annotation.Author:
        """:class:`.annotation.Author` AUTHOR record."""
        if self._author is None:
            self._author = annotation.Author()
        return self._author

    @author.setter
    def setter(self, value):
        self._setter = value

    @property
    def revision_data(self) -> annotation.RevisionData:
        """:class:`.annotation.RevisionData` REVDAT record."""
        if self._revision_data is None:
            self._revision_data = annotation.RevisionData()
        return self._revision_data

    @revision_data.setter
    def revision_data(self, value):
        self._revision_data = value

    @property
    def supersedes(self) -> annotation.Supersedes:
        """:class:`.annotation.Supersedes` SPRSDE record."""
        if self._supersedes is None:
            self._supersedes = annotation.Supersedes()
        return self._supersedes

    @supersedes.setter
    def supersedes(self, value):
        self._supersedes = value

    @property
    def journal(self) -> annotation.Journal:
        """:class:`.annotation.Journal` JRNL record."""
        if self._journal is None:
            self._journal = annotation.Journal()
        return self._journal

    @journal.setter
    def journal(self, value):
        self._journal = value

    @property
    def remark(self) -> list:
        """List of :class:`annotation.Remark` REMARK records."""
        return self._remark

    @remark.setter
    def remark(self, value):
        self._remark = value

    @property
    def database_reference(self) -> list:
        """List of :class:`primary.DatabaseReference` DBREF records."""
        return self._database_reference

    @database_reference.setter
    def database_reference(self, value):
        self._database_reference = value

    @property
    def sequence_difference(self) -> list:
        """List of :class:`primary.SequenceDifferences` SEQADV records."""
        return self._sequence_difference

    @sequence_difference.setter
    def sequence_difference(self, value):
        self._sequence_difference = value

    @property
    def sequence_residue(self) -> list:
        """List of :class:`.primary.SequenceResidues` SEQRES records."""
        return self._sequence_residue

    @sequence_residue.setter
    def sequence_residue(self, value):
        self._sequence_residue = value

    @property
    def modified_residue(self) -> list:
        """List of :class:`.primary.ModifiedResidue` MODRES records."""
        return self._modified_residue

    @modified_residue.setter
    def modified_residue(self, value):
        self._modified_residue = value

    @property
    def heterogen(self) -> list:
        """List of :class:`.heterogen.Heterogen` HET records."""
        return self._heterogen

    @heterogen.setter
    def heterogen(self, value):
        self._heterogen = value

    @property
    def heterogen_name(self) -> list:
        """List of :class:`.heterogen.HeterogenName` HETNAM records."""
        if self._heterogen_name is None:
            self._heterogen_name = heterogen.HeterogenName()
        return self._heterogen_name

    @heterogen_name.setter
    def heterogen_name(self, value):
        self._heterogen_name = value

    @property
    def heterogen_synonym(self):
        """:class:`.heterogen.HeterogenSynonym` HETSYN record."""
        if self._heterogen_synonym is None:
            self._heterogen_synonym = heterogen.HeterogenSynonym()
        return self._heterogen_synonym

    @heterogen_synonym.setter
    def heterogen_synonym(self, value):
        self._heterogen_synonym = value

    @property
    def heterogen_formula(self):
        """:class:`heterogen.Formula` FORMUL record."""
        if self._heterogen_formula is None:
            self._heterogen_formula = heterogen.Formula()

    @heterogen_formula.setter
    def heterogen_formula(self, value):
        self._heterogen_formula = value

    @property
    def helix(self) -> list:
        """List of :class:`.secondary.Helix` HELIX records."""
        return self._helix

    @helix.setter
    def helix(self, value):
        self._helix = value

    @property
    def sheet(self) -> list:
        """List of :class:`.secondary.Sheet` SHEET records."""
        return self._sheet

    @sheet.setter
    def sheet(self, value):
        self._sheet = value

    @property
    def disulfide_bond(self) -> list:
        """List of :class:`.secondary.DisulfideBond` SSBOND records."""
        return self._disulfide_bond

    @disulfide_bond.setter
    def disulfide_bond(self, value):
        self._disulfide_bond = value

    @property
    def link(self) -> list:
        """List of :class:`.secondary.Link` LINK records."""
        return self._link

    @link.setter
    def link(self, value):
        self._link = value

    @property
    def cis_peptide(self) -> list:
        """List of :class:`.secondary.CisPeptide` CISPEP records."""
        return self._cis_peptide

    @cis_peptide.setter
    def cis_peptide(self, value):
        self._cis_peptide = value

    @property
    def unit_cell(self) -> crystallography.UnitCell:
        """:class:`.crystallography.UnitCell` CRYST1 record."""
        if self._unit_cell is None:
            self._unit_cell = crystallography.UnitCell()
        return self._unit_cell

    @unit_cell.setter
    def unit_cell(self, value):
        self._unit_cell = value

    @property
    def original_transform(self) -> list:
        """List of :class:`.crystallography.OriginalTransform` ORIGX
        records."""
        return self._orig_transform

    @original_transform.setter
    def original_transform(self, value):
        self._orig_transform = value

    @property
    def frac_transform(self) -> list:
        """List of :class:`crystallography.FractionalTransform` SCALEn
        records."""
        return self._frac_transform

    @frac_transform.setter
    def frac_transform(self, value):
        self._frac_transform = value

    @property
    def noncrystal_transform(self) -> list:
        """List of :class:`.crystallography.NoncrystalTransform` MTRIXn
        records."""
        return self._noncrystal_transform

    @noncrystal_transform.setter
    def noncrystal_transform(self, value):
        self._noncrystal_transform = value

    @property
    def model(self) -> list:
        """List of :class:`.coordinates.Model` MODEL records."""
        return self._model

    @model.setter
    def model(self, value):
        self._model = value

    @property
    def connect(self) -> list:
        """List of :class:`.bookkeeping.Connection` CONECT records."""
        return self._connect

    @connect.setter
    def connect(self, value):
        self._connect = value

    @property
    def master(self) -> bookkeeping.Master:
        """:class:`.bookkeeping.Master` MASTER record."""
        if self._master is not None:
            self._master = bookkeeping.Master()
        return self._master

    def find_residue(self, chain_id, residue_id, model_num=1) -> list:
        """Find a specific residue.

        :param str chain_id:  chain ID to find
        :param int residue_id: residue ID to find
        :param int model_num:  model number to use
        :returns:  list of :class:`.coordinates.Atom`-like objects
        """
        model = self._model[model_num - 1]
        atoms = []
        for atom in model.all_atoms:
            if (atom.chain_id == chain_id) and (atom.res_seq == residue_id):
                atoms.append(atom)
        if len(atoms) == 0:
            _LOGGER.warning(
                f"Unable to find residue {residue_id} in chain {chain_id} of "
                f"model {model_num}."
            )
        return atoms

    def find_atom_by_name(
        self, chain_id, residue_id, atom_name, model_num=1
    ) -> coordinates.Atom:
        """Find a specific atom by name.

        :param str chain_id:  chain ID to find
        :param int residue_id: residue ID to find
        :param str atom_name:  name of atom to find
        :param int model_num:  model number to use
        :returns:  ATOM or HETATM object
        """
        atoms = self.find_residue(chain_id, residue_id, model_num)
        for atom in atoms:
            if atom.name == atom_name:
                return atom

    def annotate_link(self, record) -> secondary.Link:
        """Annotate LINK to indicate whether the named atoms are elements.

        Creates two new Boolean attributes in record:  ``is_element1`` and
        ``is_element2``.

        :param secondary.Link record:  record to annotate
        :returns:  annotated record
        """
        atom1 = self.find_atom_by_name(
            record.chain_id1, record.res_seq1, record.name1
        )
        try:
            record.is_element1 = atom1.name[:2] == atom1.element[:2]
        except AttributeError as exc:
            err = (
                f"Unable to find atom {record.name1} in residue "
                f"{record.res_seq1} in chain {record.chain_id1}."
            )
            raise AttributeError(err) from exc
        atom2 = self.find_atom_by_name(
            record.chain_id2, record.res_seq2, record.name2
        )
        record.is_element2 = atom2.name[:2] == atom2.element[:2]
        return record

    def __str__(self):
        strings = []
        # Title section
        for record in (
            [
                self._header,
                self._obsolete,
                self._title,
                self._split,
                self._caveat,
                self._compound,
                self._source,
                self._keyword,
                self._experimental_data,
                self._num_model,
                self._model_type,
                self.author,
            ]
            + [self._revision_data]
            + [self._supersedes]
            + self._journal
            + self._remark
        ):
            if record is not None:
                strings.append(str(record))
        # Primary structure section
        for record in (
            self._database_reference
            + self._sequence_difference
            + [self._sequence_residue]
            + self._modified_residue
        ):
            if record is not None:
                strings.append(str(record))
        # Heterogen section
        for record in self._heterogen + [
            self._heterogen_name,
            self._heterogen_synonym,
            self._heterogen_formula,
        ]:
            if record is not None:
                strings.append(str(record))
        # Secondary structure section
        for record in self._helix + self._sheet:
            if record is not None:
                strings.append(str(record))
        # Connectivity annotation section
        for record in self._disulfide_bond:
            if record is not None:
                strings.append(str(record))
        for record in self._link:
            record = self.annotate_link(record)
            if record is not None:
                strings.append(str(record))
        for record in self._cis_peptide:
            if record is not None:
                strings.append(str(record))
        # Miscellaneous section
        for record in self.site:
            if record is not None:
                strings.append(str(record))
        # Crystallographic and coordinate transformation section
        for record in (
            [self._unit_cell]
            + self._orig_transform
            + self._frac_transform
            + self._noncrystal_transform
        ):
            if record is not None:
                strings.append(str(record))
        # Coordinate section
        for record in self._model:
            if record is not None:
                strings.append(str(record))
                if len(self._model) > 1:
                    strings.append("ENDMDL")
        # Connectivity section
        for record in self._connect:
            if record is not None:
                strings.append(str(record))
        # Bookkeeping section
        if self._master is not None:
            strings.append(str(self._master))
        strings.append("END   ")
        return "\n".join(strings)

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
            if self._header:
                err = f"HEADER already exists:\n{self._header}"
                raise ValueError(err)
            self._header = annotation.Header()
            self._header.parse_line(line)
        elif name == "OBSLTE":
            if not self._obsolete:
                self._obsolete = annotation.Obsolete()
            self._obsolete.parse_line(line)
        elif name == "TITLE":
            if not self._title:
                self._title = annotation.Title()
            self._title.parse_line(line)
        elif name == "SPLIT":
            if not self._split:
                self._split = annotation.Split()
            self._split.parse_line(line)
        elif name == "CAVEAT":
            if not self._caveat:
                self._caveat = annotation.Caveat()
            self._caveat.parse_line(line)
        elif name == "COMPND":
            if not self._compound:
                self._compound = annotation.Compound()
            self._compound.parse_line(line)
        elif name == "SOURCE":
            if not self._source:
                self._source = annotation.Source()
            self._source.parse_line(line)
        elif name == "KEYWDS":
            if not self._keyword:
                self._keyword = annotation.Keywords()
            self._keyword.parse_line(line)
        elif name == "EXPDTA":
            if not self._experimental_data:
                self._experimental_data = annotation.ExperimentalData()
            self._experimental_data.parse_line(line)
        elif name == "NUMMDL":
            if self._num_model is not None:
                err = f"NUMMDL already exists:\n{self._num_model}"
                raise ValueError(err)
            self._num_model = annotation.NumModels()
            self._num_model.parse_line(line)
        elif name == "MDLTYP":
            if not self._model_type:
                self._model_type = annotation.ModelType()
            self._model_type.parse_line(line)
        elif name == "AUTHOR":
            if not self.author:
                self.author = annotation.Author()
            self.author.parse_line(line)
        elif name == "REVDAT":
            if not self._revision_data:
                self._revision_data = annotation.RevisionData()
            self._revision_data.parse_line(line)
        elif name == "SPRSDE":
            if not self._supersedes:
                self._supersedes = annotation.Supersedes()
            self._supersedes.parse_line(line)
        elif name == "JRNL":
            journal = annotation.Journal()
            journal.parse_line(line)
            self._journal.append(journal)
        elif name == "REMARK":
            remark = annotation.Remark()
            remark.parse_line(line)
            self._remark.append(remark)
        elif name == "DBREF":
            database = primary.DatabaseReference()
            database.parse_line(line)
            self._database_reference.append(database)
        elif name == "DBREF1":
            database = primary.DatabaseReference1()
            database.parse_line(line)
            self._database_reference.append(database)
        elif name == "DBREF2":
            database = primary.DatabaseReference2()
            database.parse_line(line)
            self._database_reference.append(database)
        elif name == "SEQADV":
            seqadv = primary.SequenceDifferences()
            seqadv.parse_line(line)
            self._sequence_difference.append(seqadv)
        elif name == "SEQRES":
            if not self._sequence_residue:
                self._sequence_residue = primary.SequenceResidues()
            self._sequence_residue.parse_line(line)
        elif name == "HET":
            het = heterogen.Heterogen()
            het.parse_line(line)
            self._heterogen.append(het)
        elif name == "HETNAM":
            if not self._heterogen_name:
                self._heterogen_name = heterogen.HeterogenName()
            self._heterogen_name.parse_line(line)
        elif name == "HETSYN":
            if not self._heterogen_synonym:
                self._heterogen_synonym = heterogen.HeterogenSynonym()
            self._heterogen_synonym.parse_line(line)
        elif name == "FORMUL":
            if not self._heterogen_formula:
                self._heterogen_formula = heterogen.Formula()
            self._heterogen_formula.parse_line(line)
        elif name == "HELIX":
            helix = secondary.Helix()
            helix.parse_line(line)
            self._helix.append(helix)
        elif name == "SHEET":
            sheet = secondary.Sheet()
            sheet.parse_line(line)
            self._sheet.append(sheet)
        elif name == "SSBOND":
            bond = secondary.DisulfideBond()
            bond.parse_line(line)
            self._disulfide_bond.append(bond)
        elif name == "LINK":
            link = secondary.Link()
            link.parse_line(line)
            self._link.append(link)
        elif name == "CISPEP":
            pep = secondary.CisPeptide()
            pep.parse_line(line)
            self._cis_peptide.append(pep)
        elif name == "SITE":
            site = annotation.Site()
            site.parse_line(line)
            self.site.append(site)
        elif name == "CRYST1":
            if self._unit_cell is not None:
                err = f"CRYST1 already exists:\n{self._unit_cell}"
                raise ValueError(err)
            self._unit_cell = crystallography.UnitCell()
            self._unit_cell.parse_line(line)
        elif name in ["ORIGX1", "ORIGX2", "ORIGX3"]:
            n = int(name[5])
            orig = crystallography.OriginalTransform(n)
            orig.parse_line(line)
            self._orig_transform.append(orig)
            if len(self._orig_transform) > 3:
                err = f"Too many ({len(self._orig_transform)}) transforms."
                raise ValueError(err)
        elif name in ["SCALE1", "SCALE2", "SCALE3"]:
            n = int(name[5])
            scale = crystallography.FractionalTransform(n)
            scale.parse_line(line)
            self._frac_transform.append(scale)
            if len(self._frac_transform) > 3:
                err = f"Too many ({len(self._frac_transform)}) transforms."
                raise ValueError(err)
        elif name in ["MTRIX1", "MTRIX2", "MTRIX3"]:
            n = int(name[5])
            matrix = crystallography.NoncrystalTransform(n)
            matrix.parse_line(line)
            self._noncrystal_transform.append(matrix)
            if len(self._noncrystal_transform) > 3:
                err = (
                    f"Too many ({len(self._noncrystal_transform)}) transforms."
                )
                raise ValueError(err)
        elif name == "MODEL":
            model = coordinates.Model()
            model.parse_line(line)
            self._model.append(model)
        elif name in ["ATOM", "ANISOU", "TER", "HETATM"]:
            if len(self._model) == 0:
                self._model = [coordinates.Model()]
            self._model[-1].parse_line(line)
        elif name == "CONECT":
            connect = bookkeeping.Connection()
            connect.parse_line(line)
            self._connect.append(connect)
        elif name == "MASTER":
            if self._master:
                err = f"MASTER record already exists. Got: {line}."
                raise ValueError(err)
            self._master = bookkeeping.Master()
            self._master.parse_line(line)
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
            len(self._orig_transform)
            + len(self._frac_transform)
            + len(self._noncrystal_transform)
        )

    def num_atoms(self, heavy_only=True) -> int:
        """Number of ATOM and HETATM entries in all chains in entry.

        :param bool heavy_only:  exclude hydrogen atoms from count
        """
        return self._model[0].num_atoms(heavy_only)

    def num_chains(self) -> int:
        """Number of chains in entry."""
        return self._model[0].num_chains()

    def num_residues(self, count_hetatm=False) -> int:
        """Number of residues in entry.

        :param bool count_hetam:  include heterogen residues in count
        """
        return self._model[0].num_residues(count_hetatm)

    def num_ter(self) -> int:
        """Number of TER records in entry."""
        num = 0
        for model in self._model:
            num += model.num_ter()
        return num

    def check_master(self):
        """Check the contents against internal bookkeeping records.

        :raises AssertionError:  if checks fail
        """
        master = self._master
        for field, expected, test in [
            (
                "model",
                self._num_model.model_number
                if self._num_model is not None
                else 1,
                len(self._model),
            ),
            ("REMARK", master.num_remark, len(self._remark)),
            ("HETATM", master.num_het, len(self._heterogen)),
            ("HELIX", master.num_helix, len(self._helix)),
            ("SHEET", master.num_sheet, len(self._sheet)),
            ("SITE", master.num_site, len(self.site)),
            ("transform", master.num_xform, self.num_transforms()),
            ("coordinate", master.num_coord, self.num_atoms()),
            ("TER", master.num_ter, self.num_ter()),
            ("CONECT", master.num_conect, len(self._connect)),
        ]:
            test_str = (
                f"MASTER indicates {expected} {field} records; found {test}."
            )
            try:
                assert expected == test
            except AssertionError:
                err = (
                    f"{test_str}"
                    " However, the MASTER record is hard to interpret."
                )
                _LOGGER.warning(err)
