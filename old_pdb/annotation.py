"""Classes for PDB records that provide annotation information.

.. codeauthor::  Todd Dolinsky
.. codeauthor::  Yong Huang
.. codeauthor::  Nathan Baker
"""
import logging
from collections import OrderedDict
from .general import BaseRecord, grouper, date_parse, date_format


_LOGGER = logging.getLogger(__name__)


class Author(BaseRecord):
    """AUTHOR field

    The AUTHOR record contains the names of the people responsible for the
    contents of the entry.

    +---------+--------------+---------------+-------------------------------+
    | COLUMNS | DATA TYPE    | FIELD         | DEFINITION                    |
    +=========+==============+===============+===============================+
    | 1-6     | Record name  | "AUTHOR"      |                               |
    +---------+--------------+---------------+-------------------------------+
    | 9-10    | Continuation | continuation  | Allows concatenation of       |
    |         |              |               | multiple records.             |
    +---------+--------------+---------------+-------------------------------+
    | 11-79   | List         | author_list   | List of the author names,     |
    |         |              |               | separated by commas.          |
    +---------+--------------+---------------+-------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.author_list = []

    def parse_line(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_line(line)
        self.author_list.append(line[10:79].strip())

    def __str__(self):
        strings = []
        for iline, line in enumerate(self.author_list):
            continuation = iline + 1
            if continuation > 1:
                strings += [f"AUTHOR  {continuation:>2} {line:78}"]
            else:
                strings += [f"AUTHOR    {line:79}"]
        return "\n".join(strings)


class Caveat(BaseRecord):
    """CAVEAT field

    CAVEAT warns of severe errors in an entry. Use caution when using an entry
    containing this record.

    +---------+--------------+--------------+---------------------------------+
    | COLUMNS | DATA TYPE    | FIELD        | DEFINITION                      |
    +=========+==============+==============+=================================+
    | 1-6     | Record name  | "CAVEAT"     |                                 |
    +---------+--------------+--------------+---------------------------------+
    | 9-10    | Continuation | continuation | Allows concatenation of         |
    |         |              |              | multiple records.               |
    +---------+--------------+--------------+---------------------------------+
    | 12-15   | IDcode       | id_code      | PDB ID code of this entry.      |
    +---------+--------------+--------------+---------------------------------+
    | 20-79   | String       | comment      | Free text giving the reason for |
    |         |              |              | the CAVEAT.                     |
    +---------+--------------+--------------+---------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.id_code = None
        self.comment = []

    def parse_line(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_line(line)
        self.id_code = line[11:15].strip()
        self.comment.append(line[19:70].strip())

    def __str__(self):
        strings = []
        for iline, line in enumerate(self.comment):
            continuation = iline + 1
            if continuation > 1:
                strings += [f"CAVEAT  {continuation:>2}     {line:51}"]
            else:
                strings += ["CAVEAT        {line:51}"]
        return "\n".join(strings)


class Compound(BaseRecord):
    """COMPND field

    The COMPND record describes the macromolecular contents of an entry.
    Each macromolecule found in the entry is described by a set of token:
    value pairs, and is referred to as a COMPND record component. Since the
    concept of a molecule is difficult to specify exactly, PDB staff may
    exercise editorial judgment in consultation with depositors in
    assigning these names.

    For each macromolecular component, the molecule name, synonyms, number
    assigned by the Enzyme Commission (EC), and other relevant details are
    specified.

    +---------+---------------+--------------+--------------------------------+
    | COLUMNS | DATA TYPE     | FIELD        | DEFINITION                     |
    +=========+===============+==============+================================+
    | 1-6     | Record name   | "COMPND"     |                                |
    +---------+---------------+--------------+--------------------------------+
    | 8-10    | Continuation  | continuation | Allows concatenation of        |
    |         |               |              | multiple records.              |
    +---------+---------------+--------------+--------------------------------+
    | 11-80   | Specification | compound     | Description of the molecular   |
    |         | list          |              | components.                    |
    +---------+---------------+--------------+--------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.compound = []

    def parse_line(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_line(line)
        self.compound.append(line[10:80].strip())

    def __str__(self):
        strings = []
        for iline, line in enumerate(self.compound):
            continuation = iline + 1
            if continuation > 1:
                strings += [f"COMPND {continuation:>3} {line:69}"]
            else:
                strings += [f"COMPND    {line:70}"]
        return "\n".join(strings)


class ExperimentalData(BaseRecord):
    """EXPDTA field

    The EXPDTA record identifies the experimental technique used. This may
    refer to the type of radiation and sample, or include the spectroscopic
    or modeling technique. Permitted values include:

    * ELECTRON DIFFRACTION
    * FIBER DIFFRACTION
    * FLUORESCENCE TRANSFER
    * NEUTRON DIFFRACTION
    * NMR
    * THEORETICAL MODEL
    * X-RAY DIFFRACTION

    +---------+--------------+--------------+---------------------------------+
    | COLUMNS | DATA TYPE    | FIELD        | DEFINITION                      |
    +=========+==============+==============+=================================+
    | 1-6     | Record name  | "EXPDTA"     |                                 |
    +---------+--------------+--------------+---------------------------------+
    | 9-10    | Continuation | continuation | Allows concatenation of         |
    |         |              |              | multiple records.               |
    +---------+--------------+--------------+---------------------------------+
    | 11-79   | SList        | technique    | The experimental technique(s)   |
    |         |              |              | with optional comment           |
    |         |              |              | describing the sample or        |
    |         |              |              | experiment.                     |
    +---------+--------------+--------------+---------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.technique = []

    def parse_line(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_line(line)
        self.technique.append(line[10:79].strip())

    def __str__(self):
        strings = []
        for iline, line in enumerate(self.technique):
            continuation = iline + 1
            if continuation > 1:
                strings += [f"EXPDTA {continuation:>2} {line:79}"]
            else:
                strings += [f"EXPDTA    {line:79}"]
        return "\n".join(strings)


class Header(BaseRecord):
    """HEADER field

    The HEADER record uniquely identifies a PDB entry through the id_code
    field. This record also provides a classification for the entry. Finally,
    it contains the date the coordinates were deposited at the PDB.

    +---------+-------------+----------------+--------------------------------+
    | COLUMNS | DATA TYPE   | FIELD          | DEFINITION                     |
    +=========+=============+================+================================+
    | 1-6     | Record name | "HEADER"       |                                |
    +---------+-------------+----------------+--------------------------------+
    | 11-50   | String(40)  | classification | Classifies the molecule(s).    |
    +---------+-------------+----------------+--------------------------------+
    | 51-59   | Date        | dep_date       | Deposition date. This is the   |
    |         |             |                | date the coordinates  were     |
    |         |             |                | received at the PDB.           |
    +---------+-------------+----------------+--------------------------------+
    | 63-66   | IDcode      | id_code        | This identifier is unique      |
    |         |             |                | within the PDB.                |
    +---------+-------------+----------------+--------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.classification = None
        self.dep_date = None
        self.id_code = None

    def parse_line(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_line(line)
        self.classification = line[10:50].strip()
        self.dep_date = date_parse(line[50:59].strip())
        self.id_code = line[62:66].strip()

    def __str__(self):
        return (
            f"HEADER    "
            f"{self.classification:40}{date_format(self.dep_date):9}   "
            f"{self.id_code:4}"
        )


class Journal(BaseRecord):
    """JRNL field

    The JRNL record contains the primary literature citation that describes
    the experiment which resulted in the deposited coordinate set. There is
    at most one JRNL reference per entry. If there is no primary reference,
    then there is no JRNL reference. Other references are given in REMARK 1.

    +---------+-------------+--------+----------------------------------------+
    | COLUMNS | DATA TYPE   | FIELD  | DEFINITION                             |
    +=========+=============+========+========================================+
    | 1-6     | Record name | "JRNL" |                                        |
    +---------+-------------+--------+----------------------------------------+
    | 13-79   | LString     | text   | See details in PDB specification.      |
    +---------+-------------+--------+----------------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.text = None

    def parse_line(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_line(line)
        self.text = line[12:79].strip()

    def __str__(self):
        return f"JRNL        {self.text:67}"


class Keywords(BaseRecord):
    """KEYWDS field

    The KEYWDS record contains a set of terms relevant to the entry. Terms
    in the KEYWDS record provide a simple means of categorizing entries and
    may be used to generate index files. This record addresses some of the
    limitations found in the classification field of the HEADER record. It
    provides the opportunity to add further annotation to the entry in a
    concise and computer-searchable fashion.

    +---------+--------------+--------------+---------------------------------+
    | COLUMNS | DATA TYPE    | FIELD        | DEFINITION                      |
    +=========+==============+==============+=================================+
    | 1-6     | Record name  | "KEYWDS"     |                                 |
    +---------+--------------+--------------+---------------------------------+
    | 9-10    | Continuation | continuation | Allows concatenation of records |
    |         |              |              | if necessary.                   |
    +---------+--------------+--------------+---------------------------------+
    | 11-79   | List         | keywords     | Comma-separated list of         |
    |         |              |              | keywords relevant to the entry. |
    +---------+--------------+--------------+---------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.keywords = []

    def parse_line(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_line(line)
        self.keywords.append(line[10:80].strip())

    def __str__(self):
        strings = []
        for iline, line in enumerate(self.keywords):
            continuation = iline + 1
            if continuation > 1:
                strings += [f"KEYWDS  {continuation:>2} {line:79}"]
            else:
                strings += [f"KEYWDS    {line:79}"]
        return "\n".join(strings)


class ModelType(BaseRecord):
    """MDLTYP field.

    The MDLTYP record contains additional annotation pertinent to the
    coordinates presented in the entry.

    +---------+---------------+--------------+--------------------------------+
    | COLUMNS | DATA TYPE     | FIELD        | DEFINITION                     |
    +=========+===============+==============+================================+
    | 1-6     | Record name   | "MDLTYP"     |                                |
    +---------+---------------+--------------+--------------------------------+
    | 9-10    | Continuation  | continuation | Allows concatenation of        |
    |         |               |              | multiple records.              |
    +---------+---------------+--------------+--------------------------------+
    | 11-80   | SList         | comment      | Free Text providing additional |
    |         |               |              | structural annotation.         |
    +---------+---------------+--------------+--------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.comment = []

    def parse_line(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_line(line)
        self.comment.append(line[10:80].strip())

    def __str__(self):
        strings = []
        for iline, line in enumerate(self.comment):
            continuation = iline + 1
            if continuation > 1:
                strings += [f"MDLTYP  {continuation:>2} {line:79}"]
            else:
                strings += [f"MDLTYP    {line:80}"]
        return "\n".join(strings)


class Obsolete(BaseRecord):
    """OBSLTE field

    This record acts as a flag in an entry which has been withdrawn from the
    PDB's full release. It indicates which, if any, new entries have replaced
    the withdrawn entry.

    The format allows for the case of multiple new entries replacing one
    existing entry.

    +---------+--------------+---------------------+--------------------------+
    | COLUMNS | DATA TYPE    | FIELD               | DEFINITION               |
    +=========+==============+=====================+==========================+
    | 1-6     | Record name  | "OBSLTE"            |                          |
    +---------+--------------+---------------------+--------------------------+
    | 9-10    | Continuation | continuation        | Allows concatenation of  |
    |         |              |                     | multiple records         |
    +---------+--------------+---------------------+--------------------------+
    | 12-20   | Date         | replace_date        | Date that this entry was |
    |         |              |                     | replaced.                |
    +---------+--------------+---------------------+--------------------------+
    | 22-25   | IDcode       | id_code             | ID code of this entry.   |
    +---------+--------------+---------------------+--------------------------+
    | 32-35   | IDcode       | replace_id_codes[0] | ID of entry replacing    |
    |         |              |                     | this one.                |
    +---------+--------------+---------------------+--------------------------+
    | 37-40   | IDcode       | replace_id_codes[1] | ID of entry replacing    |
    |         |              |                     | this one.                |
    +---------+--------------+---------------------+--------------------------+
    | 42-45   | IDcode       | replace_id_codes[2] | ID of entry replacing    |
    |         |              |                     | this one.                |
    +---------+--------------+---------------------+--------------------------+
    | 47-50   | IDcode       | replace_id_codes[3] | ID of entry replacing    |
    |         |              |                     | this one.                |
    +---------+--------------+---------------------+--------------------------+
    | 52-55   | IDcode       | replace_id_codes[4] | ID of entry replacing    |
    |         |              |                     | this one.                |
    +---------+--------------+---------------------+--------------------------+
    | 57-60   | IDcode       | replace_id_codes[5] | ID of entry replacing    |
    |         |              |                     | this one.                |
    +---------+--------------+---------------------+--------------------------+
    | 62-65   | IDcode       | replace_id_codes[6] | ID of entry replacing    |
    |         |              |                     | this one.                |
    +---------+--------------+---------------------+--------------------------+
    | 67-70   | IDcode       | replace_id_codes[7] | ID of entry replacing    |
    |         |              |                     | this one.                |
    +---------+--------------+---------------------+--------------------------+
    | 72-75   | IDcode       | replace_id_codes[8] | ID of entry replacing    |
    |         |              |                     | this one.                |
    +---------+--------------+---------------------+--------------------------+
    """

    def __init__(self):
        super().__init__()
        self.replace_date = None
        self.id_code = None
        self.replace_id_codes = []

    def parse_line(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_line(line)
        self.replace_date = date_parse(line[11:20].strip())
        self.id_code = line[21:25].strip()
        self.replace_id_codes = [line[31:35].strip()]
        start = 36
        end = 40
        while True:
            id_code = line[start:end].strip()
            if id_code:
                self.replace_id_codes.append(id_code)
            start += 5
            end += 5
            if start > 67:
                break

    def __str__(self):
        strings = []
        err = f"This PDB is obsolete. Use one of the following instead:"
        err += f"{self.replace_id_codes}"
        _LOGGER.error(err)
        for ichunk, chunk in enumerate(grouper(self.replace_id_codes, 8)):
            continuation = ichunk + 1
            if continuation > 1:
                string = f"OBSLTE  {continuation:>2}"
            else:
                string = "OBSLTE    "
            string += (
                f" {date_format(self.replace_date):9} {self.id_code}     "
            )
            for code in chunk:
                if code is not None:
                    string += f" {code:4}"
            strings.append(string)
        return "\n".join(strings)


class Remark(BaseRecord):
    """REMARK field

    .. todo:: REMARK fields are horrible to parse. Someday we should implement.

    REMARK records present experimental details, annotations, comments, and
    information not included in other records. In a number of cases,
    REMARKs are used to expand the contents of other record types. A new
    level of structure is being used for some REMARK records. This is
    expected to facilitate searching and will assist in the conversion to a
    relational database.
    """

    def __init__(self):
        super().__init__()
        self.remark_num = None
        self.remark_text = None

    def parse_line(self, line):
        """Initialize by parsing line.

        +---------+------+-------------+--------------------------------------+
        | COLUMNS | TYPE | FIELD       | DEFINITION                           |
        +=========+======+=============+======================================+
        | 8-10    | int  | remark_num  | Remark number. It is not an error    |
        |         |      |             | for remark n to exist in an entry    |
        |         |      |             | when remark n-1 does not.            |
        +---------+------+-------------+--------------------------------------+
        | 12-79   | str  | remark_text | Left as white space in first line of |
        |         |      |             | each new remark.                     |
        +---------+------+-------------+--------------------------------------+

        :param str line:  line with PDB class
        """
        super().parse_line(line)
        self.remark_num = int(line[7:10].strip())
        self.remark_text = line[11:79]

    def __str__(self):
        return f"REMARK {self.remark_num:3} {self.remark_text:68}"


class Revision(BaseRecord):
    """Class to store contents of a single REVDAT modification.

    +---------+--------------+-------------------+---------------------------------+
    | COLUMNS | DATA TYPE    | FIELD             | DEFINITION                      |
    +=========+==============+===================+=================================+
    | 1-6     | Record name  | "REVDAT"          |                                 |
    +---------+--------------+-------------------+---------------------------------+
    | 8-10    | Integer      | modification_num  | Modification number.            |
    +---------+--------------+-------------------+---------------------------------+
    | 11-12   | Continuation | continuation      | Allows concatenation of         |
    |         |              |                   | multiple records.               |
    +---------+--------------+-------------------+---------------------------------+
    | 14-22   | Date         | modification_date | Date of modification (or        |
    |         |              |                   | for new entries) in DD-MMM-YY   |
    |         |              |                   | format. This is not repeated on |
    |         |              |                   | continued lines.                |
    +---------+--------------+-------------------+---------------------------------+
    | 24-27   | IDCode       | modification_id   | ID code of this entry. This is  |
    |         |              |                   | not repeated on continuation    |
    |         |              |                   | lines.                          |
    +---------+--------------+-------------------+---------------------------------+
    | 32      | Integer      | modification_type | An integer identifying the type |
    |         |              |                   | of modification. For all        |
    |         |              |                   | revisions, the modification     |
    |         |              |                   | type is listed as 1             |
    +---------+--------------+-------------------+---------------------------------+
    | 40-45   | LString(6)   | record            | Modification detail.            |
    +---------+--------------+-------------------+---------------------------------+
    | 47-52   | LString(6)   | record            | Modification detail.            |
    +---------+--------------+-------------------+---------------------------------+
    | 54-59   | LString(6)   | record            | Modification detail.            |
    +---------+--------------+-------------------+---------------------------------+
    | 61-66   | LString(6)   | record            | Modification detail.            |
    +---------+--------------+-------------------+---------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.modification_num = ""
        self.modification_date = None
        self.modification_id = ""
        self.modification_type = ""
        self.records = []

    def parse_line(self, line):
        """Parse PDB-format line for specific revision.

        :param str line:  line to parse.
        """
        super().parse_line(line)
        self.modification_num = int(line[7:10].strip())
        try:
            self.modification_date = date_parse(line[13:22].strip())
        except ValueError:
            pass
        mod_id = line[23:28].strip()
        if mod_id:
            self.modification_id = mod_id
        mod_type = line[31].strip()
        if mod_type:
            self.modification_type = int(mod_type)
        for start, end in [(39, 45), (46, 52), (53, 59), (60, 66)]:
            record = line[start:end].strip()
            if record:
                self.records.append(record)

    def __str__(self):
        if len(self.records) == 0:
            return (
                f"REVDAT {self.modification_num:>3}"
                f"   {date_format(self.modification_date):9} "
                f"{self.modification_id:4}    {self.modification_type:1}"
                f"      "
            )
        strings = []
        for ichunk, chunk in enumerate(grouper(self.records, 4)):
            continuation = ichunk + 1
            string = f"REVDAT {self.modification_num:>3}"
            if continuation > 1:
                string += f"{continuation:>2}                          "
            else:
                string += (
                    f"   {date_format(self.modification_date):9} "
                    f"{self.modification_id:4}    {self.modification_type:1}"
                    f"      "
                )
            for record in chunk:
                if record is not None:
                    string += f" {record:6}"
            strings.append(string.strip())
        return "\n".join(strings)


class RevisionData(BaseRecord):
    """REVDAT field

    REVDAT records contain a history of the modifications made to an entry
    since its release.

    +---------+--------------+--------------+---------------------------------+
    | COLUMNS | DATA TYPE    | FIELD        | DEFINITION                      |
    +=========+==============+==============+=================================+
    | 1-6     | Record name  | "REVDAT"     |                                 |
    +---------+--------------+--------------+---------------------------------+
    | 8-10    | Integer      | modNum       | Modification number.            |
    +---------+--------------+--------------+---------------------------------+
    | 11-12   | Continuation | continuation | Allows concatenation of         |
    |         |              |              | multiple records.               |
    +---------+--------------+--------------+---------------------------------+
    | 14-22   | Date         | modDate      | Date of modification (or        |
    |         |              |              | for new entries) in DD-MMM-YY   |
    |         |              |              | format. This is not repeated on |
    |         |              |              | continued lines.                |
    +---------+--------------+--------------+---------------------------------+
    | 24-27   | IDCode       | modId        | ID code of this entry. This is  |
    |         |              |              | not repeated on continuation    |
    |         |              |              | lines.                          |
    +---------+--------------+--------------+---------------------------------+
    | 32      | Integer      | modType      | An integer identifying the type |
    |         |              |              | of modification. For all        |
    |         |              |              | revisions, the modification     |
    |         |              |              | type is listed as 1             |
    +---------+--------------+--------------+---------------------------------+
    | 40-45   | LString(6)   | record       | Modification detail.            |
    +---------+--------------+--------------+---------------------------------+
    | 47-52   | LString(6)   | record       | Modification detail.            |
    +---------+--------------+--------------+---------------------------------+
    | 54-59   | LString(6)   | record       | Modification detail.            |
    +---------+--------------+--------------+---------------------------------+
    | 61-66   | LString(6)   | record       | Modification detail.            |
    +---------+--------------+--------------+---------------------------------+
    """

    def __init__(self):
        super().__init__()
        self._revisions = OrderedDict()

    @property
    def revisions(self) -> OrderedDict:
        """Get revisions.

        :returns: dictionary with modifiction numbers as keys and
            :class:`Revision` objects as values
        """
        return self._revisions

    @revisions.setter
    def revisions(self, value):
        self._revisions = value

    def parse_line(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_line(line)
        mod_num = int(line[7:10].strip())
        revision = self._revisions.get(mod_num, Revision())
        revision.parse_line(line)
        self._revisions[mod_num] = revision

    def __str__(self):
        strings = []
        curr_mod = None
        continuation = 1
        for mod_num, revision in self._revisions.items():
            string = str(revision)
            if mod_num == curr_mod:
                continuation += 1
                string = string[:11] + f"{continuation:>2}" + string[12:]
            else:
                continuation = 1
                curr_mod = mod_num
            strings.append(string)
        return "\n".join(strings)


class Site(BaseRecord):
    """SITE class

    The SITE records supply the identification of groups comprising
    important sites in the macromolecule.

    +---------+--------------+-----------+------------------------------------+
    | COLUMNS | DATA TYPE    | FIELD     | DEFINITION                         |
    +=========+==============+===========+====================================+
    | 1-6     | Record name  | "SITE  "  |                                    |
    +---------+--------------+-----------+------------------------------------+
    | 8-10    | Integer      | seq_num   | Sequence number.                   |
    +---------+--------------+-----------+------------------------------------+
    | 12-14   | LString(3)   | site_id   | Site name.                         |
    +---------+--------------+-----------+------------------------------------+
    | 16-17   | Integer      | num_res   | Number of residues that compose    |
    |         |              |           | the site.                          |
    +---------+--------------+-----------+------------------------------------+
    | 19-21   | Residue name | res_name1 | Residue name for first residue     |
    |         |              |           | that creates the site.             |
    +---------+--------------+-----------+------------------------------------+
    | 23      | Character    | chain_id1 | Chain identifier for first residue |
    |         |              |           | of site.                           |
    +---------+--------------+-----------+------------------------------------+
    | 24-27   | Integer      | seq1      | Residue sequence number for first  |
    |         |              |           | residue of the site.               |
    +---------+--------------+-----------+------------------------------------+
    | 28      | AChar        | ins_code1 | Insertion code for first residue   |
    |         |              |           | of the site.                       |
    +---------+--------------+-----------+------------------------------------+
    | 30-32   | Residue name | res_name2 | Residue name for second residue    |
    |         |              |           | that creates the site.             |
    +---------+--------------+-----------+------------------------------------+
    | 34      | Character    | chain_id2 | Chain identifier for second        |
    |         |              |           | residue of the site.               |
    +---------+--------------+-----------+------------------------------------+
    | 35-38   | Integer      | seq2      | Residue sequence number for second |
    |         |              |           | residue of the site.               |
    +---------+--------------+-----------+------------------------------------+
    | 39      | AChar        | ins_code2 | Insertion code for second residue  |
    |         |              |           | of the site.                       |
    +---------+--------------+-----------+------------------------------------+
    | 41-43   | Residue name | res_name3 | Residue name for third residue     |
    |         |              |           | that creates the site.             |
    +---------+--------------+-----------+------------------------------------+
    | 45      | Character    | chain_id3 | Chain identifier for third residue |
    |         |              |           | of the site.                       |
    +---------+--------------+-----------+------------------------------------+
    | 46-49   | Integer      | seq3      | Residue sequence number for third  |
    |         |              |           | residue of the site.               |
    +---------+--------------+-----------+------------------------------------+
    | 50      | AChar        | ins_code3 | Insertion code for third residue   |
    |         |              |           | of the site.                       |
    +---------+--------------+-----------+------------------------------------+
    | 52-54   | Residue name | res_name4 | Residue name for fourth residue    |
    |         |              |           | that creates the site.             |
    +---------+--------------+-----------+------------------------------------+
    | 56      | Character    | chain_id4 | Chain identifier for fourth        |
    |         |              |           | residue of the site.               |
    +---------+--------------+-----------+------------------------------------+
    | 57-60   | Integer      | seq4      | Residue sequence number for fourth |
    |         |              |           | residue of the site.               |
    +---------+--------------+-----------+------------------------------------+
    | 61      | AChar        | ins_code4 | Insertion code for fourth residue  |
    |         |              |           | of the site.                       |
    +---------+--------------+-----------+------------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.sites = OrderedDict()
        self.seq_num = None
        self.site_id = None
        self.num_res = None
        self.res_name1 = None
        self.chain_id1 = None
        self.seq1 = None
        self.ins_code1 = ""
        self.res_name2 = ""
        self.chain_id2 = ""
        self.seq2 = ""
        self.ins_code2 = ""
        self.res_name3 = ""
        self.chain_id3 = ""
        self.seq3 = ""
        self.ins_code3 = ""
        self.res_name4 = ""
        self.chain_id4 = ""
        self.seq4 = ""
        self.ins_code4 = ""

    def parse_line(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_line(line)
        self.seq_num = int(line[7:10].strip())
        self.site_id = line[11:14].strip()
        self.num_res = int(line[15:17].strip())
        self.res_name1 = line[18:21].strip()
        self.chain_id1 = line[22].strip()
        self.seq1 = int(line[23:27].strip())
        try:
            self.ins_code1 = line[27].strip()
            self.res_name2 = line[29:32].strip()
            self.chain_id2 = line[33].strip()
            self.seq2 = int(line[34:38].strip())
        except (IndexError, ValueError):
            pass
        try:
            self.ins_code2 = line[38].strip()
            self.res_name3 = line[40:43].strip()
            self.chain_id3 = line[44].strip()
            self.seq3 = int(line[45:49].strip())
        except (IndexError, ValueError):
            pass
        try:
            self.ins_code3 = line[49].strip()
            self.res_name4 = line[51:54].strip()
            self.chain_id4 = line[55].strip()
            self.seq4 = int(line[56:60].strip())
            self.ins_code4 = line[60].strip()
        except (IndexError, ValueError):
            pass

    def __str__(self):
        return (
            f"SITE   {self.seq_num:3} {self.site_id:3} {self.num_res:2}"
            f" {self.res_name1:>3} {self.chain_id1:1}{self.seq1:4}"
            f"{self.ins_code1:1} {self.res_name2:>3} {self.chain_id2:1}"
            f"{self.seq2:4}{self.ins_code2:1} {self.res_name3:>3}"
            f" {self.chain_id3:1}{self.seq3:4}{self.ins_code3:1}"
            f" {self.res_name4:>3} {self.chain_id4:1}{self.seq4:4}"
            f"{self.ins_code4:1}"
        )


class NumModels(BaseRecord):
    """NUMMDL field

    The NUMMDL record indicates total number of models in a PDB entry.

    +---------+-------------+--------------+----------------------------------+
    | COLUMNS | DATA TYPE   | FIELD        | DEFINITION                       |
    +=========+=============+==============+==================================+
    | 1-6     | Record name | "NUMMDL"     |                                  |
    +---------+-------------+--------------+----------------------------------+
    | 11-14   | Integer     | model_number | Number of models.                |
    +---------+-------------+--------------+----------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.model_number = None

    def parse_line(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_line(line)
        self.model_number = int(line[10:14])

    def __str__(self):
        return f"NUMMDL    {self.model_number:<4}"


class Source(BaseRecord):
    """SOURCE field

    The SOURCE record specifies the biological and/or chemical source of
    each biological molecule in the entry. Sources are described by both
    the common name and the scientific name, e.g., genus and species.
    Strain and/or cell-line for immortalized cells are given when they help
    to uniquely identify the biological entity studied.

    +---------+---------------+--------------+--------------------------------+
    | COLUMNS | DATA TYPE     | FIELD        | DEFINITION                     |
    +=========+===============+==============+================================+
    | 1-6     | Record name   | "SOURCE"     |                                |
    +---------+---------------+--------------+--------------------------------+
    | 8-10    | Continuation  | continuation | Allows concatenation of        |
    |         |               |              | multiple records.              |
    +---------+---------------+--------------+--------------------------------+
    | 11-79   | Specification | source       | Identifies the source of the   |
    |         | List          |              | macromolecule in a token:      |
    |         |               |              | value format.                  |
    +---------+---------------+--------------+--------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.source = []

    def parse_line(self, line):
        """Parse a PDB-format line.

        :param str line:  line to parse
        """
        super().parse_line(line)
        self.continuation = line[7:10].strip()
        self.source.append(line[10:79].strip())

    def __str__(self):
        strings = []
        for iline, line in enumerate(self.source):
            continuation = iline + 1
            if continuation > 1:
                strings += [f"SOURCE {continuation:>3} {line:79}"]
            else:
                strings += [f"SOURCE    {line:79}"]
        return "\n".join(strings)


class Split(BaseRecord):
    """SPLIT field

    The SPLIT record is used in instances where a specific entry composes
    part of a large macromolecular complex. It will identify the PDB entries
    that are required to reconstitute a complete complex.

    +---------+--------------+--------------+---------------------------------+
    | COLUMNS | DATA TYPE    | FIELD        | DEFINITION                      |
    +=========+==============+==============+=================================+
    | 1-6     | Record  name | "SPLIT "     |                                 |
    +---------+--------------+--------------+---------------------------------+
    | 9-10    | Continuation | continuation | Allows concatenation of         |
    |         |              |              | multiple records.               |
    +---------+--------------+--------------+---------------------------------+
    | 12-15   | IDcode       | id_codes[0]  | ID code of related entry.       |
    +---------+--------------+--------------+---------------------------------+
    | 17-20   | IDcode       | id_codes[1]  | ID code of related entry.       |
    +---------+--------------+--------------+---------------------------------+
    | 22-25   | IDcode       | id_codes[2]  | ID code of related entry.       |
    +---------+--------------+--------------+---------------------------------+
    | 27-30   | IDcode       | id_codes[3]  | ID code of related entry.       |
    +---------+--------------+--------------+---------------------------------+
    | 32-35   | IDcode       | id_codes[4]  | ID code of related entry.       |
    +---------+--------------+--------------+---------------------------------+
    | 37-40   | IDcode       | id_codes[5]  | ID code of related entry.       |
    +---------+--------------+--------------+---------------------------------+
    | 42-45   | IDcode       | id_codes[6]  | ID code of related entry.       |
    +---------+--------------+--------------+---------------------------------+
    | 47-50   | IDcode       | id_codes[7]  | ID code of related entry.       |
    +---------+--------------+--------------+---------------------------------+
    | 52-55   | IDcode       | id_codes[8]  | ID code of related entry.       |
    +---------+--------------+--------------+---------------------------------+
    | 57-60   | IDcode       | id_codes[9]  | ID code of related entry.       |
    +---------+--------------+--------------+---------------------------------+
    | 62-65   | IDcode       | id_codes[10] | ID code of related entry.       |
    +---------+--------------+--------------+---------------------------------+
    | 67-70   | IDcode       | id_codes[11] | ID code of related entry.       |
    +---------+--------------+--------------+---------------------------------+
    | 72-75   | IDcode       | id_codes[12] | ID code of related entry.       |
    +---------+--------------+--------------+---------------------------------+
    | 77-80   | IDcode       | id_codes[13] | ID code of related entry.       |
    +---------+--------------+--------------+---------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.id_codes = []

    def parse_line(self, line):
        """Parse input line.

        :param str line:  PDB-format line to parse
        """
        start = 11
        end = 15
        while True:
            code = line[start:end].strip()
            if code:
                self.id_codes.append(code)
            start += 5
            end += 5
            if start > 77:
                break

    def __str__(self):
        strings = []
        for ichunk, chunk in enumerate(grouper(self.id_codes, 14)):
            string = ""
            continuation = ichunk + 1
            if continuation > 1:
                string += f"\nSITE    {continuation:>2}"
            else:
                string += "SITE      "
            for code in chunk:
                if code is not None:
                    string += f" {code:4}"
            strings += [string]
        return "\n".join(strings)


class Supersedes(BaseRecord):
    """SPRSDE field

    The SPRSDE records contain a list of the ID codes of entries that were
    made obsolete by the given coordinate entry and withdrawn from the PDB
    release set. One entry may replace many. It is PDB policy that only the
    principal investigator of a structure has the authority to withdraw it.

    +---------+--------------+----------------+-------------------------------+
    | COLUMNS | DATA TYPE    | FIELD          | DEFINITION                    |
    +=========+==============+================+===============================+
    | 1-6     | Record name  | "SPRSDE"       |                               |
    +---------+--------------+----------------+-------------------------------+
    | 9-10    | Continuation | continuation   | Allows for multiple ID codes. |
    +---------+--------------+----------------+-------------------------------+
    | 12-20   | Date         | super_date     | Date entry superseded the     |
    |         |              |                | listed entries. This field is |
    |         |              |                | not copied on continuations.  |
    +---------+--------------+----------------+-------------------------------+
    | 22-25   | IDcode       | id_code        | ID code of this entry. This   |
    |         |              |                | field is not copied on        |
    |         |              |                | continuations.                |
    +---------+--------------+----------------+-------------------------------+
    | 32-35   | IDcode       | super_id_codes | ID code of superseded entry.  |
    +---------+--------------+----------------+-------------------------------+
    | 37-40   | IDcode       | super_id_codes | ID code of superseded entry.  |
    +---------+--------------+----------------+-------------------------------+
    | 42-45   | IDcode       | super_id_codes | ID code of superseded entry.  |
    +---------+--------------+----------------+-------------------------------+
    | 47-50   | IDcode       | super_id_codes | ID code of superseded entry.  |
    +---------+--------------+----------------+-------------------------------+
    | 52-55   | IDcode       | super_id_codes | ID code of superseded entry.  |
    +---------+--------------+----------------+-------------------------------+
    | 57-60   | IDcode       | super_id_codes | ID code of superseded entry.  |
    +---------+--------------+----------------+-------------------------------+
    | 62-65   | IDcode       | super_id_codes | ID code of superseded entry.  |
    +---------+--------------+----------------+-------------------------------+
    | 67-70   | IDcode       | super_id_codes | ID code of superseded entry.  |
    +---------+--------------+----------------+-------------------------------+
    | 72-75   | IDcode       | super_id_codes | ID code of superseded entry.  |
    +---------+--------------+----------------+-------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.super_date = None
        self.id_code = None
        self.super_id_codes = []

    def parse_line(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_line(line)
        self.super_date = date_parse(line[11:20].strip())
        self.id_code = line[21:25].strip()
        self.super_id_codes = [line[31:35].strip()]
        self.super_id_codes.append(line[36:40].strip())
        self.super_id_codes.append(line[41:45].strip())
        self.super_id_codes.append(line[46:50].strip())
        self.super_id_codes.append(line[51:55].strip())
        self.super_id_codes.append(line[56:60].strip())
        self.super_id_codes.append(line[61:65].strip())
        self.super_id_codes.append(line[66:70].strip())

    def __str__(self):
        strings = []
        for ichunk, chunk in enumerate(grouper(self.super_id_codes, 8)):
            continuation = ichunk + 1
            if continuation > 1:
                string = f"SPRSDE  {continuation:>2}                    "
            else:
                string = (
                    f"SPRSDE     {date_format(self.super_date):9} "
                    f"{self.id_code:4}     "
                )
            for code in chunk:
                string += f" {code:4}"
            strings.append(string.strip())
        return "\n".join(strings)


class Title(BaseRecord):
    """TITLE field

    The TITLE record contains a title for the experiment or analysis that
    is represented in the entry. It should identify an entry in the PDB in
    the same way that a title identifies a paper.

    +---------+--------------+--------------+---------------------------------+
    | COLUMNS | DATA TYPE    | FIELD        | DEFINITION                      |
    +=========+==============+==============+=================================+
    | 1-6     | Record name  | "TITLE "     |                                 |
    +---------+--------------+--------------+---------------------------------+
    | 9-10    | Continuation | continuation | Allows concatenation of         |
    |         |              |              | multiple records.               |
    +---------+--------------+--------------+---------------------------------+
    | 11-80   | String       | title        | Title of the  experiment.       |
    +---------+--------------+--------------+---------------------------------+
    """

    def __init__(self):
        super().__init__()
        self.title = []

    def parse_line(self, line):
        """Parse PDB-format line.

        :param str line:  line to parse
        """
        super().parse_line(line)
        self.title.append(line[10:80].strip())

    def __str__(self):
        strings = []
        for iline, line in enumerate(self.title):
            continuation = iline + 1
            if continuation > 1:
                strings += [f"TITLE   {continuation:>2} {line:69}"]
            else:
                strings += [f"TITLE     {line:70}"]
        return "\n".join(strings)
