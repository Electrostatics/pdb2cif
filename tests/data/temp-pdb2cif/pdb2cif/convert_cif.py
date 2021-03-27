"""Convert mmCIF format to old PDB format."""
import logging
from datetime import datetime
from pprint import pprint, pformat
from old_pdb import pdb_entry, annotation
import old_pdb


_LOGGER = logging.getLogger(__name__)


def _skip(cif_obj, skipped) -> dict:
    """Store skipped CIF object

    :param :class:`pdbx.containers.DataCategory` cif_object:  CIF object to
        parse
    :param dict skipped:  list of skipped rows/attributes
    :returns:  modified ``skipped``
    """
    for iattr, attribute in enumerate(cif_obj.attribute_list):
        skipped[attribute] = [row[iattr] for row in cif_obj.row_list]
    return skipped


def _entry(entry, cif_obj, skipped) -> tuple:
    """Convert CIF object.

    :param :class:`old_pdb.pdb_entry.Entry` entry:  PDB entry object to update
    :param :class:`pdbx.containers.DataCategory` cif_object:  CIF object to
        parse
    :param dict skipped:  list of skipped rows/attributes
    :returns:  (modified ``entry``, modified ``skipped``)
    """
    for iattr, attribute in enumerate(cif_obj.attribute_list):
        values = [row[iattr] for row in cif_obj.row_list]
        if len(values) > 1:
            raise ValueError(f"Got {len(values)} instead of single value.")
        else:
            value = values[0]
        if attribute == "id":
            entry.header.id_code = value
        else:
            raise ValueError(f"Unable to parse attribute {attribute}.")
    return entry, skipped


def _pdbx_database_pdb_obs_spr(entry, cif_obj, skipped) -> tuple:
    """Convert CIF object.

    :param :class:`old_pdb.pdb_entry.Entry` entry:  PDB entry object to update
    :param :class:`pdbx.containers.DataCategory` cif_object:  CIF object to
        parse
    :param dict skipped:  list of skipped rows/attributes
    :returns:  (modified ``entry``, modified ``skipped``)
    """
    for iattr, attribute in enumerate(cif_obj.attribute_list):
        values = [row[iattr] for row in cif_obj.row_list]
        if len(values) > 1:
            raise ValueError(f"Got {len(values)} instead of single value.")
        else:
            value = values[0]
        if attribute in ["id", "details"]:
            if value is not None:
                skipped[attribute] = value
        elif attribute == "date":
            date = datetime.strptime(value, r"%Y-%m-%d")
            entry.obsolete.replace_date = date
            entry.supersedes.super_date = date
        elif attribute == "pdb_id":
            entry.obsolete.replace_id_codes.append(value)
            entry.supersedes.super_id_codes.append(value)
        elif attribute == "replace_pdb_id":
            entry.obsolete.id_code = value
            entry.supersedes.id_code = value
        else:
            raise ValueError(f"Unable to parse attribute {attribute}.")
    return entry, skipped


def _pdbx_database_related(entry, cif_obj, skipped) -> tuple:
    """Convert CIF object.

    :param :class:`old_pdb.pdb_entry.Entry` entry:  PDB entry object to update
    :param :class:`pdbx.containers.DataCategory` cif_object:  CIF object to
        parse
    :param dict skipped:  list of skipped rows/attributes
    :returns:  (modified ``entry``, modified ``skipped``)
    """
    for iattr, attribute in enumerate(cif_obj.attribute_list):
        values = [row[iattr] for row in cif_obj.row_list]
        if len(values) > 1:
            raise ValueError(f"Got {len(values)} instead of single value.")
        else:
            value = values[0]
        if attribute in ["db_name", "details", "content_type"]:
            if value is not None:
                skipped[attribute] = value
        elif attribute == "db_id":
            entry.split.id_codes.append(value)
        else:
            raise ValueError(f"Unable to parse attribute {attribute}.")
    return entry, skipped


def _audit_author(entry, cif_obj, skipped) -> tuple:
    """Convert CIF object.

    :param :class:`old_pdb.pdb_entry.Entry` entry:  PDB entry object to update
    :param :class:`pdbx.containers.DataCategory` cif_object:  CIF object to
        parse
    :param dict skipped:  list of skipped rows/attributes
    :returns:  (modified ``entry``, modified ``skipped``)
    """
    for iattr, attribute in enumerate(cif_obj.attribute_list):
        values = [row[iattr] for row in cif_obj.row_list]
        if attribute == "name":
            entry.author.author_list += values
        elif attribute == "pdbx_ordinal":
            skipped[attribute] = values
        else:
            raise ValueError(f"Unable to parse attribute {attribute}.")
    return entry, skipped


def _cell(entry, cif_obj, skipped) -> tuple:
    """Convert CIF object.

    :param :class:`old_pdb.pdb_entry.Entry` entry:  PDB entry object to update
    :param :class:`pdbx.containers.DataCategory` cif_object:  CIF object to
        parse
    :param dict skipped:  list of skipped rows/attributes
    :returns:  (modified ``entry``, modified ``skipped``)
    """
    for iattr, attribute in enumerate(cif_obj.attribute_list):
        values = [row[iattr] for row in cif_obj.row_list]
        if len(values) > 1:
            raise ValueError(f"Got {len(values)} instead of single value.")
        else:
            value = values[0]
        if attribute == "length_a":
            entry.unit_cell.a = value
        elif attribute == "length_b":
            entry.unit_cell.b = value
        elif attribute == "length_c":
            entry.unit_cell.c = value
        elif attribute == "angle_alpha":
            entry.unit_cell.alpha = value
        elif attribute == "angle_beta":
            entry.unit_cell.beta = value
        elif attribute == "angle_gamma":
            entry.unit_cell.gamma = value
        elif attribute == "Z_PDB":
            entry.unit_cell.z = value
        elif attribute in [
            "pdbx_unique_axis",
            "length_a_esd",
            "length_b_esd",
            "length_c_esd",
            "angle_alpha_esd",
            "angle_beta_esd",
            "angle_gamma_esd",
            "entry_id",
        ]:
            skipped[attribute] = value
        else:
            errstr = f"Unknown attribute {attribute} with value {value}."
            raise ValueError(errstr)
    return entry, skipped


def _symmetry(entry, cif_obj, skipped) -> tuple:
    """Convert CIF object.

    :param :class:`old_pdb.pdb_entry.Entry` entry:  PDB entry object to update
    :param :class:`pdbx.containers.DataCategory` cif_object:  CIF object to
        parse
    :param dict skipped:  list of skipped rows/attributes
    :returns:  (modified ``entry``, modified ``skipped``)
    """
    for iattr, attribute in enumerate(cif_obj.attribute_list):
        values = [row[iattr] for row in cif_obj.row_list]
        if len(values) > 1:
            raise ValueError(f"Got {len(values)} instead of single value.")
        else:
            value = values[0]
        if attribute == "space_group_name_H-M":
            entry.unit_cell.space_group = value
        elif attribute in [
            "entry_id",
            "pdbx_full_space_group_name_H-M",
            "cell_setting",
            "Int_Tables_number",
            "space_group_name_Hall",
        ]:
            skipped[attribute] = value
        else:
            errstr = f"Unknown attribute {attribute} with value {value}."
            raise ValueError(errstr)
    return entry, skipped


def _entity_poly_seq(entry, cif_obj, skipped) -> tuple:
    """Convert CIF object.

    :param :class:`old_pdb.pdb_entry.Entry` entry:  PDB entry object to update
    :param :class:`pdbx.containers.DataCategory` cif_object:  CIF object to
        parse
    :param dict skipped:  list of skipped rows/attributes
    :returns:  (modified ``entry``, modified ``skipped``)
    """
    if entry.sequence_residue is None:
        entry.sequence_residue = old_pdb.primary.SequenceResidues()
    residues = entry.sequence_residue.residues
    for row in cif_obj.row_list:
        entity_id = row[cif_obj.attribute_list.index("entity_id")]
        residues.get(entity_id, []).append(
            row[cif_obj.attribute_list.index("mon_id")]
        )
    entry.sequence_residue.residues = residues
    return entry, skipped


def _struct_ref(entry, cif_obj, skipped) -> tuple:
    """Convert CIF object.

    :param :class:`old_pdb.pdb_entry.Entry` entry:  PDB entry object to update
    :param :class:`pdbx.containers.DataCategory` cif_object:  CIF object to
        parse
    :param dict skipped:  list of skipped rows/attributes
    :returns:  (modified ``entry``, modified ``skipped``)
    """
    attr_list = cif_obj.attribute_list
    for row in cif_obj.row_list:
        print(attr_list)
        print(row)
        id_code = row[attr_list.index("id")]
        try:
            chain_id = row[attr_list.index("pdbx_strand_id")]
        except ValueError:
            chain_id = None
        seq_begin = row[attr_list.index("seq_align_beg")]
        ins_begin = row[attr_list.index("pdbx_seq_align_beg_ins_code")]
        seq_end = row[attr_list.index("seq_align_end")]
        ins_end = row[attr_list.index("pdbx_seq_align_end_ins_code")]
        database = row[attr_list.index("db_name")]
        database_accession = row[attr_list.index("pdbx_db_accession")]
        database_id_code = row[attr_list.index("db_code")]
        database_seq_begin = row[attr_list.index("db_align_beg")]
        database_seq_end = row[attr_list.index("db_align_end")]
        if len(database_accession) < 7:
            db_ref = old_pdb.primary.DatabaseReference()
            db_ref.id_code = id_code
            db_ref.chain_id = chain_id
            db_ref.seq_begin = seq_begin
            db_ref.ins_begin = ins_begin
            db_ref.seq_end = seq_end
            db_ref.ins_end = ins_end
            db_ref.database = database
            db_ref.database_accession = database_accession
            db_ref.database_id_code = database_id_code
            db_ref.database_seq_begin = database_seq_begin
            db_ref.database_seq_end = database_seq_end
            entry.database_reference.append(db_ref)
        else:
            db_ref1 = old_pdb.primary.DatabaseReference1()
            db_ref2 = old_pdb.primary.DatabaseReference2()
            db_ref1.id_code = id_code
            db_ref1.chain_id = chain_id
            db_ref1.seq_begin = seq_begin
            db_ref1.ins_begin = ins_begin
            db_ref1.seq_end = seq_end
            db_ref1.ins_end = ins_end
            db_ref1.database = database
            db_ref1.db_accession = database_accession
            db_ref1.db_id_code = database_id_code
            db_ref2.id_code = id_code
            db_ref2.chain_id = chain_id
            db_ref2.db_accession = database_accession
            db_ref2.seq_begin = database_seq_begin
            db_ref2.seq_end = database_seq_end
            entry.database_reference += [db_ref1, db_ref2]
    return entry, skipped


def _test(entry, cif_obj, skipped) -> tuple:
    """Convert CIF object.

    :param :class:`old_pdb.pdb_entry.Entry` entry:  PDB entry object to update
    :param :class:`pdbx.containers.DataCategory` cif_object:  CIF object to
        parse
    :param dict skipped:  list of skipped rows/attributes
    :returns:  (modified ``entry``, modified ``skipped``)
    """
    print(cif_obj.attribute_list)
    print(cif_obj.row_list)
    for iattr, attribute in enumerate(cif_obj.attribute_list):
        values = [row[iattr] for row in cif_obj.row_list]
        if len(values) > 1:
            raise ValueError(f"Got {len(values)} instead of single value.")
        else:
            value = values[0]
    raise NotImplementedError()
    return entry, skipped


def convert(containers) -> pdb_entry.Entry:
    """Convert a mmCIF entry into a old-format PDB entry.

    :param list containers:  list of :class:`pdbx.DataContainer` objects
    :returns:  old-format PDB entry
    """
    if len(containers) > 1:
        err = f"Can only handle 1 container; got {len(containers)}."
        raise NotImplementedError(err)
    entry = pdb_entry.Entry()
    container = containers[0]
    skipped_objects = {}
    entry, skipped["HEADER"] = header(entry, container)
    for name in object_names:
        cif_obj = container.get_object(name)
        if name == "entry":
            entry, skipped_objects[name] = _entry(
                entry, cif_obj, skipped_objects.get(name, {})
            )
        elif name == "pdbx_database_PDB_obs_spr":
            entry, skipped_objects[name] = _pdbx_database_pdb_obs_spr(
                entry, cif_obj, skipped_objects.get(name, {})
            )
        elif name == "pdbx_database_related":
            entry, skipped_objects[name] = _pdbx_database_related(
                entry, cif_obj, skipped_objects.get(name, {})
            )
        elif name == "audit_author":
            entry, skipped_objects[name] = _audit_author(
                entry, cif_obj, skipped_objects.get(name, {})
            )
        elif name == "cell":
            entry, skipped_objects[name] = _cell(
                entry, cif_obj, skipped_objects.get(name, {})
            )
        elif name == "symmetry":
            entry, skipped_objects[name] = _symmetry(
                entry, cif_obj, skipped_objects.get(name, {})
            )
        elif name == "entity_poly_seq":
            entry, skipped_objects[name] = _entity_poly_seq(
                entry, cif_obj, skipped_objects.get(name, {})
            )
        elif name == "struct_ref":
            entry, skipped_objects[name] = _struct_ref(
                entry, cif_obj, skipped_objects.get(name, {})
            )
        elif name in [
            "audit_conform",
            "database_2",
            "pdbx_database_status",
            "citation",
            "citation_author",
            "entity",
            "entity_name_com",
            "entity_poly",
            "entity_src_gen",
        ]:
            skipped_objects[name] = _skip(
                cif_obj, skipped_objects.get(name, {})
            )
        else:
            _LOGGER.error(f"Unexpected object:  {name}.")
            entry, skipped_objects[name] = _test(
                entry, cif_obj, skipped_objects.get(name, {})
            )
            errstr = f"Unable to parse object {name}."
            raise ValueError(errstr)

    raise NotImplementedError()
