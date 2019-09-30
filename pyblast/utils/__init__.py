from .cmd import run_cmd
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from glob import glob
from os.path import join
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import tempfile
from Bio.SeqFeature import FeatureLocation, CompoundLocation
from pyblast.constants import Constants as C
from typing import List
from uuid import uuid4
from .span import Span, SpanError


def reverse_complement(seq_str: str) -> str:
    """
    Return reverse complement of nucleotide sequence

    :param seq_str:
    :type seq_str:
    :return:
    :rtype:
    """
    return str(Seq(seq_str, generic_dna).reverse_complement())


def complement(seq_str: str) -> str:
    """
    Return the complement of the nucleotide sequence.

    :param seq_str:
    :type seq_str:
    :return:
    :rtype:
    """
    return str(Seq(seq_str, generic_dna).complement())


def new_feature_location(start: int, end: int, length: int, strand: int):
    """Makes a FeatureLocation. If necessary, makes CompoundLocation."""
    if start > end:
        if length is None:
            raise ValueError(
                "A length must be provided to create a feature with start > end."
            )
        f1 = FeatureLocation(start, length, strand)
        f2 = FeatureLocation(1, end, strand)
        if strand == -1:
            location = CompoundLocation([f2, f1])
        else:
            location = CompoundLocation([f1, f2])
    else:
        location = FeatureLocation(start, end, strand=strand)
    return location


def records_to_tmpfile(records: List[SeqRecord]) -> str:
    """
    Write SeqRecords to a temporary file.

    :param records: list of SeqRecords
    :type records: list
    :return: temporary file path
    :rtype: str
    """
    fd, tmp_path_handle = tempfile.mkstemp(suffix=".fasta")
    SeqIO.write(records, tmp_path_handle, format="fasta")
    return tmp_path_handle


def glob_fasta_to_tmpfile(dirpath: str) -> str:
    """
    Concatenate all fasta files into a temporary fasta file.

    :param dirpath: directory containing fasta files
    :type dirpath: str
    :return: the temporary filename
    :rtype: str
    """
    fasta_files = glob(join(dirpath, "*.fsa"))
    fasta_files += glob(join(dirpath, "*.fasta"))
    records = []
    for fsa in fasta_files:
        records += list(SeqIO.parse(fsa, "fasta"))
    return records_to_tmpfile(records)


def _grouped_by_id(records):
    """
    Return a dictionary of records grouped by their record_id
    :param records:
    :type records:
    :return:
    :rtype:
    """
    rec_id_dict = {}
    for r in records:
        rec_id_dict.setdefault(r.id, list()).append(r)


def force_unique_record_ids(records, use_uuid=False):
    """
    Force unique ids for all records in the list

    :param records: list of SeqRecords
    :type records: list
    :param use_uuid: If True, will append a unique identifier to the front of the record id. Else, enumerate the records.
    :type use_uuid: bool
    :return: the records
    :rtype: list
    """
    rec_id_dict = {}
    for r in records:
        rec_id_dict.setdefault(r.id, list()).append(r)
    for k, v in rec_id_dict.items():
        if len(v) > 1:
            for i, r in enumerate(v):
                if use_uuid:
                    i = str(uuid4())
                r.id = "({})__{}".format(i, r.id)
    return records


def load_glob(
    path: str, format_str: str, recursive=False, force_unique_ids=False
) -> List[SeqRecord]:
    """
    Load SeqRecords from a glob-like path.

    :param path: glob-like filepath
    :type path: str
    :param format_str: Fileformat. e.g. 'genbank', 'fasta'
    :type format_str: basestring
    :param recursive: whether to search recusively in the glob path.
    :type recursive: bool
    :param force_unique_ids: If True, force records to have unique ids
    :type force_unique_ids: bool
    :return: the list of SeqRecords
    :rtype: list
    """
    records = []
    filenames = glob(path, recursive=recursive)
    for f in filenames:
        new_records = list(SeqIO.parse(f, format=format_str))
        for rec in new_records:
            rec.from_file = f
            records.append(rec)

    if force_unique_ids:
        force_unique_record_ids(records)
    return records


def load_fasta_glob(path, recursive=False, force_unique_ids=False):
    """
    Load SeqRecords from a glob-like path containing fasta files.

    :param path: glob-like filepath
    :type path: str
    :param recursive: whether to search recusively in the glob path.
    :type recursive: bool
    :param force_unique_ids: If True, force records to have unique ids
    :type force_unique_ids: bool
    :return: the list of SeqRecords
    :rtype: list
    """
    return load_glob(
        path, "fasta", recursive=recursive, force_unique_ids=force_unique_ids
    )


def load_genbank_glob(path, recursive=False, force_unique_ids=False):
    """
    Load SeqRecords from a glob-like path containing genbank files.

    :param path: glob-like filepath
    :type path: str
    :param recursive: whether to search recusively in the glob path.
    :type recursive: bool
    :param force_unique_ids: If True, force records to have unique ids
    :type force_unique_ids: bool
    :return: the list of SeqRecords
    :rtype: list
    """
    records = load_glob(
        path, "genbank", recursive=recursive, force_unique_ids=force_unique_ids
    )
    for r in records:
        if C.TOPOLOGY in r.annotations and r.annotations[C.TOPOLOGY] == C.CIRCULAR:
            r.annotations[C.CIRCULAR] = True
        else:
            r.annotations[C.CIRCULAR] = False
    return records


def clean_records(records: List[SeqRecord]) -> None:
    """
    Remove features for records that have no location.
    """
    for r in records:
        r.features = [f for f in r.features if f.location is not None]


def is_circular(record: SeqRecord) -> bool:
    """
    Returns whether a "topology" annotation was found in the SeqRecord and
    whether this topolgy indicates that the sequence is circular.

    :param record: the SeqRecord
    :type record: SeqRecord
    :return: whether the SeqRecord is circular
    :rtype: bool
    """
    annotations = {k.lower(): v for k, v in record.annotations.items()}
    return (
        annotations.get(C.CIRCULAR.lower(), False) is True
        or annotations.get(C.LINEAR.lower(), True) is False
        or annotations.get(C.TOPOLOGY.lower(), C.LINEAR).lower() == C.CIRCULAR.lower()
    )


def make_linear(records: List[SeqRecord]) -> List[SeqRecord]:
    """
    Annotates the SeqRecords as linear by adding the 'topology' annotation and
    setting it to 'linear'

    :param records:
    :type records:
    :return: the records
    :rtype: list
    """
    for r in records:
        if C.CIRCULAR in r.annotations:
            r.annotations[C.CIRCULAR] = False
        elif C.LINEAR in r.annotations:
            r.annotations[C.LINEAR] = True
        else:
            r.annotations[C.TOPOLOGY] = C.LINEAR
    return records


def make_circular(records: List[SeqRecord]) -> List[SeqRecord]:
    """
    Annotates the SeqRecords as linear by adding the 'topology' annotation and
    setting it to 'linear'

    :param records:
    :type records:
    :return: the records
    :rtype: list
    """
    for r in records:
        if C.CIRCULAR in r.annotations:
            r.annotations[C.CIRCULAR] = True
        elif C.LINEAR in r.annotations:
            r.annotations[C.LINEAR] = False
        else:
            r.annotations[C.TOPOLOGY] = C.CIRCULAR
    return records
