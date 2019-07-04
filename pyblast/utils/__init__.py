"""utils"""

from .cmd import run_cmd
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from glob import glob
from os.path import join
from Bio import SeqIO
import tempfile
from Bio.SeqFeature import FeatureLocation, CompoundLocation


def reverse_complement(seq_str):
    return str(Seq(seq_str, generic_dna).reverse_complement())


def complement(seq_str):
    return str(Seq(seq_str, generic_dna).complement())


def new_feature_location(start, end, length, strand):
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


def records_to_tmpfile(records):
    fd, tmp_path_handle = tempfile.mkstemp(suffix=".fasta")
    SeqIO.write(records, tmp_path_handle, format="fasta")
    return tmp_path_handle


def glob_fasta_to_tmpfile(dir):
    fasta_files = glob(join(dir, "*.fsa"))
    fasta_files += glob(join(dir, "*.fasta"))
    records = []
    for fsa in fasta_files:
        records += list(SeqIO.parse(fsa, "fasta"))
    return records_to_tmpfile(records)


def load_genbank_glob(path):
    records = []
    for f in glob(path):
        records += list(SeqIO.parse(f, format="genbank"))
    for r in records:
        if "topology" in r.annotations and r.annotations["topology"] == "circular":
            r.annotations["circular"] = True
        else:
            r.annotations["circular"] = False
    return records


def clean_records(records):
    for r in records:
        r.features = [f for f in r.features if f.location is not None]
