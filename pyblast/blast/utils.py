from glob import glob
from os.path import join
from Bio import SeqIO
import tempfile
from Bio.SeqFeature import FeatureLocation, CompoundLocation


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
