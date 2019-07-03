from glob import glob
from os.path import join
from Bio import SeqIO
import tempfile


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
