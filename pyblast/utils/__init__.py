"""utils"""

from .cmd import run_cmd
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq


def reverse_complement(seq_str):
    return str(Seq(seq_str, generic_dna).reverse_complement())


def complement(seq_str):
    return str(Seq(seq_str, generic_dna).complement())
