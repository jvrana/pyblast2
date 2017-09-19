import pytest
import os
from pyblast import *

def test_circular():
    directory = os.path.abspath("tests/data/test_data/seqio_plasmids")
    assert dna_at_path_is_circular(os.path.join(directory, "circular_test_dna.gb")) == True
    assert dna_at_path_is_circular(os.path.join(directory, "linear_test_dna.gb")) == True

def test_open_sequence():
    directory = os.path.abspath("tests/data/test_data/seqio_plasmids")
    with pytest.raises(ValueError):
        open_sequence(os.path.join(directory, "alsdkjfaldsf"))
    print(open_sequence(os.path.join(directory, "linear_test_dna.gb")))

def test_sanitize_filename():
    directory = os.path.abspath("tests/data/test_data/seqio_plasmids/to_be_sanitized")
    sanitize_filenames(directory)
    for f in os.listdir(directory):
        assert not ' ' in f
        sanitize_filenames(directory, replacements=[('_', ' ')])
    for f in os.listdir(directory):
        assert not '_' in f

def test_concat_seqs():
    directory = os.path.abspath("tests/data/test_data/templates")
    out = os.path.abspath("tests/data/test_data/seqio_plasmids/concatenated/concat.fasta")
    concat_seqs(directory, out, savemeta=True)