import os

import pytest

# from pyblast.seqio import open_sequence, sanitize_filenames, concat_seqs
from Bio.SeqRecord import SeqRecord
from pyblast.pysequence import PySequence


@pytest.fixture(scope="module")
def get_path(here):
    def wrapped(name):
        return os.path.join(here, 'data/test_data/seqio_plasmids', name)
    return wrapped


def test_circular(get_path):
    circular_path = get_path("circular_test_dna.gb")
    linear_path = get_path("linear_test_dna.gb")

    assert PySequence.dna_at_path_is_circular(circular_path)
    assert not PySequence.dna_at_path_is_circular(linear_path)


def test_parse_error():
    with pytest.raises(FileNotFoundError):
        PySequence.parse("alsdkjfaldsf")


def test_parse_linear(get_path):
    in_path = get_path("linear_test_dna.gb")
    seqs = PySequence.parse(in_path)
    assert len(seqs) == 1
    assert isinstance(seqs[0], PySequence)
    assert seqs[0].filename == in_path
    assert not seqs[0].circular


def test_parse_circular(get_path):
    in_path = get_path("circular_test_dna.gb")
    seqs = PySequence.parse(in_path)
    assert len(seqs) == 1
    assert isinstance(seqs[0], PySequence)
    assert seqs[0].filename == in_path
    assert seqs[0].circular




def test_sanitize_filename(here):
    directory = os.path.join(here, 'data/test_data/seqio_plasmids/to_be_sanitized')
    PySequence.sanitize_filenames(directory)
    for f in os.listdir(directory):
        assert ' ' not in f

    PySequence.sanitize_filenames(directory, replacements=[('_', ' ')])
    for f in os.listdir(directory):
        assert '_' not in f
