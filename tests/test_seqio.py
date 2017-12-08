import os

import pytest

from pyblast.seqio import open_sequence, sanitize_filenames, concat_seqs


@pytest.fixture(scope="module")
def get_path(here):
    def wrapped(name):
        return os.path.join(here, 'data/test_data/seqio_plasmids', name)
    return wrapped


def test_circular():
    assert get_path("circular_test_dna.gb")
    assert get_path("linear_test_dna.gb")


def test_open_sequence(get_path):
    with pytest.raises(ValueError):
        open_sequence(get_path("alsdkjfaldsf"))
    print(open_sequence(get_path("linear_test_dna.gb")))


def test_sanitize_filename(here):
    directory = os.path.join(here, 'data/test_data/seqio_plasmids/to_be_sanitized')
    sanitize_filenames(directory)
    for f in os.listdir(directory):
        assert not ' ' in f
        sanitize_filenames(directory, replacements=[('_', ' ')])
    for f in os.listdir(directory):
        assert not '_' in f


def test_concat_seqs(here):
    directory = os.path.join(here, "data/test_data/templates")
    out = os.path.join(here, "data/test_data/seqio_plasmids/concatenated/concat.fasta")
    concat_seqs(directory, out, savemeta=True)