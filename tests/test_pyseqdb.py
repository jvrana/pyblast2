from pyblast import PySeqDB
import os
from Bio import SeqIO


def test_init(here):
    directory = os.path.join(here, "data/test_data/templates")
    db = PySeqDB()
    db.add_from_directory(directory)
    assert len(db.db) > 0


def test_concat(here, tmpdir):
    directory = os.path.join(here, "data/test_data/templates")

    exp_num_seqs = len(os.listdir(directory))

    out = os.path.join(tmpdir.mkdir('concat'), 'concat.fasta')
    db = PySeqDB()
    db.add_from_directory(directory)
    db.concatenate_and_save(out)

    seqs = []
    with open(out, 'r') as handle:
        seqs += SeqIO.parse(handle, 'fasta')

    assert len(seqs) > 0
    assert len(seqs) == len(db.db)
    assert len(db.db) == exp_num_seqs


def test_ensure_db_returns_copy(here):
    directory = os.path.join(here, "data/test_data/templates")
    db = PySeqDB()
    db.add_from_directory(directory)

    assert len(db.db) > 0

    d = db.db
    d['foo'] = 'bar'
    assert d['foo'] == 'bar'
    assert 'foo' not in db.db