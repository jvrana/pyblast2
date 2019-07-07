import pytest
from pyblast.blast.seqdb import SeqRecordDB
from lobio import SeqUtils
from Bio.SeqRecord import SeqRecord


def random_record():
    r = SeqRecord(SeqUtils.random_dna(10), annotations={"topology": "circular"})
    return r


def test_init():
    SeqRecordDB()


def test_add_record():
    db = SeqRecordDB()
    db.add(random_record())


def test_add_records():
    db = SeqRecordDB()
    db.add_many([random_record(), random_record()])


def test_get_record():
    db = SeqRecordDB()
    r = random_record()
    k = db.add(r)
    assert k
    assert db.get(k) is r


def test_get_records():
    db = SeqRecordDB()
    keys = db.add_many([random_record(), random_record()])
    assert keys
    assert db.get_many(keys)
