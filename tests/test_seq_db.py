from pyblast.blast.seqdb import SeqRecordDB
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pyblast.utils import make_circular
from copy import deepcopy


def random_record():
    r = SeqRecord("ACGTCGTATGTATGTATTGATGATG", annotations={"topology": "circular"})
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


def test_add_with_transform():
    """We expect that when we add a record with a transformation, that the
    we can obtain the origin record and transformed record."""

    db = SeqRecordDB()
    from uuid import uuid4
    from pyblast.constants import Constants as C

    def pseudocircularize(r):
        r2 = r + r
        r2.name = C.PSEUDOCIRCULAR + "__" + r.name
        r2.id = str(uuid4())
        return r2

    record = SeqRecord(
        Seq("ACGTTCGTGATTGTGCTGTGTGTATGGTATGATTATAGTGATGTAGTGATGATGTAGTAGTATA")
    )
    records = make_circular([record])

    keys = db.add_many_with_transformations(
        records, pseudocircularize, C.PSEUDOCIRCULAR
    )

    key = keys[0]
    origin_key = db.get_origin_key(key)
    origin = db.get_origin(key)
    transformed = db.get(key)

    assert origin is record
    assert origin is not transformed
    assert origin_key is not key
    assert len(transformed) == 2 * len(record)


def test_add_same_transformation():
    """We expect that when we add a record with a transformation, that the
    we can obtain the origin record and transformed record."""

    db = SeqRecordDB()
    from uuid import uuid4
    from pyblast.constants import Constants as C

    def pseudocircularize(r):
        r2 = r + r
        r2.name = C.PSEUDOCIRCULAR + "__" + r.name
        r2.id = str(uuid4())
        return r2

    record = SeqRecord(
        Seq("ACGTTCGTGATTGTGCTGTGTGTATGGTATGATTATAGTGATGTAGTGATGATGTAGTAGTATA")
    )
    records = make_circular([record, record])

    keys = db.add_many_with_transformations(
        records, pseudocircularize, C.PSEUDOCIRCULAR
    )
    assert len(set(keys)) == 1
    assert len(db) == 2


def test_add_multiple_transformations():
    """We expect that when we add a record with a transformation, that the
    we can obtain the origin record and transformed record."""

    db = SeqRecordDB()
    from uuid import uuid4
    from pyblast.constants import Constants as C

    def pseudocircularize(r):
        pseudor = r + r
        pseudor.name = C.PSEUDOCIRCULAR + "__" + r.name
        pseudor.id = str(uuid4())
        return pseudor

    record = SeqRecord(
        Seq("ACGTTCGTGATTGTGCTGTGTGTATGGTATGATTATAGTGATGTAGTGATGATGTAGTAGTATA")
    )
    r1 = deepcopy(record)
    r2 = deepcopy(record)
    r1.id = "record1"
    r2.id = "record2"
    records = make_circular([r1, r2])

    keys = db.add_many_with_transformations(
        records, pseudocircularize, C.PSEUDOCIRCULAR
    )
    assert len(keys) == 2
    assert len(db) == 4
