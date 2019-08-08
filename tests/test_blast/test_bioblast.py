from os.path import join

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pyblast import BioBlast
from pyblast.exceptions import PyBlastException
from pyblast.utils import (
    make_linear,
    load_genbank_glob,
    force_unique_record_ids,
    make_circular,
)
import json


def test_basic_run():
    junk1 = "atgctatgctgatgctgctgtgctgatgctgatgtgtattgctgtatcgcgcgagttagc"
    junk2 = "g" * 30
    frag = "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta"

    query = SeqRecord(seq=Seq(frag), annotations={"circular": False})
    subject = SeqRecord(seq=Seq(junk1 + frag + junk2), annotations={"circular": False})

    # print(type(query))
    # print(type(subject))
    blaster = BioBlast([subject], [query])
    blaster.quick_blastn()
    alignments = blaster.alignments()
    print(alignments)


def test_perfect_match():
    s = "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta"
    r1 = SeqRecord(seq=Seq(s), annotations={"topology": "linear"})
    r2 = SeqRecord(seq=Seq(s), annotations={"topology": "linear"})

    blaster = BioBlast([r1], [r2])
    blaster.quick_blastn()
    result = blaster.results[0]
    assert result["query"]["start"] == 1
    assert result["query"]["end"] == len(s)
    assert result["subject"]["start"] == 1
    assert result["subject"]["end"] == len(s)


# def test_valid_results(new_bio_blast):
#     blast = new_bio_blast()
#     blast.quick_blastn()
#     Validator.validate_blaster_results(blast)


def test_short_blastn(new_primer_blast):
    blast = new_primer_blast()
    blast.quick_blastn()

    blast.results


def test_blast_with_circular(new_circular_bio_blast):
    blast = new_circular_bio_blast()
    blast.quick_blastn()
    results = blast.results


def test_raises_pyblast_when_not_unique(here):
    subjects = load_genbank_glob(join(here, "data/test_data/genbank/templates/*.gb"))
    queries = load_genbank_glob(join(here, "data/test_data/genbank/designs/*.gb"))
    print("n_queres: {}".format(len(queries)))
    print("n_subjects: {}".format(len(subjects)))
    with pytest.raises(PyBlastException):
        BioBlast(subjects, queries)


def test_not_raise_pyblast_when_unique(here):
    subjects = load_genbank_glob(join(here, "data/test_data/genbank/templates/*.gb"))
    queries = load_genbank_glob(join(here, "data/test_data/genbank/designs/*.gb"))

    force_unique_record_ids(subjects + queries)
    print("n_queres: {}".format(len(queries)))
    BioBlast(subjects, queries)


def test_multiquery_blast(here):
    subjects = load_genbank_glob(
        join(here, "data/test_data/genbank/templates/*.gb"), force_unique_ids=True
    )
    queries = load_genbank_glob(
        join(here, "data/test_data/genbank/designs/*.gb"), force_unique_ids=True
    )
    print("n_queres: {}".format(len(queries)))
    print("n_subjects: {}".format(len(subjects)))
    bioblast = BioBlast(subjects, queries)

    results = bioblast.quick_blastn()
    recids = set()
    for res in results:
        recid = res["query"]["origin_record_id"]
        recids.add(recid)
    print("n_records: {}".format(len(results)))
    assert len(recids) == len(queries)


def test_unnamed_queries_raises_duplicate_error(here):
    subjects = load_genbank_glob(
        join(here, "data/test_data/genbank/templates/*.gb"), force_unique_ids=True
    )
    queries = [
        SeqRecord(Seq(str(subjects[0][:1000].seq))),
        SeqRecord(Seq(str(subjects[1][:1000].seq))),
        # SeqRecord(Seq(str(subjects[1][:1000]))),
    ]
    make_linear(queries)
    with pytest.raises(PyBlastException):
        BioBlast(subjects, queries)


def test_unnamed_queries(here):
    subjects = load_genbank_glob(
        join(here, "data/test_data/genbank/templates/*.gb"), force_unique_ids=True
    )
    queries = [
        SeqRecord(Seq(str(subjects[0][:1000].seq))),
        SeqRecord(Seq(str(subjects[1][:1000].seq))),
        # SeqRecord(Seq(str(subjects[1][:1000]))),
    ]
    force_unique_record_ids(make_linear(queries))

    bioblast = BioBlast(subjects, queries)
    results = bioblast.quick_blastn()
    recids = set()
    for res in results:
        recid = res["query"]["origin_record_id"]
        recids.add(recid)
    print("n_records: {}".format(len(results)))
    assert len(recids) == len(queries)


def test_self_blast(here):
    subjects = load_genbank_glob(
        join(here, "data/test_data/genbank/templates/*.gb"), force_unique_ids=True
    )
    queries = [
        SeqRecord(Seq(str(subjects[0][:1000].seq))),
        # SeqRecord(Seq(str(subjects[1][:1000]))),
    ]
    force_unique_record_ids(make_linear(queries))

    bioblast = BioBlast(queries, queries)
    results = bioblast.quick_blastn()
    print(json.dumps(results[0], indent=2))


# TODO: move these types of test to a new file
# TODO: ensure proper start and end for queries and subjects for circular and linear plasmids and rc
def test_circular_over_query_origin():
    seq = "AAGGTCGGGCAGGAAGAGGGCCTATTTCCCATGATTCCTTCATATTTGCATATACGATACAAGGCTGTTAGAGAGATAATTAGAATTAATTTGACTGTAAACACAAAGATATTAGTACAAAATACGTGACGTAGAAAGTAATAATTTCTTGGGTAGTTTGCAGTTTTAAAATTATGTTTTAAAATGGACTATCATATGCTTACCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCTTGTGGAAAGGACG"

    queries = [SeqRecord(Seq(seq))]

    subjects = [SeqRecord(Seq(seq[-100:] + seq[:100]))]

    assert str(subjects[0].seq) in seq + seq

    bioblast = BioBlast(make_linear(subjects), make_circular(queries))
    assert len(bioblast.seq_db) == 4
    results = bioblast.quick_blastn()
    query_ends = []
    subject_ends = []

    for r in results:
        query_ends.append((r["query"]["start"], r["query"]["end"]))
        subject_ends.append((r["subject"]["start"], r["subject"]["end"]))
    print(bioblast.raw_results)

    subject_lengths = [x[1] - x[0] for x in subject_ends]
    assert len(subjects[0].seq) in subject_lengths
