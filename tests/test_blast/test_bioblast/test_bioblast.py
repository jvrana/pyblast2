import json
from os.path import join

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pyblast import BioBlast
from pyblast.exceptions import PyBlastException
from pyblast.utils import force_unique_record_ids
from pyblast.utils import load_genbank_glob
from pyblast.utils import make_linear


def test_basic_run():
    junk1 = "atgctatgctgatgctgctgtgctgatgctgatgtgtattgctgtatcgcgcgagttagc"
    junk2 = "g" * 30
    frag = "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta"

    query = SeqRecord(seq=Seq(frag), annotations={"circular": False})
    subject = SeqRecord(seq=Seq(junk1 + frag + junk2), annotations={"circular": False})

    # print(type(query))
    # print(type(subject))
    blaster = BioBlast([subject], [query])
    blaster.blastn()
    alignments = blaster.results
    print(alignments)


def test_basic_run_reverse_complement():
    junk1 = "atgctatgctgatgctgctgtgctgatgctgatgtgtattgctgtatcgcgcgagttagc"
    junk2 = "g" * 30
    frag = "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta"

    query = SeqRecord(
        seq=Seq(frag), annotations={"circular": False}
    ).reverse_complement()
    subject = SeqRecord(seq=Seq(junk1 + frag + junk2), annotations={"circular": False})

    make_linear([query])
    # print(type(query))
    # print(type(subject))
    blaster = BioBlast([subject], [query])
    blaster.blastn()
    alignments = blaster.results
    for a in alignments:
        print(json.dumps(a, indent=2))
        assert a["subject"]["strand"] == -1


def test_run_bioblast_twice():
    junk1 = "atgctatgctgatgctgctgtgctgatgctgatgtgtattgctgtatcgcgcgagttagc"
    junk2 = "g" * 30
    frag = "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta"

    query = SeqRecord(seq=Seq(frag), annotations={"circular": False})
    subject = SeqRecord(seq=Seq(junk1 + frag + junk2), annotations={"circular": False})

    blaster = BioBlast([subject], [query])
    blaster.blastn()
    blaster.blastn()
    alignments = blaster.results
    print(alignments)


def test_perfect_match():
    s = "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta"
    r1 = SeqRecord(seq=Seq(s), annotations={"topology": "linear"})
    r2 = SeqRecord(seq=Seq(s), annotations={"topology": "linear"})

    blaster = BioBlast([r1], [r2])
    blaster.blastn()
    result = blaster.results[0]
    assert result["query"]["start"] == 1
    assert result["query"]["end"] == len(s)
    assert result["subject"]["start"] == 1
    assert result["subject"]["end"] == len(s)


# def test_valid_results(new_bio_blast):
#     blast = new_bio_blast()
#     blast.blastn()
#     Validator.validate_blaster_results(blast)


def test_short_blastn(new_primer_blast):
    blast = new_primer_blast()
    blast.blastn()

    assert blast.results


def test_blast_with_circular(new_circular_bio_blast):
    blast = new_circular_bio_blast()
    blast.blastn()
    assert blast.results


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

    results = bioblast.blastn()
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

    seqstr1 = str(subjects[0].seq)[:1000]
    seqstr2 = str(subjects[1].seq)[:1000]

    queries = [
        SeqRecord(Seq(seqstr1)),
        SeqRecord(Seq(seqstr2))
        # SeqRecord(Seq(str(subjects[1][:1000]))),
    ]
    make_linear(queries)
    with pytest.raises(PyBlastException):
        BioBlast(subjects, queries)


def test_unnamed_queries(here):
    subjects = load_genbank_glob(
        join(here, "data/test_data/genbank/templates/*.gb"), force_unique_ids=True
    )

    seqstr1 = str(subjects[0].seq)[:1000]
    seqstr2 = str(subjects[1].seq)[:1000]

    queries = [
        SeqRecord(Seq(seqstr1)),
        SeqRecord(Seq(seqstr2))
        # SeqRecord(Seq(str(subjects[1][:1000]))),
    ]
    force_unique_record_ids(make_linear(queries))

    bioblast = BioBlast(subjects, queries)
    results = bioblast.blastn()
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
    results = bioblast.blastn()
    assert not results


def test_with_flags(new_circular_bio_blast):
    blast = new_circular_bio_blast()
    blast.update_config({"ungapped": None})
    blast.blastn()
    assert blast.results


def test_ungapped():
    frag = "GtctaaaggtgaagaattattcactggtgttgtcccaattttggttgaattagatggtgatgttaatggtcacaaattttctgtctccggtgaaggtgaaggtgatgctacttacggtaaattgaccttaaaatttatttgtactactggtaaattgccagttccatggccaaccttagtcactactttcggttatggtgttcaatgttttgcgagatacccagatcatatgaaacaacatgactttttcaagtctgccatgccagaaggttatgttcaagaaagaactatttttttcaaagatgacggtaactacaagaccagagctgaagtcaagtttgaaggtgataccttagttaatagaatcgaattaaaaggtattgattttaaagaagatggtaacattttaggtcacaaattggaatacaactataactctcacaatgtttacatcatggctgacaaacaaaagaatggtatcaaagttaacttcaaaattagacacaacattgaagatggttctgttcaattagctgaccattatcaacaaaatactccaattggtgatggtccagtcttgttaccagacaaccattacttatccactcaatctgccttatccaaagatccaaacgaaaagagagaccacatggtcttgttagaatttgttactgctgctggtattacccatggtatggatgaattgtacaaaTAGTGATACCGTCGACCTCGAGTCAattagttatgtcacgcttacattcacgccctccccccacatccgctctaaccgaaaaggaaggagttagacaacctgaagtctaggtccctatttatttttttatagttatgttagtattaagaacgttatttatatttcaaatttttcttt"

    query = SeqRecord(seq=Seq(frag), annotations={"circular": False})
    subject = SeqRecord(
        seq=Seq(frag[:400] + "atgctatgctgatgctgctgtgctgat" + frag[400:]),
        annotations={"circular": False},
    )

    # print(type(query))
    # print(type(subject))
    blaster = BioBlast([subject], [query])
    blaster.update_config({"ungapped": None})
    blaster.blastn()
    alignments = blaster.results
    print(alignments)
