import json
import random

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pyblast import BioBlast, BioBlastFactory
from pyblast.utils import make_linear, make_circular, force_unique_record_ids, Span


def random_sequence(length):
    bases = "ACGT"
    seq = ""
    for i in range(length):
        seq += bases[random.randint(0, 3)]
    return seq


def ns(nbases):
    return SeqRecord(Seq("N" * nbases))


def to_record(seq, linear=True):
    record = SeqRecord(Seq(seq))
    if linear:
        make_linear([record])
    else:
        make_circular([record])
    return record


def rand_record(nbases, linear=True):
    return to_record(random_sequence(nbases), linear=linear)


def compare_result(result, qs, qe, ss, se):
    x = [
        result["query"]["start"],
        result["query"]["end"],
        result["subject"]["start"],
        result["subject"]["end"],
    ]

    if ss is None:
        ss = result["subject"]["start"]
    if se is None:
        se = result["subject"]["end"]
    if qs is None:
        qs = result["subject"]["start"]
    if qe is None:
        qe = result["subject"]["end"]
    y = [qs, qe, ss, se]
    assert x == y
    return x == y


def test_perfect_alignment():
    record = rand_record(1000)
    queries = [record[:]]
    subjects = [record[:]]

    queries = make_linear(queries)
    subjects = make_linear(subjects)

    bioblast = BioBlast(subjects, queries)
    results = bioblast.quick_blastn()
    assert len(results) == 1
    compare_result(results[0], 1, len(record), 1, len(record))


def test_simple_alignment():
    record = rand_record(1000)
    queries = [record[:]]
    subjects = [record[10:-10]]

    queries = make_linear(queries)
    subjects = make_linear(subjects)

    bioblast = BioBlast(subjects, queries)
    results = bioblast.quick_blastn()
    assert len(results) == 1
    compare_result(results[0], 11, len(record) - 10, 1, len(record) - 10 - 10)


def test_align_Ns():
    record = rand_record(1000)
    nseq = SeqRecord(Seq("N" * 500))
    queries = [record[:]]
    subjects = [nseq + record + nseq]

    queries = make_linear(queries)
    subjects = make_linear(subjects)

    bioblast = BioBlast(subjects, queries)
    results = bioblast.quick_blastn()
    print(results)


@pytest.mark.parametrize("left_spacer", [0, 100, 200])
@pytest.mark.parametrize("ij", [(100, 300), (101, 300), (447, 888)])
def test_partial_alignment(left_spacer, ij):
    record = rand_record(1000)
    queries = [record[:]]
    subjects = [ns(left_spacer) + record[ij[0] : ij[1]]]

    queries = make_linear(queries)
    subjects = make_linear(subjects)

    bioblast = BioBlast(subjects, queries)
    results = bioblast.quick_blastn()
    assert len(results) == 1

    compare_result(
        results[0], ij[0] + 1, ij[1], 1 + left_spacer, ij[1] - ij[0] + left_spacer
    )


@pytest.mark.parametrize("left_spacer", [0, 100, 200])
@pytest.mark.parametrize("ij", [(100, 300), (101, 300), (447, 888)])
def test_partial_alignment_reverse_complement(left_spacer, ij):
    record = rand_record(1000)
    queries = [record[:]]
    subjects = [record[ij[0] : ij[1]]]
    subjects[0] = ns(left_spacer) + subjects[0].reverse_complement()

    queries = make_linear(queries)
    subjects = make_linear(subjects)

    bioblast = BioBlast(subjects, queries)
    results = bioblast.quick_blastn()
    assert len(results) == 1

    compare_result(results[0], ij[0] + 1, ij[1], len(subjects[0].seq), left_spacer + 1)


def test_reverse_alignment_simple():
    record = rand_record(1000)
    query = record
    subject = record[10:990].reverse_complement()

    subjects = make_linear([subject])
    queries = make_linear([query])

    bioblast = BioBlast(subjects, queries)
    results = bioblast.quick_blastn()
    for k, v in bioblast.seq_db.records.items():
        print(k)
        print(v)
    print(json.dumps(results, indent=2))
    assert results[0]["query"]["start"] == 10 + 1
    assert results[0]["query"]["end"] == 990
    assert results[0]["subject"]["start"] == 980
    assert results[0]["subject"]["end"] == 1


def test_reverse_alignment_simple2():
    record = rand_record(1000)
    query = record
    subject = record[10:-10].reverse_complement()

    subjects = make_linear([subject])
    queries = make_linear([query])

    factory = BioBlastFactory()
    factory.add_records(subjects, "subjects")
    factory.add_records(queries, "queries")

    bioblast = factory("subjects", "queries")

    results = bioblast.quick_blastn()
    for k, v in bioblast.seq_db.records.items():
        print(k)
        print(v)
    print(json.dumps(results, indent=2))
    assert results[0]["subject"]["strand"] == -1
    assert results[0]["subject"]["start"] == 980
    assert results[0]["subject"]["end"] == 1
    assert results[0]["query"]["start"] == 11
    assert results[0]["query"]["end"] == 990


class TestCircular:
    def test_circular_over_query(self):
        record = rand_record(1000)
        queries = [record]
        subjects = [record[-100:] + record[:100]]

        queries = make_circular(queries)
        subjects = make_linear(subjects)

        bioblast = BioBlast(subjects, queries)
        results = bioblast.quick_blastn()

        result = results[0]

        result_seq = str(
            (
                record[result["query"]["start"] - 1 :]
                + record[: result["query"]["end"]]
            ).seq
        )
        expected_seq = str(subjects[0].seq)
        assert result_seq == expected_seq

        compare_result(results[0], 1000 - 100 + 1, 100, 1, 200)

    def test_circular_over_query_revcomp(self):
        record = rand_record(1000)
        queries = [record]
        subjects = [record[-100:] + record[:100]]
        subjects[0] = subjects[0].reverse_complement()

        queries = make_circular(queries)
        subjects = make_linear(subjects)

        bioblast = BioBlast(subjects, queries)
        results = bioblast.quick_blastn()

        result = results[0]

        result_seq = str(
            (
                record[result["query"]["start"] - 1 :]
                + record[: result["query"]["end"]]
            ).seq
        )
        expected_seq = str(subjects[0].reverse_complement().seq)
        assert result_seq == expected_seq

        compare_result(results[0], 901, 100, 200, 1)

    def test_circular_over_subject(self):
        record = rand_record(1000)
        queries = [record]
        subjects = [record[200:300] + ns(500) + record[100:200]]

        queries = make_linear(queries)
        subjects = make_circular(subjects)

        bioblast = BioBlast(subjects, queries)
        results = bioblast.quick_blastn()

        compare_result(results[0], 101, 300, 601, 100)

    def test_circular_complete_subject(self):
        record = rand_record(1000)
        queries = [record]
        subjects = [record[500:] + record[:400]]

        queries = make_circular(queries)
        subjects = make_linear(subjects)

        bioblast = BioBlast(subjects, queries)
        results = bioblast.quick_blastn()

        result = results[0]
        print(result)
        assert result["subject"]["start"] == 1
        assert result["subject"]["end"] == 900
        assert result["query"]["start"] == 501
        assert result["query"]["end"] == 400

    def test_circular_complete_query_1(self):
        """In this situation, the subject is completely aligned with a circular query
        starting at index 500 (starting index = 0). Note that the
         pyblast results start at index 1."""
        record = rand_record(1000)
        queries = [record]
        subjects = [ns(100) + record[500:] + record[:500] + ns(100)]

        queries = make_circular(queries)
        subjects = make_linear(subjects)

        bioblast = BioBlast(subjects, queries)
        results = bioblast.quick_blastn()
        result = results[0]

        assert result["query"]["start"] == 501
        assert result["query"]["raw_end"] == 1500
        assert result["subject"]["start"] == 101
        assert result["subject"]["end"] == 1100

    def test_circular_complete_query_2(self):
        """In this situation, the subject is wraps around the query for 10 extra bases on the
        left site.
        Note that pyblast results start at index 1."""
        record = rand_record(1000)
        queries = [record]
        subjects = [ns(100) + record[-10 + 500 :] + record[:500] + ns(100)]

        queries = make_circular(queries)
        subjects = make_linear(subjects)

        bioblast = BioBlast(subjects, queries)
        results = bioblast.quick_blastn()
        result = results[0]

        assert result["query"]["start"] == 491
        assert result["query"]["raw_end"] == 1500
        assert result["subject"]["start"] == 101
        assert result["subject"]["end"] == 1110

    def test_circular_complete_query_3(self):
        """In this situation, the subject is wraps around the query for 10 extra bases on the
        right site.
        Note that pyblast results start at index 1."""
        record = rand_record(1000)
        queries = [record]
        subjects = [ns(100) + record[500:] + record[: 500 + 10] + ns(100)]

        queries = make_circular(queries)
        subjects = make_linear(subjects)

        bioblast = BioBlast(subjects, queries)
        results = bioblast.quick_blastn()
        result = results[0]

        assert result["query"]["start"] == 501
        assert result["query"]["raw_end"] == 1510
        assert result["subject"]["start"] == 101
        assert result["subject"]["end"] == 1110

    def test_circular_complete_query_4(self):
        """In this situation, the subject is wraps around the query for 10 extra bases on the
        left and right site.
        Note that pyblast results start at index 1."""
        record = rand_record(1000)
        queries = [record]
        subjects = [ns(100) + record[-10 + 500 :] + record[: 500 + 10] + ns(100)]

        queries = make_circular(queries)
        subjects = make_linear(subjects)

        bioblast = BioBlast(subjects, queries)
        results = bioblast.quick_blastn()
        result = results[0]

        assert result["query"]["start"] == 491
        assert result["query"]["raw_end"] == 1510
        assert result["subject"]["start"] == 101
        assert result["subject"]["end"] == 1120

    @pytest.mark.parametrize("extra_right", [0, 1, 10])
    @pytest.mark.parametrize("extra_left", [0, 1, 10])
    def test_circular_complete_query_parametrized(self, extra_right, extra_left):
        record = rand_record(1000)
        queries = [record]
        subjects = [
            ns(100)
            + record[(500 - extra_left) :]
            + record[: (500 + extra_right)]
            + ns(100)
        ]

        queries = make_circular(queries)
        subjects = make_linear(subjects)

        bioblast = BioBlast(subjects, queries)
        results = bioblast.quick_blastn()

        result = results[0]
        print(json.dumps(result, indent=2))

        assert result["query"]["start"] == 501 - extra_left
        assert result["query"]["raw_end"] == 1500 + extra_right
        assert result["subject"]["start"] == 101
        assert result["subject"]["end"] == 1100 + extra_right + extra_left

        # to spans
        query_span = bioblast.parse_result_to_span(result["query"], output_index=0)
        subject_span = bioblast.parse_result_to_span(result["subject"], output_index=0)

        assert len(subject_span) == len(query_span) == 1000 + extra_right + extra_left
        assert query_span.a == 500 - extra_left
        assert query_span.b == 500 + extra_right

        assert subject_span.a == 100
        assert subject_span.b == 1100 + extra_right + extra_left

    @pytest.mark.parametrize("extra_right", [0, 1, 10])
    @pytest.mark.parametrize("extra_left", [0, 1, 10])
    def test_circular_complete_query_parametrized_rc(self, extra_right, extra_left):
        record = rand_record(1000)
        queries = [record]
        subjects = [
            ns(100)
            + record[(500 - extra_left) :]
            + record[: (500 + extra_right)]
            + ns(100)
        ]

        subjects = [subjects[0].reverse_complement()]

        queries = make_circular(queries)
        subjects = make_linear(subjects)

        bioblast = BioBlast(subjects, queries)
        results = bioblast.quick_blastn()

        result = results[0]
        print(json.dumps(result, indent=2))

        assert result["query"]["start"] == 501 - extra_left
        assert result["query"]["raw_end"] == 1500 + extra_right
        assert result["subject"]["start"] == 1100 + extra_right + extra_left
        assert result["subject"]["end"] == 101

        # to spans
        query_span = bioblast.parse_result_to_span(result["query"], output_index=0)
        subject_span = bioblast.parse_result_to_span(result["subject"], output_index=0)

        assert len(subject_span) == len(query_span) == 1000 + extra_right + extra_left
        assert query_span.a == 500 - extra_left
        assert query_span.b == 500 + extra_right

        assert subject_span.a == 100
        assert subject_span.b == 1100 + extra_right + extra_left


# TODO: fix very long repeats
# def test_very_long_repeat(self):
#     record = rand_record(500)
#
#     queries = [record]
#     subjects = [record + record[:100]]
#
#     queries = make_circular(queries)
#     subjects = make_circular(subjects)
#
#     bioblast = BioBlast(subjects, queries)
#     results = bioblast.quick_blastn()
#
#     result = results[0]
#
#     compare_result(results[0], 1, 0, 1, 0)


def test_interaction_network():
    """We expect self alignments to be removed from the results."""
    records = [None, None, None, None]

    records[0] = rand_record(500)
    records[1] = rand_record(100) + records[0][:-100] + rand_record(1000)
    records[2] = rand_record(200) + records[1][:700] + rand_record(500)
    records[3] = records[2][-500:] + rand_record(500)

    force_unique_record_ids(records)

    queries = make_linear(records)

    bioblast = BioBlast(queries, queries)
    results = bioblast.quick_blastn()
    assert results
    for r in results:
        k1 = r["query"]["origin_key"]
        k2 = r["subject"]["origin_key"]
        print(k1, k2)
        assert not k1 == k2


def test_interaction_network_from_factory():
    """We expect self alignments to be removed from the results."""
    records = [None, None, None, None]

    records[0] = rand_record(500)
    records[1] = rand_record(100) + records[0][:-100] + rand_record(1000)
    records[2] = rand_record(200) + records[1][:700] + rand_record(500)
    records[3] = records[2][-500:] + rand_record(500)

    force_unique_record_ids(records)

    queries = make_linear(records)

    factory = BioBlastFactory()
    factory.add_records(queries, "queries")
    bioblast = factory("queries", "queries")

    results = bioblast.quick_blastn()
    assert results
    for r in results:
        k1 = r["query"]["origin_key"]
        k2 = r["subject"]["origin_key"]
        print(k1, k2)
        assert not k1 == k2


# def test_interaction_network_from_data(here):
#     from pyblast.utils import load_genbank_glob
#     from os.path import join
#
#     queries = load_genbank_glob(
#         join(here, "data/test_data/genbank/designs/*.gb"), force_unique_ids=True
#     )
#
#     templates = make_linear(queries[:1])
#     queries = make_linear([templates[0].reverse_complement()])
#
#     factory = BioBlastFactory()
#     factory.add_records(queries, "queries")
#     factory.add_records(templates, "templates")
#
#     blaster = factory("templates", "templates")
#
#     results = blaster.quick_blastn()
#     assert results
#     for r in results:
#         k1 = r["query"]["origin_key"]
#         k2 = r["subject"]["origin_key"]
#         print(k1, k2)
#         assert not k1 == k2
