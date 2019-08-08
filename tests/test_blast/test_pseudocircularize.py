import json

import pytest

from pyblast.blast import JSONBlast
from pyblast.utils import reverse_complement


class TestPseudocircularize:
    @pytest.fixture
    def seq_example(self):
        return "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta"

    @pytest.fixture
    def x1(self):
        return -30

    @pytest.fixture
    def x2(self):
        return 28

    @pytest.fixture
    def frag1(self, seq_example, x1):
        return seq_example[x1:]

    @pytest.fixture
    def frag2(self, seq_example, x2):
        return seq_example[:x2]

    @pytest.fixture
    def seqs(self, seq_example, frag1, frag2):

        query = {
            "id": "myquery",
            "bases": seq_example,
            "name": "myseq",
            "circular": True,
        }

        subject = {
            "id": "mysubject",
            "bases": frag1 + frag2,
            "name": "myseq",
            "circular": True,
        }

        return [subject], query

    @pytest.fixture
    def alignments_span_origin_false(self, seqs):
        subjects, query = seqs
        j = JSONBlast(subjects, query, span_origin=False)
        j.quick_blastn()
        results = j.get_perfect()
        return results

    @pytest.fixture
    def alignments_span_origin_true(self, seqs):
        subjects, query = seqs
        j = JSONBlast(subjects, query, span_origin=True)
        j.quick_blastn()
        results = j.get_perfect()
        return results

    def test_alignments(self, seqs):
        """Expect alignments to be returned"""
        subjects, query = seqs
        j = JSONBlast(subjects, query, span_origin=False)
        j.quick_blastn()
        assert j.results

    def test_pseudocircularized_alignments(self, seqs):
        """Expect alignments to be returned"""
        subjects, query = seqs
        j = JSONBlast(subjects, query, span_origin=True)
        j.quick_blastn()
        assert j.results

    @pytest.fixture
    def alignments_span_origin_true_query_not_circular(self, seqs):
        subjects, query = seqs
        query["circular"] = False
        j = JSONBlast(subjects, query, span_origin=False)
        j.quick_blastn()
        results = j.get_perfect()
        return results

    def test_sequence_size_doesnt_change(self, seqs):
        subjects, query = seqs
        j = JSONBlast(subjects, query, span_origin=False)
        j.quick_blastn()

        assert len(j.subject_json[0]["bases"]) == len(subjects[0]["bases"])
        assert len(j.query_json["bases"]) == len(query["bases"])
        assert j.results[0]["query"]["length"] == len(query["bases"])
        assert j.results[0]["subject"]["length"] == len(subjects[0]["bases"])

    def test_span_origin_false(self, alignments_span_origin_false):
        assert alignments_span_origin_false[0]["meta"]["span_origin"] is False

    def test_span_origin_true(self, alignments_span_origin_true):
        assert alignments_span_origin_true[0]["meta"]["span_origin"] is True

    def test_json_pseudocircularize_is_false(
        self, alignments_span_origin_false, frag1, frag2, seqs
    ):
        assert 2 == len(alignments_span_origin_false)
        assert (
            alignments_span_origin_false[0]["subject"]["bases"].upper() == frag1.upper()
        )
        assert (
            alignments_span_origin_false[1]["subject"]["bases"].upper() == frag2.upper()
        )

        subjects, query = seqs
        alignment_lengths = [
            a["meta"]["alignment length"] for a in alignments_span_origin_false
        ]
        assert len(subjects[0]["bases"]) not in alignment_lengths

    def test_json_pseudocircularize_is_true_but_sequences_are_not_circular(
        self, alignments_span_origin_true_query_not_circular, frag1, frag2, seqs
    ):
        """Expect same results as pseudocircularize==False if sequences are not circular"""
        assert len(alignments_span_origin_true_query_not_circular) == 2
        assert (
            alignments_span_origin_true_query_not_circular[0]["subject"][
                "bases"
            ].upper()
            == frag1.upper()
        )
        assert (
            alignments_span_origin_true_query_not_circular[1]["subject"][
                "bases"
            ].upper()
            == frag2.upper()
        )

        subjects, query = seqs
        alignment_lengths = [
            a["meta"]["alignment length"]
            for a in alignments_span_origin_true_query_not_circular
        ]
        assert len(subjects[0]["bases"]) not in alignment_lengths

    def test_json_pseudocircularize_is_true(self, alignments_span_origin_true, seqs):
        subjects, query = seqs
        # assert len(j.results) == 3

        # make sure size reflects the pseudocircularized sequences
        for align in alignments_span_origin_true:
            assert align["subject"]["length"] == len(subjects[0]["bases"]) * 2
            assert align["query"]["length"] == len(query["bases"]) * 2

        # assert full subject exists even if it spans over the origin
        alignment_lengths = [
            a["meta"]["alignment length"] for a in alignments_span_origin_true
        ]
        assert len(subjects[0]["bases"]) in alignment_lengths

    def test_assert_query_and_subject_bases_are_equal(
        self, alignments_span_origin_true
    ):
        for align in alignments_span_origin_true:
            assert align["query"]["bases"] == align["subject"]["bases"]


class TestPseudocircularize_SwitchQueryAndSubject:
    @pytest.fixture
    def seq_example(self):
        return "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta"

    @pytest.fixture
    def x1(self):
        return -30

    @pytest.fixture
    def x2(self):
        return 28

    @pytest.fixture
    def frag1(self, seq_example, x1):
        return seq_example[x1:]

    @pytest.fixture
    def frag2(self, seq_example, x2):
        return seq_example[:x2]

    @pytest.fixture
    def seqs(self, seq_example, frag1, frag2):

        subject = {
            "id": "myquery",
            "bases": seq_example,
            "name": "myseq",
            "circular": True,
        }

        query = {
            "id": "mysubject",
            "bases": frag1 + frag2,
            "name": "myseq",
            "circular": True,
        }

        return [subject], query

    def test_alignments(self, seqs):
        """Expect alignments to be returned"""
        subjects, query = seqs
        j = JSONBlast(subjects, query, span_origin=False)
        j.quick_blastn()
        assert j.results

    def test_pseudocircularized_alignments(self, seqs):
        """Expect alignments to be returned"""
        subjects, query = seqs
        j = JSONBlast(subjects, query, span_origin=True)
        j.quick_blastn()
        assert j.results

    @pytest.fixture
    def alignments_span_origin_false(self, seqs):
        subjects, query = seqs
        j = JSONBlast(subjects, query, span_origin=False)
        j.quick_blastn()
        results = j.get_perfect()
        return results

    @pytest.fixture
    def alignments_span_origin_true(self, seqs):
        subjects, query = seqs
        j = JSONBlast(subjects, query, span_origin=True)
        j.quick_blastn()
        results = j.get_perfect()
        return results

    @pytest.fixture
    def alignments_span_origin_true_query_not_circular(self, seqs):
        subjects, query = seqs
        query["circular"] = False
        j = JSONBlast(subjects, query, span_origin=False)
        j.quick_blastn()
        results = j.get_perfect()
        return results

    def test_sequence_size_doesnt_change(self, seqs):
        subjects, query = seqs
        j = JSONBlast(subjects, query, span_origin=False)
        j.quick_blastn()

        assert len(j.subject_json[0]["bases"]) == len(subjects[0]["bases"])
        assert len(j.query_json["bases"]) == len(query["bases"])
        assert j.results[0]["query"]["length"] == len(query["bases"])
        assert j.results[0]["subject"]["length"] == len(subjects[0]["bases"])

    def test_json_pseudocircularize_is_false(
        self, alignments_span_origin_false, frag1, frag2, seqs
    ):
        assert 2 == len(alignments_span_origin_false)
        assert (
            alignments_span_origin_false[0]["subject"]["bases"].upper() == frag1.upper()
        )
        assert (
            alignments_span_origin_false[1]["subject"]["bases"].upper() == frag2.upper()
        )

        subjects, query = seqs
        alignment_lengths = [
            a["meta"]["alignment length"] for a in alignments_span_origin_false
        ]
        assert len(subjects[0]["bases"]) not in alignment_lengths

    def test_json_pseudocircularize_is_true_but_sequences_are_not_circular(
        self, alignments_span_origin_true_query_not_circular, frag1, frag2, seqs
    ):
        """Expect same results as pseudocircularize==False if sequences are not circular"""
        assert len(alignments_span_origin_true_query_not_circular) == 2
        assert (
            alignments_span_origin_true_query_not_circular[0]["subject"][
                "bases"
            ].upper()
            == frag1.upper()
        )
        assert (
            alignments_span_origin_true_query_not_circular[1]["subject"][
                "bases"
            ].upper()
            == frag2.upper()
        )

        subjects, query = seqs
        alignment_lengths = [
            a["meta"]["alignment length"]
            for a in alignments_span_origin_true_query_not_circular
        ]
        assert len(subjects[0]["bases"]) not in alignment_lengths

    def test_json_pseudocircularize_is_true(
        self, alignments_span_origin_true, seqs, x1, x2
    ):
        subjects, query = seqs
        # assert len(j.results) == 3

        # make sure size reflects the pseudocircularized sequences
        for align in alignments_span_origin_true:
            assert align["subject"]["length"] == len(subjects[0]["bases"]) * 2
            assert align["query"]["length"] == len(query["bases"]) * 2

        # assert full subject exists even if it spans over the origin
        alignment_lengths = [
            a["meta"]["alignment length"] for a in alignments_span_origin_true
        ]
        assert len(query["bases"]) in alignment_lengths

    def test_assert_query_and_subject_bases_are_equal(
        self, alignments_span_origin_true
    ):
        for align in alignments_span_origin_true:
            assert align["query"]["bases"] == align["subject"]["bases"]


class TestPseudocirculariseWithLongSeqs:
    @pytest.fixture
    def frag1(self):
        return "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta"

    @pytest.fixture
    def frag2(self, frag1):
        return frag1[-40:]

    @pytest.fixture
    def long_seqs(self, frag1, frag2):
        """Alignments that wraps around the query more than onces"""

        junk1 = "atgctatgctgatgctgctgtgctgatgctgatgtgtattgctgtatcgcgcgagttagc"
        junk2 = "g" * 30

        query = {"id": "myquery", "bases": frag1, "name": "myseq", "circular": True}

        subject = {
            "id": "mysubject",
            "bases": junk1 + reverse_complement(frag2 + frag1) + junk2,
            "name": "myseq",
            "circular": True,
        }

        return [subject], query

    @pytest.fixture
    def long_seqs_alignment(self, long_seqs):
        subjects, query = long_seqs
        j = JSONBlast(
            subjects,
            query,
            span_origin=True,
            gapopen=3,
            gapextend=3,
            penalty=-5,
            reward=1,
        )
        j.quick_blastn()
        results = j.get_perfect()
        return results

    def test_json_pseudocircularize_is_true_long_seqs(self, long_seqs_alignment):
        results = long_seqs_alignment
        print(results)
        for i, align in enumerate(results):
            print(">>> Alignment {}".format(i))
            print(json.dumps(align, indent=4))

    def test_long_seqs_query_bases_length_equals_alignment_length(
        self, long_seqs_alignment
    ):
        for align in long_seqs_alignment:
            assert align["meta"]["alignment length"] == len(align["query"]["bases"])

    def test_long_seqs_subject_bases_length_equals_alignment_length(
        self, long_seqs_alignment
    ):
        for align in long_seqs_alignment:
            assert align["meta"]["alignment length"] == len(align["subject"]["bases"])

    def test_long_alignment_found(self, long_seqs_alignment, frag1, frag2):
        expected_len = len(frag1)
        lens = [a["meta"]["alignment length"] for a in long_seqs_alignment]
        assert expected_len in lens

    def test_expect_max_alignment(self, long_seqs_alignment, frag1, frag2):
        expected_len = len(frag1 + frag2)
        lens = [a["meta"]["alignment length"] for a in long_seqs_alignment]
        assert expected_len == max(lens)

    def test_expected_query_sequence(self, long_seqs_alignment, frag1):
        expected_len = len(frag1)
        passed = False
        for align in long_seqs_alignment:
            if align["meta"]["alignment length"] == expected_len:
                assert align["query"]["bases"].upper() == frag1.upper()
                passed = True
        assert passed

    def test_expected_subject_sequence(self, long_seqs_alignment, frag1):
        expected_len = len(frag1)
        passed = False
        for align in long_seqs_alignment:
            if align["meta"]["alignment length"] == expected_len:
                assert align["subject"]["bases"].upper() == frag1.upper()
                passed = True
        assert passed
