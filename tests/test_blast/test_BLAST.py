import os
import pytest

from pyblast import BlastBase, TmpBlast


def pytest_namespace(here):
    return {
        "b": BlastBase(
            "db",
            os.path.join(here, "data/test_data/db.fsa"),
            os.path.join(here, "data/test_data/query.fsa"),
            os.path.join(here, "data/blast_results"),
            os.path.join(here, "data/blast_results/results.out"),
        )
    }


@pytest.fixture
def b(here):
    return BlastBase(
        "db",
        os.path.join(here, "data/test_data/db.fsa"),
        os.path.join(here, "data/test_data/query.fsa"),
        os.path.join(here, "data/blast_results"),
        os.path.join(here, "data/blast_results/results.out"),
    )


def test_makedb(b):
    b.makedb()


def test_blasn(b):
    b.makedb()
    b.blastn()


def test_parse_results(b):
    b.makedb()
    b.blastn()
    b.parse_results()

    results = b.results
    res = results[0]
    assert "query" in res
    assert "subject" in res
    assert "meta" in res

    query = res["query"]
    assert "sequence_id" in query
    assert "start" in query
    assert "end" in query
    assert "length" in query

    subject = res["subject"]
    assert "sequence_id" in subject
    assert "start" in subject
    assert "end" in subject
    assert "length" in subject
    assert subject["strand"] in [1, -1]

    meta = res["meta"]
    assert "score" in meta
    assert "evalue" in meta
    assert "bit score" in meta
    assert "identical" in meta
    assert "gap opens" in meta
    assert "gaps" in meta
    assert "alignment length" in meta


def test_quick_blastn(b):
    b.quick_blastn()
    results = b.raw_results


class TestAligner:
    @pytest.fixture
    def aligner(self, here):
        template_dictionary = os.path.join(here, "data/test_data/db.fsa")
        query_path = os.path.join(here, "data/test_data/query.fsa")
        db_name = "db"

        a = TmpBlast(db_name, template_dictionary, query_path)
        a.quick_blastn()
        return a

    def test_query(self, aligner):
        results = aligner.results

    def test_subject(self, aligner):
        results = aligner.results
        assert len(results)
        for res in results:
            assert res["subject"]["strand"] in [1, -1]

    def test_example(self):
        a = TmpBlast.use_test_data()
        a.quick_blastn()
        print(a.results[0])

    # def test_get_metadata():
    #     a = TmpBlast.use_test_data()
    #     a.quick_blastn()
    #     a.get_filename(a.input_sequences[0].id)
    #     a.get_is_circular(a.input_sequences[0].id)
