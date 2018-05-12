import os

import pytest

from pyblast import Blast, Aligner


def pytest_namespace(here):
    return {'b':
                Blast('db',
                      os.path.join(here, 'data/test_data/db.fsa'),
                      os.path.join(here, 'data/test_data/designs/pmodkan-ho-pact1-z4-er-vpr.gb'),
                      os.path.join(here, 'data/blast_results'),
                      os.path.join(here, 'data/blast_results/results.out')
                      )
            }


@pytest.fixture
def b(here):
    return Blast('db',
                 os.path.join(here, 'data/test_data/db.fsa'),
                 os.path.join(here, 'data/test_data/designs/pmodkan-ho-pact1-z4-er-vpr.gb'),
                 os.path.join(here, 'data/blast_results'),
                 os.path.join(here, 'data/blast_results/results.out')
                 )


def test_makedb(b):
    b.makedb()
    b.blastn()
    b.parse_results()

    results = b.results
    res = results[0]
    assert 'query' in res
    assert 'subject' in res
    assert 'meta' in res

    query = res['query']
    assert 'acc' in query
    assert 'start' in query
    assert 'end' in query
    assert 'length' in query
    assert query['name'] is None
    assert query['filename'] is None
    assert query['circular'] is None

    subject = res['subject']
    assert 'acc' in subject
    assert 'start' in subject
    assert 'end' in subject
    assert 'length' in subject
    assert subject['name'] is None
    assert subject['filename'] is None
    assert subject['circular'] is None
    assert subject['strand'] in ['plus', 'minus']

    meta = res['meta']
    assert 'score' in meta
    assert 'evalue' in meta
    assert 'bit_score' in meta
    assert 'identical' in meta
    assert 'gaps_open' in meta
    assert 'gaps' in meta
    assert 'alignment_length' in meta


def test_quick_blastn(b):
    b.quick_blastn()
    results = b.raw_results


class TestAligner:

    @pytest.fixture
    def aligner(self, here):
        template_dictionary = os.path.join(here, 'data/test_data/templates')
        query_path = os.path.join(here, 'data/test_data/designs/pmodkan-ho-pact1-z4-er-vpr.gb')
        db_name = 'db'

        a = Aligner(db_name, template_dictionary, query_path)
        a.quick_blastn()
        return a

    def test_query(self, aligner):
        results = aligner.results

        for res in results:
            assert res['query']['filename'] is not None
            assert res['query']['circular'] is not None
            assert res['query']['name'] is not None

    def test_subject(self, aligner):
        results = aligner.results
        for res in results:
            assert res['subject']['filename'] is not None
            assert res['subject']['circular'] is not None
            assert res['subject']['name'] is not None
            assert res['subject']['strand'] in ['plus', 'minus']

    def test_example(self):
        a = Aligner.use_test_data()
        a.quick_blastn()

    def test_perfect_matches(self):
        a = Aligner.use_test_data()
        a.find_perfect_matches()
    # def test_get_metadata():
    #     a = Aligner.use_test_data()
    #     a.quick_blastn()
    #     a.get_filename(a.input_sequences[0].id)
    #     a.get_is_circular(a.input_sequences[0].id)
