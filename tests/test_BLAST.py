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


def test_quick_blastn(b):
    b.quick_blastn()
    b.raw_results


def test_aligner(here):
    template_dictionary = os.path.join(here, 'data/test_data/templates')
    query_path = os.path.join(here, 'data/test_data/designs/pmodkan-ho-pact1-z4-er-vpr.gb')
    db_name = 'db'

    a = Aligner(db_name, template_dictionary, query_path)
    a.quick_blastn()

    results = a.results

    for res in results:
        assert res['query_filename'] is not None
        assert res['subject_filename'] is not None
        assert res['query_circular'] is not None
        assert res['subject_circular'] is not None
        assert res['query_name'] is not None
        assert res['subject_name'] is not None





def test_example():
    a = Aligner.use_test_data()
    a.quick_blastn()

    # def test_get_metadata():
    #     a = Aligner.use_test_data()
    #     a.quick_blastn()
    #     a.get_filename(a.input_sequences[0].id)
    #     a.get_is_circular(a.input_sequences[0].id)
