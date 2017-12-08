import sys
from pyblast import Blast, Aligner
import pytest
import os

def pytest_namespace(here):
    return {'b':
                Blast('db',
                      os.path.join(here, 'data/test_data/templates'),
                      os.path.join(here, 'data/test_data/designs/pmodkan-ho-pact1-z4-er-vpr.gb'),
                      os.path.join(here, 'data/blast_results'),
                      os.path.join(here, 'data/blast_results/results.out')
                      )
            }

@pytest.fixture
def b(here):
    return Blast('db',
                      os.path.join(here, 'data/test_data/templates'),
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

def test_aligner(b):
    a = Aligner(b.name, b.path_to_input_dir, b.path_to_query)
    a.quick_blastn()
    print(a.raw_results)

def test_example():
    a = Aligner.use_test_data()
    a.quick_blastn()
    print(a.results)

def test_example():
    a = Aligner.use_test_data()
    a.quick_blastn()

def test_get_metadata():
    a = Aligner.use_test_data()
    a.quick_blastn()
    a.get_filename(a.input_sequences[0].id)
    a.get_is_circular(a.input_sequences[0].id)