import sys
from core import Blast, Aligner
import pytest

def pytest_namespace():
    return {'b':
                Blast('db',
                      'tests/data/test_data/templates',
                      'tests/data/test_data/designs/pmodkan-ho-pact1-z4-er-vpr.gb',
                      'tests/data/blast_results',
                      'tests/data/blast_results/results.out')
            }

@pytest.fixture
def create_blast():
    return  Blast('db',
                      'tests/data/test_data/templates',
                      'tests/data/test_data/designs/pmodkan-ho-pact1-z4-er-vpr.gb',
                      'tests/data/blast_results',
                      'tests/data/blast_results/results.out')

def test_makedb():
    b = create_blast()
    b.makedb()
    b.blastn()
    b.parse_results()

def test_quick_blastn():
    b = create_blast()
    b.quick_blastn()
    b.raw_results

def test_aligner():
    b = create_blast()
    a = Aligner(b.name, b.path_to_input_dir, b.path_to_query)
    a.quick_blastn()
    print(a.raw_results)