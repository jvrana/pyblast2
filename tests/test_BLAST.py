import sys
from core import BLAST

def test_makedb():
    b = BLAST('db',
              'tests/data/test_data/templates',
              'tests/data/test_data/designs/pmodkan-ho-pact1-z4-er-vpr.gb',
              'tests/data/blast_results',
              'tests/data/blast_results/results.out')
    b.makedb()
    print(b)
    b.run()

    import re
    r = b.results