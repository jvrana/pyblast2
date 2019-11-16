from os.path import abspath
from os.path import dirname
from os.path import join

import pytest

from pyblast import BioBlast
from pyblast.blast import BlastBase
from pyblast.utils import is_circular
from pyblast.utils import load_fasta_glob
from pyblast.utils import load_genbank_glob
from pyblast.utils import make_circular
from pyblast.utils import make_linear


@pytest.fixture(scope="module")
def here():
    return dirname(abspath(__file__))


@pytest.fixture(scope="module")
def new_blast(here):
    def make_blast():
        return BlastBase(
            "db",
            join(here, "data/test_data/db.fsa"),
            join(here, "data/test_data/query.fsa"),
            join(here, "data/blast_results"),
            join(here, "data/blast_results/results.out"),
        )

    return make_blast


@pytest.fixture(scope="module")
def new_primer_blast(here):
    def make_blast():

        subjects = load_fasta_glob(
            join(here, "data/test_data/primers/primers.fasta"), force_unique_ids=True
        )
        subjects = make_linear(subjects)
        queries = load_genbank_glob(
            join(here, "data/test_data/genbank/designs/pmodkan-ho-pact1-z4-er-vpr.gb"),
            force_unique_ids=True,
        )
        return BioBlast(subjects, queries)

    return make_blast


@pytest.fixture(scope="module")
def new_circular_bio_blast(here):
    def make_blast():

        subjects = load_genbank_glob(
            join(here, "data/test_data/genbank/templates/*.gb"), force_unique_ids=True
        )
        queries = load_genbank_glob(
            join(here, "data/test_data/genbank/designs/pmodkan-ho-pact1-z4-er-vpr.gb"),
            force_unique_ids=True,
        )
        queries = make_circular(queries)
        assert is_circular(queries[0])
        return BioBlast(subjects, queries)

    return make_blast
