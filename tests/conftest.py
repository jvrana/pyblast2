import pytest
from os.path import dirname, abspath, join
from pyblast import BlastBase, BioBlast
from pyblast.utils import load_genbank_glob, load_fasta_glob, make_linear, is_circular, make_circular


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

        subjects = load_fasta_glob(join(here, "data/test_data/primers/primers.fasta"))
        subjects = make_linear(subjects)
        queries = load_genbank_glob(
            join(
                here,
                "data/test_data/genbank/designs/pmodkan-ho-pact1-z4-er-vpr_pseudocircular.gb",
            )
        )
        return BioBlast(subjects, queries)

    return make_blast


@pytest.fixture(scope="module")
def new_bio_blast(here):
    def make_blast():

        subjects = load_genbank_glob(
            join(here, "data/test_data/genbank/templates/*.gb")
        )
        queries = load_genbank_glob(
            join(
                here,
                "data/test_data/genbank/designs/pmodkan-ho-pact1-z4-er-vpr_pseudocircular.gb",
            )
        )
        return BioBlast(subjects, queries)

    return make_blast

@pytest.fixture(scope="module")
def new_circular_bio_blast(here):
    def make_blast():

        subjects = load_genbank_glob(
            join(here, "data/test_data/genbank/templates/*.gb")
        )
        queries = load_genbank_glob(
            join(
                here,
                "data/test_data/genbank/designs/pmodkan-ho-pact1-z4-er-vpr.gb",
            )
        )
        queries = make_circular(queries)
        assert is_circular(queries[0])
        return BioBlast(subjects, queries)

    return make_blast
