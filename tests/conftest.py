import pytest
import os
from pyblast import Blast


@pytest.fixture(scope="module")
def here():
    return os.path.dirname(os.path.abspath(__file__))


@pytest.fixture(scope="module")
def new_blast(here):
    def make_blast():
        return Blast(
            "db",
            os.path.join(here, "data/test_data/db.fsa"),
            os.path.join(here, "data/test_data/query.fsa"),
            os.path.join(here, "data/blast_results"),
            os.path.join(here, "data/blast_results/results.out"),
        )

    return make_blast
