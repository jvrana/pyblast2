import pytest
import os


@pytest.fixture(scope="module")
def here():
    return os.path.dirname(os.path.abspath(__file__))
