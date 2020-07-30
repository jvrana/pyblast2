from pyblast.cli import is_installed
from pyblast.cli import is_locally_installed


def test_check_installation():
    assert is_installed() or is_locally_installed()
