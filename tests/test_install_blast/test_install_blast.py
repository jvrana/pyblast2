from pyblast.blast_bin.blast_wrapper import BlastWrapper


def test_check_installation():

    b = BlastWrapper()
    if not b.is_installed():
        b.install("justin.vrana@gmail.com", "mac")
    print(b.is_installed())
    # b.install("justin.vrana@gmail.com", "mac")