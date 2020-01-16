import os

from pyblast.blast.blast_parser import BlastResultParser


def test_parser():
    filepath = os.path.join(
        os.path.abspath(os.path.dirname(__file__)), "example_file.txt"
    )

    with open(filepath, "r") as f:
        raw_txt = f.read()
        assert BlastResultParser.raw_results_to_json(raw_txt)
