from pyblast.blast import BlastBase
import pytest
import os

from multiprocessing import Process


def test_threading(new_blast):
    def run_blast():
        b = new_blast()
        b.quick_blastn()
        return b

    p = Process(target=run_blast)
    p.start()
