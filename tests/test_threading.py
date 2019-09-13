from pyblast.blast import BlastBase
import pytest
import os
import time

from multiprocessing import Pool


def run(params):
    blast = params
    blast.makedb()
    blast.quick_blastn()


def test_threading(new_circular_bio_blast, new_primer_blast):

    blasts = [new_circular_bio_blast() for _ in range(10)]
    for b in blasts:
        b.makedb()
    # for b in blasts:
    #     b.quick_blastn()

    with Pool(processes=10) as pool:
        pool.map(run, blasts)
