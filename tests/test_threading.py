from multiprocessing import Pool


def run(params):
    blast = params
    blast.makedb()
    blast.blastn()


def test_threading(new_circular_bio_blast, new_primer_blast):

    blasts = [new_circular_bio_blast() for _ in range(10)]
    for b in blasts:
        b.makedb()
    # for b in blasts:
    #     b.blastn()

    with Pool(processes=10) as pool:
        pool.map(run, blasts)
