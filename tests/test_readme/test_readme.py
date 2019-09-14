def test_example1():
    from pyblast import BioBlast
    from pyblast.utils import make_linear
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    import json

    queries = [
        SeqRecord(
            Seq(
                "ACGTGATTCGTCGTGTAGTTGAGTGTTACGTTGCATGTCGTACGTGTGTAGTGTCGTGTAGTGCTGATGCTACGTGATCG"
            )
        )
    ]
    subjects = [SeqRecord(Seq("ACGTGATTCGTCGTGTAGTTGAGTGTTACGTTGCATGTCGTTACGTGATCG"))]

    # pyblast requires a 'topology' annotation on the SeqRecords.
    # we can make records circular or linear using `make_linear` or `make_circular` methods
    subjects = make_linear(subjects)
    queries = make_linear(queries)

    blast = BioBlast(subjects, queries)
    results = blast.quick_blastn()
    print(json.dumps(results, indent=2))


def test_example2():
    from pyblast import BioBlast
    from pyblast.utils import make_linear, make_circular
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    import json

    seq = "ACGTTGTAGTGTAGTTGATGATGATGTCTGTGTCGTGTGATGTGCTAGGGGTTGATGTGAGTAGTTAGTGGTAGTGTTTAGGGGCGGCGCGGAGTATGCTG"
    queries = [SeqRecord(Seq(seq))]

    subjects = [SeqRecord(Seq(seq[-20:] + seq[:30]))]

    # pyblast requires a 'topology' annotation on the SeqRecords.
    # we can make records circular or linear using `make_linear` or `make_circular` methods
    subjects = make_linear(subjects)
    queries = make_circular(queries)

    blast = BioBlast(subjects, queries)
    results = blast.quick_blastn()
    print(json.dumps(results, indent=2))
