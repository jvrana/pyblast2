from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pyblast import BioBlast


def compute_length(s, e, l, circular):
    length = e - s
    if length < 0 and not circular:
        raise Exception("not circular")
    elif length < 0 and circular:
        length = length + l + 1
    return length


def location_length(record):
    return compute_length(
        record.features[0].location.start,
        record.features[0].location.end,
        len(record),
        record.annotations["circular"],
    )


def region_length(region):
    return compute_length(
        region["start"], region["end"], region["length"], region["circular"]
    )


def validate_results(results):
    for r in results:
        assert region_length(r["query"]) == region_length(r["subject"])


def validate_blaster_results(blaster):
    for r, alignment in zip(blaster.results, blaster.alignments()):
        assert location_length(alignment[0]) == location_length(alignment[1])


def test_basic_run():
    junk1 = "atgctatgctgatgctgctgtgctgatgctgatgtgtattgctgtatcgcgcgagttagc"
    junk2 = "g" * 30
    frag = "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta"

    query = SeqRecord(seq=Seq(frag), annotations={"circular": False})
    subject = SeqRecord(seq=Seq(junk1 + frag + junk2), annotations={"circular": False})

    # print(type(query))
    # print(type(subject))
    blaster = BioBlast([subject], [query])
    blaster.quick_blastn()
    alignments = blaster.alignments()
    print(alignments)


def test_main(new_bio_blast):
    blast = new_bio_blast()
    blast.quick_blastn()
    validate_blaster_results(blast)
