from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pyblast import BioBlast


class Validator(object):
    @staticmethod
    def compute_length(s, e, l, circular):
        length = e - s
        if length < 0 and not circular:
            raise Exception("not circular")
        elif length < 0 and circular:
            length = length + l + 1
        return length

    @classmethod
    def location_length(cls, record):
        return cls.compute_length(
            record.features[0].location.start,
            record.features[0].location.end,
            len(record),
            record.annotations["circular"],
        )

    @classmethod
    def region_length(cls, region):
        return cls.compute_length(
            region["start"], region["end"], region["length"], region["circular"]
        )

    @classmethod
    def validate_blaster_results(cls, blaster):
        errors = []
        for r, alignment in zip(blaster.results, blaster.alignments()):
            if not cls.location_length(alignment[0]) == cls.location_length(
                alignment[1]
            ):
                errors.append(
                    (
                        "Alignments have different lengths: {}".format(alignment),
                        alignment,
                    )
                )
            if not cls.region_length(r["query"]) == cls.region_length(r["subject"]):
                errors.append(("Regions have different lengths: {}".format(r), r))
        assert not errors


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


def test_valid_results(new_bio_blast):
    blast = new_bio_blast()
    blast.quick_blastn()
    Validator.validate_blaster_results(blast)
