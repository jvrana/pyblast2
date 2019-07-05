from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pyblast import BioBlast
from lobio.models import Region, Context


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

    @staticmethod
    def blast_result_to_region(data, gaps=0):
        strand = data["strand"]
        s, e = data["start"], data["end"]
        e += gaps
        if strand == -1:
            s, e = e, s
        r1 = Region(
            s,
            e,
            Context(data["length"] + gaps, data.get("circular", False), start_index=1),
            strand,
        )
        return r1

    @classmethod
    def validate_blaster_results(cls, blaster):
        errors = []
        for r in blaster.results:
            q = r["query"]
            s = r["subject"]
            try:
                gaps = r["meta"]["gaps"]
                r1 = cls.blast_result_to_region(q, 0)
                r2 = cls.blast_result_to_region(s, gaps)
                l1 = len(r1)
                l2 = len(r2)
                if not l1 == l2:
                    errors.append(
                        ("Regions have different lengths {} {}".format(l1, l2), q, s)
                    )
            except Exception as e:
                errors.append((str(e), q, s))
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


def test_perfect_match():

    s = "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta"
    r1 = SeqRecord(seq=Seq(s), annotations={"topology": "linear"})
    r2 = SeqRecord(seq=Seq(s), annotations={"topology": "linear"})

    blaster = BioBlast([r1], [r2])
    blaster.quick_blastn()
    result = blaster.results[0]
    assert result["query"]["start"] == 1
    assert result["query"]["end"] == len(s)
    assert result["subject"]["start"] == 1
    assert result["subject"]["end"] == len(s)


def test_valid_results(new_bio_blast):
    blast = new_bio_blast()
    blast.quick_blastn()
    Validator.validate_blaster_results(blast)
