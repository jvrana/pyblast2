from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from uuid import uuid4
import re
from pyblast.utils import new_feature_location


class BlastResultParser(object):
    @staticmethod
    def str_to_f_to_i(v):
        """"""
        try:
            v = float(v)
        except ValueError:
            pass
        try:
            v = int(v)
        except ValueError:
            pass
        return v

    @staticmethod
    def _extract_metadata(r, delim):
        """Extracts information from the raw text file BLAST produces"""
        g = re.search(
            "#\s*(?P<blast_ver>.+)\n"
            + "# Query:\s*(?P<query>.*)\n"
            + "# Database:\s*(?P<database>.+)\n"
            + "(?:# Fields:\s*(?P<fields>.+))?",
            r,
        )
        metadata = g.groupdict()
        if metadata["fields"] is None:
            return metadata
        fields_array = re.split("\s*{}\s*".format(delim), metadata["fields"])
        metadata["fields"] = fields_array
        return metadata

    @staticmethod
    def _get_alignment_rows(r):
        """Split text into alignment rows"""
        return re.findall("\n([^#].*)", r)

    @classmethod
    def _validate_matches(cls, raw_matches, fields):
        """Create a dictionary from the fields and rows"""
        match_dicts = []
        for m in raw_matches:
            values = [cls.str_to_f_to_i(v) for v in m.split("\t")]
            match_dicts.append(dict(list(zip(fields, values))))
        return match_dicts

    @staticmethod
    def __convert_strand_label(strand_lable):
        if strand_lable.strip().lower() != "plus":
            return -1
        return 1

    @classmethod
    def __clean_json(cls, data_list):
        for data in data_list:
            query = {}
            subject = {}

            query["start"] = data["q. start"]
            query["end"] = data["q. end"]
            query["bases"] = data["query seq"]
            query["strand"] = cls.__convert_strand_label(
                data.get("query strand", "plus")
            )
            query["length"] = data["query length"]
            query["sequence_id"] = data["query acc."]

            subject["start"] = data["s. start"]
            subject["end"] = data["s. end"]
            subject["bases"] = data["subject seq"]
            subject["strand"] = cls.__convert_strand_label(data["subject strand"])
            subject["length"] = data["subject length"]
            subject["sequence_id"] = data["subject acc."]

            meta = dict(data)
            yield {"query": query, "subject": subject, "meta": meta}

    @classmethod
    def raw_results_to_json(cls, raw_text, delim=","):
        """
        Converts raw BLAST text into a flatten dictionary

        :param raw_text: raw text from BLAST results
        :type raw_text: basestring
        :param delim: delimiter for parsing
        :type delim: basestring
        :return: flattened dictionary
        :rtype: dict
        """

        # print(results)
        if raw_text.strip() == "":
            return []
        meta = cls._extract_metadata(raw_text, delim)
        fields = meta["fields"]
        if fields is None:
            return []
        alignment_rows = cls._get_alignment_rows(raw_text)
        match_dicts = cls._validate_matches(alignment_rows, tuple(fields))
        data = list(cls.__clean_json(match_dicts))
        return data

    @staticmethod
    def alignment_to_feature(align: dict, type: str, id: str):
        return SeqFeature(
            location=new_feature_location(
                align["start"] - 1,
                align["end"],
                strand=align["strand"],
                length=align["length"],
            ),
            type=type,
            id=id,
        )

    @classmethod
    def alignment_to_seq_record(cls, align, name, id, seq_dict: dict):

        origin_record = seq_dict[align["sequence_id"]]

        new_record = SeqRecord(
            id=id,
            seq=origin_record.seq,
            name=name,
            features=[cls.alignment_to_feature(align, "alignment", id)],
            annotations={},
        )

    @staticmethod
    def alignment_to_seqrecord(align: dict, seq_dict: dict):

        query_record = seq_dict[align["query"]["sequence_id"]]
        subject_record = seq_dict[align["subject"]["sequence_id"]]

        query_alignment = SeqRecord(
            seq=query_record.seq,
            id=str(uuid4()),
            name="QUERY__{}__alignedto__{}".format(
                subject_record.name, query_record.name
            ),
            features=[
                SeqFeature(
                    location=new_feature_location(
                        align["query"]["start"] - 1,
                        align["query"]["end"],
                        strand=1,
                        length=align["query"]["length"],
                    ),
                    type="alignment",
                    id="query alignment",
                )
            ],
            annotations={
                "alignment": align,
                "circular": query_record.annotations["circular"],
            },
        )
        if align["subject"]["strand"] != "plus":
            strand = -1
        else:
            strand = 1
        subject_alignment = SeqRecord(
            seq=subject_record.seq,
            id=str(uuid4()),
            name="SUBJECT__{}__alignedto__{}".format(
                subject_record.name, query_record.name
            ),
            features=[
                SeqFeature(
                    location=new_feature_location(
                        align["subject"]["start"] - 1,
                        align["subject"]["end"],
                        strand=strand,
                        length=align["subject"]["length"],
                    ),
                    type="alignment",
                    id="subject alignment",
                )
            ],
            annotations={
                "alignment": align,
                "circular": subject_record.annotations["circular"],
            },
        )

        return query_alignment, subject_alignment

    @staticmethod
    def get_perfect(data):
        """
        Returns only exact matches.

        :return:
        :rtype:
        """
        no_gaps = lambda x: x["meta"]["gaps"] == 0
        no_gap_opens = lambda x: x["meta"]["gap opens"] == 0
        identical = lambda x: x["meta"]["identical"] == x["meta"]["alignment length"]
        perfect = lambda x: all([no_gaps(x), no_gap_opens(x), identical(x)])
        return [r for r in data if perfect(r)]

    @staticmethod
    def get_with_perfect_subjects(data):
        """
        Returns only parsed alignments with 100% of the subject aligning to the query

        :return: perfect alignments
        :rtype:
        """
        f = lambda x: x["meta"]["alignment_length"] == x["subject"]["length"]
        return [r for r in data if f(r)]
