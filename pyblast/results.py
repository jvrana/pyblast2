"""AlignmentResults"""

import re
import json
from pyblast.schema import (
    QuerySchema,
    SubjectSchema,
    AlignmentSchema,
    AlignmentMetaSchema,
)


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


class AlignmentResults(object):
    """Alignments results container"""

    def __init__(self, alignments):
        self._alignments = tuple(alignments)

    def __len__(self):
        return len(self._alignments)

    @property
    def alignments(self):
        return self._alignments

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
            values = [str_to_f_to_i(v) for v in m.split("\t")]
            match_dicts.append(dict(list(zip(fields, values))))
        return match_dicts

    @classmethod
    def _cleanup_raw_results(cls, raw_text, delim):
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
            return {}
        meta = cls._extract_metadata(raw_text, delim)
        fields = meta["fields"]
        if fields is None:
            return [{}]
        alignment_rows = cls._get_alignment_rows(raw_text)
        match_dicts = cls._validate_matches(alignment_rows, tuple(fields))
        return match_dicts

    @staticmethod
    def _serialize_data(data, context):
        """
        Serializes a flat dictionary of data into Query, Subject, AlginmentMeta, and Alignment models

        :param data: flatten dictionary
        :type data: dict
        :return: serialized data
        :rtype: dict
        """

        query_schema = QuerySchema(many=True)
        subject_schema = SubjectSchema(many=True)
        meta_schema = AlignmentMetaSchema(many=True)
        alignment_schema = AlignmentSchema(many=True)
        alignment_schema.context = context

        # load marshall alignment
        queries = query_schema.load(data)
        subjects = subject_schema.load(data)
        metas = meta_schema.load(data)

        serialized = alignment_schema.dump(
            [
                {"query": q, "subject": s, "meta": m}
                for q, s, m in zip(queries, subjects, metas)
            ]
        )
        return serialized

    @classmethod
    def parse_results(cls, data, delim, context=None):
        """Convert raw BLAST results into Alignment model"""
        # first clean up the data and get the flattened dictionary
        cleaned = cls._cleanup_raw_results(data, delim)
        if cleaned == [{}]:
            return cls(alignments=list())
        # force data to deserialized according to the schemas
        # if context is None:
        #     context = {'db': {}}
        serialized = cls._serialize_data(cleaned, context)

        # force data back into models
        alignment_schema = AlignmentSchema(many=True)
        alignment_schema.context = context
        deserialized = alignment_schema.load(serialized)

        return cls(alignments=deserialized)

    def get_perfect(self):
        """
        Returns only exact matches.

        :return:
        :rtype:
        """
        no_gaps = lambda x: x["meta"]["gaps"] == 0
        no_gap_opens = lambda x: x["meta"]["gaps_open"] == 0
        identical = lambda x: x["meta"]["identical"] == x["meta"]["alignment_length"]
        perfect = lambda x: all([no_gaps(x), no_gap_opens(x), identical(x)])
        return self.__class__([r for r in self.alignments if perfect(r)])

    def get_with_perfect_subjects(self):
        """
        Returns only parsed alignments with 100% of the subject aligning to the query

        :return: perfect alignments
        :rtype:
        """
        f = lambda x: x["meta"]["alignment_length"] == x["subject"]["length"]
        return self.__class__([r for r in self.alignments if f(r)])

    def dump_to_json(self, path):
        with open(path, "w") as out:
            json.dump(self.alignments, out)
