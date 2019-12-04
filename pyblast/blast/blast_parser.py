"""blast_parser."""
import re


class BlastResultParser:
    """Parses blast results."""

    @staticmethod
    def str_to_f_to_i(v):
        try:
            return int(v)
        except ValueError:
            try:
                return float(v)
            except ValueError:
                pass
        return v

    @staticmethod
    def _extract_metadata(r, delim):
        """Extracts information from the raw text file BLAST produces."""
        g = re.search(
            "#\\s*(?P<blast_ver>.+)\n"
            + "# Query:\\s*(?P<query>.*)\n"
            + "# Database:\\s*(?P<database>.+)\n"
            + r"(?:# Fields:\s*(?P<fields>.+))?",
            r,
        )
        metadata = g.groupdict()
        if metadata["fields"] is None:
            return metadata
        fields_array = re.split(r"\s*{}\s*".format(delim), metadata["fields"])
        metadata["fields"] = fields_array
        return metadata

    @staticmethod
    def _get_alignment_rows(r):
        """Split text into alignment rows."""
        return re.findall("\n([^#].*)", r)

    @classmethod
    def _validate_matches(cls, raw_matches, fields):
        """Create a dictionary from the fields and rows."""
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
            query = dict()
            subject = dict()

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
        """Converts raw BLAST text into a flatten dictionary.

        :param raw_text: raw text from BLAST results
        :type raw_text: basestring
        :param delim: delimiter for parsing
        :type delim: basestring
        :return: flattened dictionary
        :rtype: dict
        """
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
    def get_perfect(data):
        """Returns only exact matches.

        :return:
        :rtype:
        """

        def no_gaps(x):
            return x["meta"]["gaps"] == 0

        def no_gap_opens(x):
            return x["meta"]["gap opens"] == 0

        def identical(x):
            return x["meta"]["identical"] == x["meta"]["alignment length"]

        def perfect(x):
            return all([no_gaps(x), no_gap_opens(x), identical(x)])

        return [r for r in data if perfect(r)]

    @staticmethod
    def get_with_perfect_subjects(data):
        """Returns only parsed alignments with 100% of the subject aligning to
        the query.

        :return: perfect alignments
        :rtype:
        """

        def f(x):
            return x["meta"]["alignment_length"] == x["subject"]["length"]

        return [r for r in data if f(r)]
