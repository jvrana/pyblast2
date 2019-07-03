"""
blast.py

More information: https://www.ncbi.nlm.nih.gov/guide/howto/run-blast-local/
BLAST User Manual: https://www.ncbi.nlm.nih.gov/books/NBK1762/
"""

import json
import os
import re
import tempfile
from uuid import uuid4
from pyblast.utils import run_cmd
from pyblast.utils import json_to_fasta_tempfile, concat_fasta_to_tempfile
from pyblast.results import AlignmentResults
from pyblast.utils.seq_parser import dump_sequence_jsons, load_sequence_jsons
from pyblast.utils import reverse_complement
from pyblast.exceptions import PyBlastException
from pyblast.blast_bin import BlastWrapper
from marshmallow import ValidationError
from copy import copy
from Bio.SeqFeature import FeatureLocation, SeqFeature, CompoundLocation
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import typing
from itertools import chain
from copy import deepcopy


class SeqRecordValidationError(Exception):
    """An error if a Bio.SeqRecord instance is invalid."""


class BlastBase(object):
    """
    A Blast initializer for running blast searches against subjects contained in a directory. Sequences
    are kept on the local machine.

    Note: This is not intended to be the endpoint. You should only use temporary files to run the blast search.
    """

    outfmt = [
        "7",
        "qacc",
        "sacc",
        "score",
        "evalue",
        "bitscore",
        "length",
        "nident",
        "gapopen",
        "gaps",
        "qlen",
        "qstart",
        "qend",
        "slen",
        "sstart",
        "send",
        "sstrand",
        "qseq",
        "sseq",
    ]

    blast_config = {"outfmt": '"{0}"'.format(" ".join(outfmt))}

    def __init__(
        self,
        db_name,
        subject_path,
        query_path,
        db_output_directory,
        results_out_path,
        output_formatter=None,
        **config
    ):
        """
        A Blast initializer for running blast searches.

        :param db_name: Name for database file structure. This name will be appended to all db
        files that blast creates.
        :param subject_path: single fsa file containing the subject sequences
        :param query_path: Location of the fasta or genbank file containing the query
        :param db_output_directory: Location to store database related files
        :param results_out_path: Path to store the results.out file. Path can be absolute or relative.
        :param config: Additional configurations to run for the blast search (see
        https://www.ncbi.nlm.nih.gov/books/NBK279682/)
        """

        # add executable to the system path
        blast_wrapper = BlastWrapper()
        if not blast_wrapper.is_installed():
            blast_wrapper.ask_to_install()

        # name of the database
        self.name = db_name

        self.query_path = query_path
        self.subject_path = subject_path

        # build the configuration
        if output_formatter is None:
            self.outfmt = BlastBase.outfmt
        self.config = {"outfmt": '"{0}"'.format(" ".join(self.outfmt))}
        if config is not None:
            self.config.update(config)

        # path to the directory holding the blast database
        self.db_output_directory = os.path.abspath(db_output_directory)

        # path to the blast database
        self.db = os.path.join(db_output_directory, db_name)

        # list of expected fields in the result
        self.fields = ()

        # the raw results from running blast from the command line
        self.raw_results = None

        # the parsed results as a dictionary
        self.results = None  #

        # the path to the saved results
        self.results_out_path = os.path.abspath(results_out_path)

    def validate_files(self):
        """
        Validate the directories and query files

        :return: None
        :rtype: None
        :raises: PyBlastException if filepaths are invalid or if query sequence has more than one sequence
        """

        def _is_file(myfile):
            return os.path.isfile(os.path.abspath(myfile))

        def _is_dir(mydir):
            return os.path.isdir(os.path.abspath(mydir))

        outdir = os.path.dirname(self.results_out_path)
        errors = []
        for file_ in [self.query_path, self.subject_path]:
            if not _is_file(file_):
                errors.append("File not found: {}".format(file_))
        for dir_ in [outdir, self.db_output_directory]:
            if not _is_dir(dir_):
                errors.append("Directory not found {}".format(dir_))
        if len(errors) > 0:
            raise PyBlastException("\n".join(errors))

        # seq = PySequence.open(self.query_path)
        # if len(seq) == 0:
        #     raise PyBlastException("Query path \"{}\" has no sequences".format(self.query_path))
        # elif len(seq) > 1:
        #     raise PyBlastException("Query path \"{}\" has more than one sequence.".format(self.query_path))

    def create_config(self):
        """Create a configuration dictionary"""
        config_dict = {
            "db": self.db,
            "out": self.results_out_path,
            "query": self.query_path,
        }
        config_dict.update(self.config)
        return config_dict

    def quick_blastn(self):
        """Produce database, run blastn protocol, & parse results """
        self.makedb()
        self.blastn()
        return self.parse_results()

    def quick_blastn_short(self):
        """Produce database, run blastn short protocol, & parse results """
        self.makedb()
        self.blastn_short()
        return self.parse_results()

    def blastn(self):
        return self._run("blastn")

    def blastn_short(self):
        return self._run("blastn", task="blastn-short")

    def _run(self, cmd, **config_opts):
        """Run the blastn using the current configuration"""
        config = self.create_config()
        config.update(config_opts)
        self.run_cmd(cmd, **config)
        with open(self.results_out_path, "rU") as handle:
            self.raw_results = handle.read()
        return self.raw_results

    # Wrapper for the util.run_cmd
    def run_cmd(self, cmd, **kwargs):
        """Wrapper for utils.run_cmd"""
        run_cmd(cmd, **kwargs)

    # def concat_templates(self):
    #     """
    #     Gathers all of the sequences in the input dir and concatenates them into a single fasta file
    #
    #     :return: [concatenated fasta file, sequences used to make fassta file, metadata for the sequences
    #     :rtype: list
    #     """
    #     out = self.db + '.fsa'
    #     seqs = PySequence.concat_seqs(self.path_to_input_dir, out)
    #
    #     self.input_sequences = seqs
    #     return out, seqs

    def _fasta_details(self, path):
        seqs = list(SeqIO.parse(path, format="fasta"))
        tot_bps = sum([len(s) for s in seqs])
        return {"num_sequence": len(seqs), "total_bps": tot_bps}

    def makedb(self, verbose=False):
        """Creates a blastdb from sequences grabbed from the input directory"""
        self.validate_files()

        if verbose:
            query_details = self._fasta_details(self.query_path)
            subject_details = self._fasta_details(self.subject_path)

            print("Making BLAST database from:")
            print(json.dumps(subject_details, indent=2))

            print("Query:")
            print(json.dumps(query_details, indent=2))

        self.run_cmd(
            "makeblastdb",
            dbtype="nucl",
            title=self.name,
            out=self.db,
            **{"in": self.subject_path}
        )
        return self.db

    def parse_results(self, delim=","):
        """
        Parses the raw blast result to a JSON

        :param save_as_json: whether to save the JSON. Results will be saved in same directory and name as
        results_out_path but with a .json extension.
        :type save_as_json: bool
        :param delim: delimiter to parse
        :type delim: str
        :return:
        :rtype:
        """
        self.results = BlastParser.results_to_json(self.raw_results, delim=delim)
        return self.results

    def get_perfect(self):
        return BlastParser.get_perfect(self.results)

    def __str__(self):
        return "{}".format(self.create_config())


class Aligner(BlastBase):
    """
    A Blast object that stores the database files in a hidden temporary directory. Use
    "quick_blastn" for returning results as a python object.
    """

    def __init__(self, db_name, subject_path, query_path, **config):
        """
        Aligner: A Blast object that stores the database files in a hidden temporary directory.
        :param db_name: Name for database file structure. This name will be appended to all db files that blast creates.
        :param subject_path: Input directory or file containing a list of subjects to align against the query
        :param query_path: Location of the fasta or genbank file containing the query
        :param config: Additional configurations to run for the blast search (see
        https://www.ncbi.nlm.nih.gov/books/NBK279682/)
        """

        db_output_directory = tempfile.mkdtemp()
        fv, out = tempfile.mkstemp(dir=db_output_directory)

        if os.path.isdir(subject_path):
            subject_path = concat_fasta_to_tempfile(subject_path)

        super(Aligner, self).__init__(
            db_name=db_name,
            subject_path=subject_path,
            query_path=query_path,
            db_output_directory=db_output_directory,
            results_out_path=out,
            **config
        )

        self.fields = self.fields + (
            "subject_name",
            "subject_filename",
            "subject_circular",
            "query_name",
            "query_filename",
            "query_circular",
        )

    @classmethod
    def use_test_data(cls):
        """Create a Blast instance using predefined data located in tests"""
        dir_path = os.path.dirname(os.path.realpath(__file__))

        return cls(
            "db",
            os.path.join(dir_path, "..", "tests/data/test_data/db.fsa"),
            os.path.join(dir_path, "..", "tests/data/test_data/query.fsa"),
        )

    def makedb(self):
        super(Aligner, self).makedb()


class Constants(object):

    PARENT = "origin_record"
    CIRCULAR = "circular"
    LINEAR = "linear"
    TRANSFORMATION = "modification"
    OLD_KEY = "original_record_id"
    PSEUDOCIRCULAR = "pseudocircularized"
    COPY_RECORD = "copy_record"

    class NULL(object):
        pass


class BlastParser(object):
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
    def __clean_json(data_list):
        for data in data_list:
            query = {}
            subject = {}

            query["start"] = data["q. start"]
            query["end"] = data["q. end"]
            query["bases"] = data["query seq"]
            query["strand"] = data.get("query strand", "plus")
            query["length"] = data["query length"]
            query["sequence_id"] = data["query acc."]

            subject["start"] = data["s. start"]
            subject["end"] = data["s. end"]
            subject["bases"] = data["subject seq"]
            subject["strand"] = data.get("subject strand", "plus")
            subject["length"] = data["subject length"]
            subject["sequence_id"] = data["subject acc."]

            meta = dict(data)
            yield {"query": query, "subject": subject, "meta": meta}

    @classmethod
    def results_to_json(cls, raw_text, delim=","):
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
        data = list(cls.__clean_json(match_dicts))
        return data

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
                    location=FeatureLocation(
                        align["query"]["start"], align["query"]["end"], strand=1
                    ),
                    type="alignment",
                    id="query alignment",
                )
            ],
            annotations=dict(align),
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
                    location=FeatureLocation(
                        align["subject"]["start"],
                        align["subject"]["end"],
                        strand=strand,
                    ),
                    type="alignment",
                    id="subject alignment",
                )
            ],
            annotations=dict(align),
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


class SeqRecordDB(object):
    def __init__(self):
        self.records = {}
        self.mapping = {}

    def incr(self):
        return str(uuid4())

    @classmethod
    def validate_records(cls, records):
        for r in records:
            cls.validate_circular(r)

    @staticmethod
    def validate_circular(r):
        if (
            Constants.LINEAR not in r.annotations
            and Constants.CIRCULAR not in r.annotations
        ):
            raise SeqRecordValidationError(
                "SeqRecord {} is missing a '{linear}' or '{circular}' annotation. This must be provided.".format(
                    r, circular=Constants.CIRCULAR, linear=Constants.LINEAR
                )
            )
        if Constants.LINEAR in r.annotations and Constants.CIRCULAR in r.annotations:
            if r.annotations[Constants.LINEAR] is r.annotations[Constants.CIRCULAR]:
                raise SeqRecordValidationError(
                    "SeqRecord {} has conflicting topology definitions for '{linear}' and '{circular}'".format(
                        r, circular=Constants.CIRCULAR, linear=Constants.LINEAR
                    )
                )

    @staticmethod
    def is_circular(r):
        return (
            r.annotations.get(Constants.CIRCULAR, False) is True
            or r.annotations.get(Constants.LINEAR, True) is False
        )

    def key(self, record):
        k = self.mapping.get(id(record), None)
        return k

    def post_transform_hook(self, key, *args, **kwargs):
        record = self.get(key)
        record.annotations[Constants.OLD_KEY] = record.id
        record.id = key
        return record

    def transform(self, key, transform, transform_label):
        record = self.get(key)
        new_record = transform(record)
        parent_key = self.add_one(record)
        new_key = self.add_one(new_record)
        new_record.annotations[Constants.PARENT] = parent_key
        new_record.annotations[Constants.TRANSFORMATION] = transform_label
        self.post_transform_hook(new_key)
        return new_key

    def add_with_transformation(self, record, transform, transform_label):
        key = self.add_one(record)
        return self.transform(key, transform, transform_label)

    def add_many_with_transformations(self, records, transform, transform_label):
        return [
            self.add_with_transformation(r, transform, transform_label) for r in records
        ]

    def get_many(self, keys, default=Constants.NULL):
        records = []
        for k in keys:
            records.append(self.get(k, default=default))
        return records

    def get(self, key, default=Constants.NULL):
        if default is not Constants.NULL:
            return self.records.get(key, default)
        return self.records[key]

    def get_origin(self, key, blacklist=None):
        r = self.get(key)
        if r:
            if Constants.PARENT in r.annotations:
                if r.annotations[Constants.TRANSFORMATION] in blacklist:
                    return r
                else:
                    return self.get_origin(r.annotations[Constants.PARENT])
        return r

    def add_one(self, record, validate=True):
        if self.key(record):
            return self.key(record)
        else:
            if validate:
                self.validate_records([record])
            key = self.incr()
            self.records[key] = record
            self.mapping[id(record)] = key
            return key

    def add_many(self, records):
        self.validate_records(records)
        keys = []
        for r in records:
            keys.append(self.add_one(r, validate=False))
        return keys


# TODO: """This should be the only entrypoint to the program."""
class SeqRecordBlast(Aligner):
    def __init__(
        self,
        subjects: typing.Sequence[SeqRecord],
        queries: typing.Sequence[SeqRecord],
        span_origin=False,
        **config
    ):
        self.seq_db = SeqRecordDB()
        subjects = self.add_records(subjects)
        queries = self.add_records(queries)
        self.subjects = subjects
        self.queries = queries
        self.__span_origin = span_origin
        db_name = str(uuid4())
        subject_path = self.tmp_fasta(subjects)
        query_path = self.tmp_fasta(queries)
        super().__init__(
            db_name=db_name, subject_path=subject_path, query_path=query_path, **config
        )

    def add_records(self, records):
        def copy_record(r):
            return deepcopy(r)

        def pseudocircularize(r):
            r2 = r + r
            r2.name = Constants.PSEUDOCIRCULAR + "__" + r.name
            r2.id = str(uuid4())
            return r2

        circular = [r for r in records if self.seq_db.is_circular(r)]
        linear = [r for r in records if not self.seq_db.is_circular(r)]
        keys = self.seq_db.add_many_with_transformations(
            circular, pseudocircularize, Constants.PSEUDOCIRCULAR
        )
        keys += self.seq_db.add_many_with_transformations(
            linear, copy_record, Constants.COPY_RECORD
        )
        return self.seq_db.get_many(keys)

    def parse_results(self, delim=","):
        parsed_results = BlastParser.results_to_json(self.raw_results, delim=delim)

        # TODO: resolve with sequence dictionary, resolving pseudocircularized constructs
        for v in parsed_results:
            for x in ["query", "subject"]:
                record = self.seq_db.get_origin(
                    v[x]["sequence_id"], blacklist=[Constants.COPY_RECORD]
                )
                v[x]["circular"] = self.seq_db.is_circular(record)
                v[x]["name"] = record.name
                v[x]["origin_sequence_id"] = record.id
        self.results = parsed_results
        return self.results

    @staticmethod
    def tmp_fasta(records):
        fd, tmp_path_handle = tempfile.mkstemp(suffix=".fasta")
        SeqIO.write(records, tmp_path_handle, format="fasta")
        return tmp_path_handle

    def alignments(self):
        results_json = BlastParser.results_to_json(self.raw_results)
        alignments = []
        for align in results_json:
            alignments.append(
                BlastParser.alignment_to_seqrecord(align, self.seq_db.records)
            )
        return alignments


class JSONParser(object):
    @staticmethod
    def clean_data(data):
        return {k: v for k, v in data.items() if v is not None}

    @staticmethod
    def JSON_to_SeqFeature(data, length=None):
        start = data["start"]
        end = data["end"]
        if start > end:
            if length is None:
                raise ValueError(
                    "A length must be provided to create a feature with start > end."
                )
            pass
            f1 = FeatureLocation(start, length, data["strand"])
            f2 = FeatureLocation(0, end, data["strand"])
            if data["strand"] == -1:
                location = CompoundLocation([f2, f1])
            else:
                location = CompoundLocation([f1, f2])
        else:
            location = FeatureLocation(
                data["start"], data["end"], strand=data["strand"]
            )
        return SeqFeature(location=location, type=data["type"], id=data["name"])

    @classmethod
    def __JSON_to_SeqRecord(cls, data):
        data = cls.clean_data(data)
        return SeqRecord(
            seq=Seq(data["bases"]),
            annotations={"circular": data.get("circular", False)},
            features=[
                cls.JSON_to_SeqFeature(f, len(data["bases"]))
                for f in data.get("features", [])
            ],
            description=data.get("description", ""),
            name=data["name"],
            id=data.get("id", str(uuid4())),
        )

    @classmethod
    def JSON_to_SeqRecords(cls, datalist):
        records = []
        if isinstance(datalist, list):
            for data in datalist:
                records.append(cls.__JSON_to_SeqRecord(data))
        else:
            records.append(cls.__JSON_to_SeqRecord(datalist))
        return records


class JSONBlast(SeqRecordBlast):
    """Object that runs blast starting from JSON inputs and outputs"""

    def __init__(self, subject_json, query_json, **config):
        """
        Initialize JSONBlast

        :param subject_json: subject sequences as serialized json. Schema may be found in :class:`pyblast.schema.SequenceSchema`
        :type subject_json: list of dict, mapping, or Object
        :param query_json: query sequence as serialized json
        :type query_json: dict, mapping, or Object
        :param preloaded: whether the data is preloaded into the SequenceSchema or not. Practially, this means the "bases"
        for the sequence if refered to as "sequence" in serialized data and as "bases" after data is loaded (preloaded=True)
        :param span_origin: default False. If True, circular sequences will be pseudocircularized for alignment over the origin
        :type span_origin: boolean
        :param config: optional arguments
        :type config: dict
        """
        srecords = JSONParser.JSON_to_SeqRecords(subject_json)
        qrecords = JSONParser.JSON_to_SeqRecords(query_json)
        self.subject_json = subject_json
        self.query_json = query_json
        super().__init__(srecords, qrecords)

    @classmethod
    def use_test_data(cls):
        """Create a Blast instance using predefined data located in tests"""
        dir_path = os.path.dirname(os.path.realpath(__file__))
        with open(
            os.path.join(dir_path, "..", "tests/data/test_data/templates.json"), "r"
        ) as f:
            subject = json.load(f)
        with open(
            os.path.join(dir_path, "..", "tests/data/test_data/query.json"), "r"
        ) as f:
            query = json.load(f)
        return cls(subject_json=subject, query_json=query)
