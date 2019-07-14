"""
blast.py

More information: https://www.ncbi.nlm.nih.gov/guide/howto/run-blast-local/
BLAST User Manual: https://www.ncbi.nlm.nih.gov/books/NBK1762/
"""

import json
import os
import tempfile
from uuid import uuid4
from pyblast.utils import run_cmd
from pyblast.exceptions import PyBlastException
from pyblast.blast_bin import BlastWrapper
from Bio.SeqRecord import SeqRecord
import typing
from copy import deepcopy
from Bio import SeqIO
from pyblast.constants import Constants as C
from .blast_parser import BlastResultParser
from .json_parser import JSONParser
from pyblast.utils import glob_fasta_to_tmpfile, records_to_tmpfile, clean_records
from .seqdb import SeqRecordDB
from more_itertools import unique_everseen


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
        **additional_config
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
        self._config = {}

        # build the configuration
        if output_formatter is None:
            self.outfmt = BlastBase.outfmt
        else:
            self.outfmt = output_formatter

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

    @property
    def config(self):
        out_config = {"outfmt": '"{0}"'.format(" ".join(self.outfmt))}
        out_config.update(self._config)
        return out_config

    def update_config(self, config):
        self._config.update(config)

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

    def _unique_results(self, results):
        return list(
            unique_everseen(
                results,
                key=lambda x: (
                    x["query"]["sequence_id"],
                    x["query"]["start"],
                    x["query"]["end"],
                    x["subject"]["sequence_id"],
                    x["subject"]["start"],
                    x["subject"]["end"],
                ),
            )
        )

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
        results = BlastResultParser.raw_results_to_json(self.raw_results, delim=delim)
        results = self._unique_results(results)
        self.results = results
        return self.results

    def get_perfect(self):
        return BlastResultParser.get_perfect(self.results)

    def __str__(self):
        return "{}".format(self.create_config())


class TmpBlast(BlastBase):
    """
    A Blast object that stores the database files in a hidden temporary directory. Use
    "quick_blastn" for returning results as a python object.
    """

    def __init__(self, db_name, subject_path, query_path, **config):
        """
        TmpBlast: A Blast object that stores the database files in a hidden temporary directory.
        :param db_name: Name for database file structure. This name will be appended to all db files that blast creates.
        :param subject_path: Input directory or file containing a list of subjects to align against the query
        :param query_path: Location of the fasta or genbank file containing the query
        :param config: Additional configurations to run for the blast search (see
        https://www.ncbi.nlm.nih.gov/books/NBK279682/)
        """

        db_output_directory = tempfile.mkdtemp()
        fv, out = tempfile.mkstemp(dir=db_output_directory)

        if os.path.isdir(subject_path):
            subject_path = glob_fasta_to_tmpfile(subject_path)

        super(TmpBlast, self).__init__(
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
            os.path.join(dir_path, "..", "..", "tests/data/test_data/db.fsa"),
            os.path.join(dir_path, "..", "..", "tests/data/test_data/query.fsa"),
        )

    def makedb(self):
        super(TmpBlast, self).makedb()


class BioBlast(TmpBlast):
    def __init__(
        self,
        subjects: typing.Sequence[SeqRecord],
        queries: typing.Sequence[SeqRecord],
        seq_db=None,
        span_origin=True,
        **config
    ):
        """
        If a seq_db is provided, it is assumed subjects and queries already exist in the seq_db.

        :param subjects:
        :param queries:
        :param seq_db:
        :param span_origin:
        :param config:
        """

        self.span_origin = span_origin
        if seq_db is None:
            self.seq_db = SeqRecordDB()
            subjects = self.add_records(subjects)
            queries = self.add_records(queries)
        else:
            self.seq_db = seq_db
        db_name = str(uuid4())
        subject_path = records_to_tmpfile(subjects)
        query_path = records_to_tmpfile(queries)
        super().__init__(
            db_name=db_name, subject_path=subject_path, query_path=query_path, **config
        )

    def add_records(self, records):
        clean_records(records)

        def copy_record(r):
            return deepcopy(r)

        def pseudocircularize(r):
            r2 = r + r
            r2.name = C.PSEUDOCIRCULAR + "__" + r.name
            r2.id = str(uuid4())
            return r2

        circular = [r for r in records if self.seq_db.is_circular(r)]
        linear = [r for r in records if not self.seq_db.is_circular(r)]
        if self.span_origin:
            keys = self.seq_db.add_many_with_transformations(
                circular, pseudocircularize, C.PSEUDOCIRCULAR
            )
        else:
            keys = self.seq_db.add_many_with_transformations(
                circular, copy_record, C.PSEUDOCIRCULAR
            )
        keys += self.seq_db.add_many_with_transformations(
            linear, copy_record, C.COPY_RECORD
        )
        return self.seq_db.get_many(keys)

    # def _merge_results(self, results):
    #
    #     def mergable(d1, d2):
    #         if d1['end'] == d1['length']:
    #             if d1['circular']:
    #                 return True
    #             else:
    #                 return False
    #         elif d1['end'] + 1 == d2['start']:
    #             return True
    #         return False
    #
    #     def alignment_mergable(v1, v2):
    #         return all([mergable(v1[x], v2[x]) for x in ['query', 'subject']])

    # sorted_by_query_ends = sorted(results, key=lambda x: x['query']['end'])
    # grouped_by_qs = groupby_transform(results,
    #                                   keyfunc=lambda x: x['query']['sequence_id'] + "__" + x['subject']['sequence_id'])
    # for _, group in grouped_by_qs:
    #     grouped_by_query_starts = dict(groupby_transform(group, keyfunc=lambda x: x['query']['start']))
    #     for v1 in group:
    #         e = v1['query']['end']
    #         if e == v1['query']['length']:
    #             qs = grouped_by_query_starts.get(0, list())
    #             for v2 in grouped_by_query_starts.get(0, list()):
    #                 if alignment_mergable(v1, v2):
    #

    def parse_results(self, delim=","):
        parsed_results = BlastResultParser.raw_results_to_json(
            self.raw_results, delim=delim
        )

        # # TODO: resolve with sequence dictionary, resolving pseudocircularized constructs
        for v in parsed_results:
            if v:
                for x in ["query", "subject"]:

                    origin_key = self.seq_db.get_origin_key(v[x]["sequence_id"])
                    record = self.seq_db.get(origin_key)
                    is_circular = self.seq_db.is_circular(record)
                    v[x]["circular"] = is_circular
                    v[x]["name"] = record.name
                    v[x]["origin_key"] = origin_key
                    v[x]["origin_record_id"] = record.id
                    v[x]["origin_sequence_length"] = len(record.seq)
                    # v[x]["length"] = len(record.seq)
                v["meta"]["span_origin"] = self.span_origin

        parsed_results = self._unique_results(parsed_results)
        self.results = parsed_results
        return self.results

    def alignments(self):
        alignments = []
        for align in self.results:
            alignments.append(
                BlastResultParser.alignment_to_seqrecord(align, self.seq_db.records)
            )
        return alignments


class JSONBlast(BioBlast):
    """Object that runs blast starting from JSON inputs and outputs"""

    def __init__(self, subject_json, query_json, span_origin=True, seq_db=None, **config):
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
        super().__init__(srecords, qrecords, seq_db=seq_db, span_origin=span_origin)

    @classmethod
    def use_test_data(cls):
        """Create a Blast instance using predefined data located in tests"""
        dir_path = os.path.dirname(os.path.realpath(__file__))
        with open(
            os.path.join(dir_path, "..", "..", "tests/data/test_data/templates.json"),
            "r",
        ) as f:
            subject = json.load(f)
        with open(
            os.path.join(dir_path, "..", "..", "tests/data/test_data/query.json"), "r"
        ) as f:
            query = json.load(f)
        return cls(subject_json=subject, query_json=query)
