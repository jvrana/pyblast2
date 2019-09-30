"""
blast.py

More information: https://www.ncbi.nlm.nih.gov/guide/howto/run-blast-local/
BLAST User Manual: https://www.ncbi.nlm.nih.gov/books/NBK1762/
"""

import json
import os
import shutil
import tempfile
import typing
from uuid import uuid4

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from more_itertools import unique_everseen

from pyblast.blast_bin import BlastWrapper
from pyblast.constants import Constants as C
from pyblast.exceptions import PyBlastException, PyBlastWarning
from pyblast.utils import Span
from pyblast.utils import (
    glob_fasta_to_tmpfile,
    records_to_tmpfile,
    clean_records,
    force_unique_record_ids,
)
from pyblast.utils import run_cmd
from .blast_parser import BlastResultParser
from .json_parser import JSONParser
from .seqdb import SeqRecordDB


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
        config=None,
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

        # name of the database
        self.db_name = None
        self.db_output_directory = None
        self.results_out_path = None
        self.set_paths(db_name, db_output_directory, results_out_path)

        # add executable to the system path
        blast_wrapper = BlastWrapper()
        if not blast_wrapper.is_installed():
            blast_wrapper.ask_to_install()

        self.query_path = query_path
        self.subject_path = subject_path

        # setup config
        if config is None:
            config = {}
        self._config = dict(config)

        # build the configuration
        if output_formatter is None:
            self.outfmt = BlastBase.outfmt
        else:
            self.outfmt = output_formatter

        # list of expected fields in the result
        self.fields = ()

        # the raw results from running blast from the command line
        self.raw_results = None

        # the parsed results as a dictionary
        self.results = None

    def set_paths(self, db_name, db_out_dir, results_out_path):
        self.db_name = db_name
        if db_out_dir:
            self.db_output_directory = os.path.abspath(db_out_dir)
        else:
            self.db_output_directory = db_out_dir
        if results_out_path:
            self.results_out_path = os.path.abspath(results_out_path)
        else:
            self.results_out_path = results_out_path

    @property
    def db(self):
        return os.path.join(self.db_output_directory, self.db_name)

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

    def blastn(self):
        """Alias of 'quick_blastn'"""
        self.quick_blastn()

    def quick_blastn(self):
        """Produce database, run blastn protocol, & parse results """
        self.makedb()
        self._run_blastn()
        self.closedb()
        return self.parse_results()

    def blastn_short(self):
        """Alias of 'quick_blastn_short'"""
        self.quick_blastn_short()

    def quick_blastn_short(self):
        """Produce database, run blastn short protocol, & parse results """
        self.makedb()
        self._run_blastn_short()
        self.closedb()
        return self.parse_results()

    def _run_blastn(self):
        return self._run("blastn")

    def _run_blastn_short(self):
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
    @staticmethod
    def run_cmd(cmd, **kwargs):
        """Wrapper for utils.run_cmd"""
        run_cmd(cmd, **kwargs)

    @staticmethod
    def _fasta_details(path):
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
            title=self.db_name,
            out=self.db,
            **{"in": self.subject_path}
        )
        return self.db

    def closedb(self):
        os.remove(self.results_out_path)
        shutil.rmtree(self.db_output_directory)

    @staticmethod
    def _filter_unique_results(results):
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

    def _filter_results(self, results):
        results = self._filter_unique_results(results)
        return results

    def parse_results(self, delim=","):
        """
        Parses the raw blast result to a JSON

        :param delim: delimiter to parse
        :type delim: str
        :return:
        :rtype:
        """
        results = BlastResultParser.raw_results_to_json(self.raw_results, delim=delim)
        results = self._filter_results(results)
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

    def __init__(self, db_name, subject_path, query_path, config=None):
        """
        TmpBlast: A Blast object that stores the database files in a hidden temporary directory.
        :param db_name: Name for database file structure. This name will be appended to all db files that blast creates.
        :param subject_path: Input directory or file containing a list of subjects to align against the query
        :param query_path: Location of the fasta or genbank file containing the query
        :param config: Additional configurations to run for the blast search (see
        https://www.ncbi.nlm.nih.gov/books/NBK279682/)
        """
        self.subject_records_path = subject_path

        super(TmpBlast, self).__init__(
            db_name=db_name,
            subject_path=None,
            query_path=query_path,
            db_output_directory=None,
            results_out_path=None,
            config=config,
        )

        self.fields = self.fields + (
            "subject_name",
            "subject_filename",
            "subject_circular",
            "query_name",
            "query_filename",
            "query_circular",
        )

    def makedb(self, **kwargs):
        db_output_directory = tempfile.mkdtemp()
        _, out = tempfile.mkstemp(dir=db_output_directory)

        if os.path.isdir(self.subject_records_path):
            self.subject_path = glob_fasta_to_tmpfile(self.subject_records_path)
        else:
            self.subject_path = self.subject_records_path
        self.set_paths(self.db_name, db_output_directory, out)
        super(TmpBlast, self).makedb(**kwargs)

    @classmethod
    def use_test_data(cls):
        """Create a Blast instance using predefined data located in tests"""
        dir_path = os.path.dirname(os.path.realpath(__file__))

        return cls(
            "db",
            os.path.join(dir_path, "..", "..", "tests/data/test_data/db.fsa"),
            os.path.join(dir_path, "..", "..", "tests/data/test_data/query.fsa"),
        )


class BioBlast(TmpBlast):

    DEFAULT_REINDEX = 1
    REINDEX_KEY = "reindex"

    def __init__(
        self,
        subjects: typing.List[SeqRecord],
        queries: typing.List[SeqRecord],
        seq_db=None,
        span_origin=True,
        force_unique_ids=False,
        config=None,
    ):
        """
        Initialize a new BioBlast.

        :param subjects: list of SeqRecords to use as subjects. Subjects get aligned to the queries.
        :type subjects: list
        :param queries: list of SeqRecords to use as queries. Subjects get aligned to the queries.
        :type queries: list
        :param seq_db: optional SeqRecordDB to use
        :type seq_db: SeqRecordDB
        :param span_origin: if False, will treat all sequences as linear (default: True)
        :type span_origin: bool
        :param force_unique_ids: if True, will force any records with the same record_id to be unique.
        :type force_unique_ids: bool
        :param config: additional blast config
        :type config: dict
        """
        if force_unique_ids:
            force_unique_record_ids(subjects)
            force_unique_record_ids(queries)
        self._check_records(subjects, queries)
        self.span_origin = span_origin
        if seq_db is None:
            self.seq_db = SeqRecordDB()
            _, subjects = self.add_records(
                subjects, seq_db=self.seq_db, span_origin=self.span_origin
            )
            _, queries = self.add_records(
                queries, seq_db=self.seq_db, span_origin=self.span_origin
            )
        else:
            self.seq_db = seq_db
        db_name = str(uuid4())
        subject_path = records_to_tmpfile(subjects)
        query_path = records_to_tmpfile(queries)
        super().__init__(
            db_name=db_name,
            subject_path=subject_path,
            query_path=query_path,
            config=config,
        )

        self.parse_options = {self.REINDEX_KEY: self.DEFAULT_REINDEX}

    def _check_records(self, subjects, queries):
        self._check_empty_records(queries, subjects)
        self._check_duplicate_records(queries, subjects)

    @staticmethod
    def _check_duplicate_records(queries, subjects):
        # Check SeqRecords have unique record_ids
        subject_rec_ids = {}
        for s in subjects:
            subject_rec_ids.setdefault(s.id, list()).append(s)
        query_rec_ids = {}
        for s in queries:
            query_rec_ids.setdefault(s.id, list()).append(s)
        duplicate_subjects = {k: v for k, v in subject_rec_ids.items() if len(v) > 1}
        duplicate_queries = {k: v for k, v in query_rec_ids.items() if len(v) > 1}
        suggestion = (
            "This may be due to clipping of the record_id that occurs when saving or loading"
            " certain file formats. Use `{}` or if using `pyblast.utils.load_glob` set the"
            "`force_unique_ids=True`.".format(force_unique_record_ids.__name__)
        )
        if duplicate_subjects:
            raise PyBlastException(
                "One or more subjects have the same record_id:\n{}.\n{}".format(
                    list(duplicate_subjects.keys()), suggestion
                )
            )
        if duplicate_queries:
            raise PyBlastException(
                "One or more queries have the smae record_id:\n{}\n{}".format(
                    list(duplicate_queries.keys()), suggestion
                )
            )

    @staticmethod
    def _check_empty_records(queries, subjects):
        # Check SeqRecords exist
        if not subjects:
            raise ValueError("Subjects is empty.")
        if not queries:
            raise ValueError("Queries is empty.")

    @staticmethod
    def _filter_remove_same_results(results):
        """
        Removes alignments that aligned to themselves.

        :param results:
        :type results:
        :return:
        :rtype:
        """
        new_results = []
        for r in results:
            query = r["query"]
            subj = r["subject"]
            k1 = query["origin_key"]
            k2 = subj["origin_key"]
            qends = (query["start"], query["end"])
            sends = (subj["start"], subj["end"])
            if k1 == k2 and qends == sends:
                continue
            else:
                new_results.append(r)
        return new_results

    @staticmethod
    def _filter_unique_results(results):
        return list(
            unique_everseen(
                results,
                key=lambda x: (
                    x["query"]["sequence_id"],
                    x["query"]["start"],
                    x["query"]["end"],
                    x["query"]["raw_end"],
                    x["subject"]["sequence_id"],
                    x["subject"]["start"],
                    x["subject"]["end"],
                    x["subject"]["raw_end"],
                ),
            )
        )

    def _filter_results(self, results):
        results = self._filter_unique_results(results)
        results = self._filter_remove_same_results(results)
        return results

    @classmethod
    def add_records(
        cls, records: typing.List[SeqRecord], seq_db: SeqRecordDB, span_origin=True
    ) -> typing.Tuple[typing.List[str], typing.List[SeqRecord]]:
        """
        Adds records to the local sequence database (SeqRecordDb). If `self.span_origin=True`,
        then in addition to adding the original SequenceRecord to the database, sequences will
        be pseudocircularized (concatenated with itself) if they are detected to be circular (via the 'topology'
        annotation') so that Blast can properly align sequences over the origin. Only sequences that were linear or
        the pseudociruclarized sequences are used in the alignment. After running blast and obtaining results, span
        indices are automatically corrected to account for the original pseudocircularization.

        :param records: list of SeqRecords annotated with 'topology' being 'linear' or 'circular'
        :type records: list
        :param seq_db: optional SeqRecordDB to use
        :type seq_db: SeqRecordDB
        :param span_origin: whether to search the sequence as a circular plasmid
        :type span_origin: bool
        :return: list of keys to use in the alignment procedure
        :rtype: list
        """
        clean_records(records)

        def copy_record(r):
            c = r[:]
            c.annotations = dict(r.annotations)
            return c

        def pseudocircularize(r):
            r2 = r + r
            r2.name = C.PSEUDOCIRCULAR + "__" + r.name
            r2.id = str(uuid4())
            return r2

        circular = [r for r in records if seq_db.is_circular(r)]
        linear = [r for r in records if not seq_db.is_circular(r)]
        if span_origin:
            keys = seq_db.add_many_with_transformations(
                circular, pseudocircularize, C.PSEUDOCIRCULAR
            )
        else:
            keys = seq_db.add_many_with_transformations(
                circular, copy_record, C.PSEUDOCIRCULAR
            )
        keys += seq_db.add_many_with_transformations(linear, copy_record, C.COPY_RECORD)
        transformed_records = seq_db.get_many(keys)
        return keys, transformed_records

    def parse_result_to_span(
        self, data, inclusive=True, input_index=1, output_index=None
    ):

        if output_index is None:
            output_index = input_index

        s, e, l = data["start"], data["raw_end"], data["origin_sequence_length"]
        if data["strand"] == -1:
            s, e = e, s
        c = data["circular"]
        if inclusive:
            e += 1
        span = Span(s, e, l, cyclic=c, ignore_wrap=False, index=input_index)
        if input_index != output_index:
            span = span.reindex(output_index)
        return span

    @classmethod
    def parse_to_span(cls, v, reindex=0):
        """Convert a JSON blast result to a Span. Optionally reindex"""

        def make_span(data):
            cls.parse_result_to_span(
                data=data,
                inclusive=v["meta"]["inclusive"],
                input_index=v["meta"]["start_index"],
                output_index=reindex,
            )

        return {"query": make_span(v["query"]), "subject": make_span(v["subject"])}

    def _parse__annotate_meta(self, meta):
        """
        Annotate default meta information from the raw JSON result.

        :param meta:
        :return:
        """
        meta["start_index"] = 1
        meta["inclusive"] = True
        meta["span_origin"] = self.span_origin

    def _parse__annotate_record(self, data):
        """
        Annotate record information from the raw JSON result using the sequence database and sequence_id.

        :param data:
        :return:
        """
        origin_key = self.seq_db.get_origin_key(data["sequence_id"])
        record = self.seq_db.get(origin_key)
        is_circular = self.seq_db.is_circular(record)
        data["circular"] = is_circular
        data["name"] = record.name
        data["origin_key"] = origin_key
        data["origin_record_id"] = record.id
        data["origin_sequence_length"] = len(record.seq)

    def _parse__correct_endpoints(
        self, data, inclusive=True, input_index=1, output_index=None
    ):
        """
        Correct endpoints of JSON results.

        :param data: The JSON result. E.g. `{"start": 0, "end": 100, "origin_sequence_length": 200}`
        :param inclusive: If True, the input JSON result is considered to be inclusive at the endpoint.
        :param input_index: The starting index assumed in the input JSON result.
        :param output_index: The index to reindex the result.
        :return:
        """

        if output_index is None:
            output_index = input_index
        data["raw_end"] = data["end"]
        s, e = data["start"], data["end"]
        reverse = data["strand"] != 1
        if reverse:
            s, e = e, s
        if e - s < 0:
            raise ValueError("End cannot be less than start")
        elif data["circular"]:
            span = self.parse_result_to_span(
                data=data,
                inclusive=inclusive,
                input_index=input_index,
                output_index=output_index,
            )
            data["start"] = span.a
            data["end"] = span.b - 1
            data["raw_end"] = span.c - 1
            if reverse:
                data["start"], data["end"], data["raw_end"] = (
                    data["end"],
                    data["start"],
                    data["start"],
                )

    def _parse_bioblast_results(self, v, output_index):
        self._parse__annotate_meta(v["meta"])
        for x in ["query", "subject"]:
            self._parse__annotate_record(v[x])
            self._parse__correct_endpoints(
                v[x],
                inclusive=v["meta"]["inclusive"],
                input_index=v["meta"]["start_index"],
                output_index=output_index,
            )
        v["meta"]["start_index"] = output_index
        return v

    def parse_results(self, delim=",", reindex=None):
        """
        Parses the blast output to a digestable JSON output of the following format. Results
        can be found in `self.results` or returned from this function. Starting and ending
        positions are inclusive and starting index is 1.

        ..code-block::

            [
                {
                  "query": {
                    "start": 1,
                    "end": 4219,
                    "bases": "TTTAGTATATAT...",
                    "strand": 1,
                    "length": 18816,
                    "sequence_id": "pMODKan-HO-pACT1-Z4-"
                  },
                  "subject": {
                    "start": 1,
                    "end": 4219,
                    "bases": "TCGCGCGTTTCGGTGATGACGG...",
                    "strand": 1,
                    "length": 7883,
                    "sequence_id": "pMODKan-HO-pACT1-ZEV4"
                  },
                  "meta": {
                    "query acc.": "pMODKan-HO-pACT1-Z4-",
                    "subject acc.": "pMODKan-HO-pACT1-ZEV4",
                    "score": 4219,
                    "evalue": 0,
                    "bit score": 7792,
                    "alignment length": 4219,
                    "identical": 4219,
                    "gap opens": 0,
                    "gaps": 0,
                    "query length": 18816,
                    "q. start": 1,
                    "q. end": 4219,
                    "subject length": 7883,
                    "s. start": 1,
                    "s. end": 4219,
                    "subject strand": "plus",
                    "query seq": "TTTAGTATATAT...",
                    "subject seq": "TCGCGCGTTTCGGTGATGACGG..."
                  }
                }
            ]

        :param delim:
        :type delim:
        :return:
        :rtype:
        """
        if reindex is None:
            reindex = self.parse_options.get(self.REINDEX_KEY, self.DEFAULT_REINDEX)
        parsed_results = BlastResultParser.raw_results_to_json(
            self.raw_results, delim=delim
        )

        # # TODO: resolve with sequence dictionary, resolving pseudocircularized constructs
        for v in parsed_results:
            if v:
                self._parse_bioblast_results(v, output_index=reindex)
        parsed_results = self._filter_results(parsed_results)
        self.results = parsed_results
        return self.results


class JSONBlast(BioBlast):
    """Object that runs blast starting from JSON inputs and outputs"""

    def __init__(
        self, subject_json, query_json, span_origin=True, seq_db=None, config=None
    ):
        """
        Initialize JSONBlast

        :param subject_json: subject sequences as serialized json. Schema may be found in :class:`pyblast.schema.SequenceSchema`
        :type subject_json: list of dict, mapping, or Object
        :param query_json: query sequence as serialized json
        :type query_json: dict, mapping, or Object
        :param span_origin: default False. If True, circular sequences will be pseudocircularized for alignment over the origin
        :type span_origin: boolean
        :param seq_db: optional SeqRecordDB to use
        :type seq_db: SeqRecordDB
        :param config: optional arguments
        :type config: dict
        """
        srecords = JSONParser.JSON_to_SeqRecords(subject_json)
        qrecords = JSONParser.JSON_to_SeqRecords(query_json)
        self.subject_json = subject_json
        self.query_json = query_json
        super().__init__(
            srecords, qrecords, seq_db=seq_db, span_origin=span_origin, config=config
        )

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
