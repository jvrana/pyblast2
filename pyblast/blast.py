"""
blast.py

More information: https://www.ncbi.nlm.nih.gov/guide/howto/run-blast-local/
BLAST User Manual: https://www.ncbi.nlm.nih.gov/books/NBK1762/
"""

import json
import os
import re
import shutil
import tempfile
from uuid import uuid4
from blast_bin import install_blast
from pyblast.utils import run_cmd
from pyblast.utils import json_to_fasta_tempfile, concat_fasta_to_tempfile
from pyblast.results import AlignmentResults
from pyblast.utils.seq_parser import dump_sequence_jsons, load_sequence_jsons
from pyblast.utils import reverse_complement
from pyblast.exceptions import PyBlastException
from marshmallow import ValidationError
from copy import copy


class Blast(object):
    """
    A Blast initializer for running blast searches against subjects contained in a directory. Sequences
    are kept on the local machine.

    Note: This is not intended to be the endpoint. You should only use temporary files to run the blast search.
    """

    outfmt = [
        "7", "qacc", "sacc",
        "score", "evalue", "bitscore", "length", "nident",
        "gapopen", "gaps", "qlen", "qstart", "qend",
        "slen", "sstart", "send", "sstrand",
        "qseq", "sseq"
    ]

    blast_config = {
        "outfmt": "\"{0}\"".format(' '.join(outfmt))
    }

    def __init__(self, db_name, subject_path, query_path, db_output_directory, results_out_path,
                 output_formatter=None, **config):
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
        Blast.add_to_sys_paths()
        self.executable_path = Blast.get_executable()
        print("BLAST executable: {}".format(self.executable_path))

        # name of the database
        self.name = db_name

        self.query_path = query_path
        self.subject_path = subject_path

        # build the configuration
        if output_formatter is None:
            self.outfmt = Blast.outfmt
        self.config = {"outfmt": "\"{0}\"".format(' '.join(self.outfmt))}
        if config is not None:
            self.config.update(config)

        # path to the directory holding the blast database
        self.db_output_directory = os.path.abspath(db_output_directory)

        # path to the blast database
        self.db = os.path.join(db_output_directory, db_name)

        #list of expected fields in the result
        self.fields = ()

        # the raw results from running blast from the command line
        self.raw_results = None

        # the parsed results as a dictionary
        self.results = None  #

        # the path to the saved results
        self.results_out_path = os.path.abspath(results_out_path)


    @staticmethod
    def add_to_sys_paths():
        """Add the path located in blast_bin/_paths.txt to the environment in an attempt to run blast"""
        if not Blast.has_executable():
            install_blast.add_paths_to_environment()
        if not Blast.has_executable():
            error_message = "BLAST executables not found in path. Be sure BLAST is correctly installed."
            help_message = "Please run 'install_pyblast <youremail> <yourplatform>' in your terminal." \
                           "Run 'install_pyblast -h' for help."
            raise PyBlastException(error_message + "\n" + "*" * 50 + "\n" + help_message)

    @staticmethod
    def has_executable():
        """Whether blast is installed and executable"""
        return install_blast.has_executable()

    @staticmethod
    def get_executable():
        """Find the executable path for 'makeblastdb'"""
        return os.path.dirname(shutil.which("makeblastdb"))

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
        self.parse_results()

    def quick_blastn_short(self):
        """Produce database, run blastn short protocol, & parse results """
        self.makedb()
        self.blastn_short()
        self.parse_results()

    def blastn(self):
        self._run("blastn")

    def blastn_short(self):
        self._run("blastn", task="blastn-short")

    def _run(self, cmd, **config_opts):
        """Run the blastn using the current configuration"""
        config = self.create_config()
        config.update(config_opts)
        self.run_cmd(cmd, **config)
        with open(self.results_out_path, 'rU') as handle:
            self.raw_results = handle.read()

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

    def makedb(self):
        """Creates a blastdb from sequences grabbed from the input directory"""
        self.validate_files()
        self.run_cmd("makeblastdb", dbtype="nucl", title=self.name, out=self.db, **{"in": self.subject_path})
        return self.db

    def parse_results(self, save_as_json=True, delim=','):
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

        self.results = AlignmentResults.parse_results(self.raw_results, delim=delim)
        if save_as_json:
            path = os.path.join(self.results_out_path + ".json")
            self.results.dump_to_json(path)
        return self.results

    def __str__(self):
        return "{}".format(self.create_config())


class Aligner(Blast):
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
            results_out_path=out, **config)

        self.fields = self.fields + (
            'subject_name',
            'subject_filename',
            'subject_circular',
            'query_name',
            'query_filename',
            'query_circular')

    @classmethod
    def use_test_data(cls):
        """Create a Blast instance using predefined data located in tests"""
        dir_path = os.path.dirname(os.path.realpath(__file__))

        return cls('db',
                   os.path.join(dir_path, '..', 'tests/data/test_data/db.fsa'),
                   os.path.join(dir_path, '..', 'tests/data/test_data/query.fsa'))

    def makedb(self):
        super(Aligner, self).makedb()


class JSONBlast(Aligner):
    """Object that runs blast starting from JSON inputs and outputs"""

    def __init__(self, subject_json, query_json, preloaded=False, span_origin=False, **config):
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
        dbname = str(uuid4())

        # force json to sequence schema
        if not preloaded:
            try:
                query_json = load_sequence_jsons(query_json)
            except ValidationError as e:
                raise PyBlastException("Validation error while deserializing query sequence data.\n"
                                       " Received a {datatype} and attempted to load using schema.\n"
                                       " Are you sure you shouldn't add 'preloaded=True'?\n"
                                       " ValidationError: {error}".format(datatype=type(query_json), error=e.messages))
            try:
                subject_json = load_sequence_jsons(subject_json)
            except ValidationError as e:
                raise PyBlastException("Validation error while deserializing subject sequence data.\n"
                                       " Received a {datatype} and attempted to load using schema.\n"
                                       " Are you sure you shouldn't add 'preloaded=True'?\n"
                                       " ValidationError: {error}".format(datatype=type(query_json), error=e.messages))

        try:
            self.query = dump_sequence_jsons(query_json)
        except ValidationError as e:
            raise PyBlastException("There was a parsing error while parsing the query sequence.\n{}\n{}".format(e.messages, query_json))
        try:
            self.subjects = dump_sequence_jsons(subject_json)
        except ValidationError as e:
            raise PyBlastException("There was a parsing error while parsing the subject sequences.\n{}\n{}".format(e.messages, subject_json))

        _subjects = self.subjects
        _query = self.query
        self.__span_origin = span_origin
        if self.__span_origin:
            _subjects = [self.pseudocircularize(s) for s in self.subjects]
            _query = self.pseudocircularize(self.query)

        # create temporary files
        subject_path = json_to_fasta_tempfile(_subjects, id="id")
        query_path = json_to_fasta_tempfile(_query, id="id")

        # seq_db
        keys = [x['id'] for x in self.subjects]
        keys.append(self.query['id'])
        seqs = self.subjects + [self.query]
        self.seq_dict = dict(zip(keys, seqs))

        super(JSONBlast, self).__init__(db_name=dbname,
                                        subject_path=subject_path,
                                        query_path=query_path,
                                        **config)

    @classmethod
    def use_test_data(cls):
        """Create a Blast instance using predefined data located in tests"""
        dir_path = os.path.dirname(os.path.realpath(__file__))
        with open(os.path.join(dir_path, '..', 'tests/data/test_data/templates.json'), "r") as f:
            subject = json.load(f)
        with open(os.path.join(dir_path, '..', 'tests/data/test_data/query.json'), "r") as f:
            query = json.load(f)
        return cls(subject_json=subject,
                   query_json=query)

    def pseudocircularize(self, seq):
        if not seq['circular']:
            return seq
        seq_copy = copy(seq)
        seq_copy['sequence'] = seq['sequence']*2
        seq_copy['size'] = seq['size']*2
        return seq_copy

    def parse_results(self, save_as_json=True, delim=','):
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

        # TODO: is replacing the raw text the right thing to do?
        raw = re.sub('Query_1', self.query['id'], self.raw_results)
        results = AlignmentResults.parse_results(raw, delim, context={"db": self.seq_dict})

        if self.__span_origin:
            for alignment in results.alignments:
                alignment['meta']['span_origin'] = True

        if save_as_json:
            path = os.path.join(self.results_out_path + ".json")
            results.dump_to_json(path)

        self.results = results
        return self.results

    def find_perfect_matches(self, min_match, filter=None):
        raise DeprecationWarning("This method is depreciated. Please use blastn with --task blastn-short")
        """Finding perfect matches using python (i.e. not BLAST)"""

        # get the query sequence as a string

        query_seq = self.query["sequence"].upper()
        for subj in self.subjects:
            subj_seq = subj["sequence"].upper()
            q_seq = query_seq
            q_rc_seq = reverse_complement(query_seq)

            for match in re.finditer(subj_seq, q_seq):
                print(subj['name'])
                print(match)
                pass
            for rc_match in re.finditer(subj_seq, q_rc_seq):
                print(subj['name'])
                print(rc_match)
                pass
