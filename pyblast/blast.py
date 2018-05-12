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

from blast_bin import install_blast
from .pysequence import PySequence, PySeqDB
from .utils import run_cmd
from .parser import parse_results

from pyblast.models import QuerySchema, SubjectSchema, AlignmentMetaSchema, AlignmentSchema

class PyBlastException(Exception):
    """A generic exception for pyBlast"""


class Blast(object):
    """
    A Blast initializer for running blast searches against subjects contained in a directory.
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

        # validate the files to make sure they exist
        self.validate_files()


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

        seq = PySequence.open(self.query_path)
        if len(seq) == 0:
            raise PyBlastException("Query path \"{}\" has no sequences".format(self.query_path))
        elif len(seq) > 1:
            raise PyBlastException("Query path \"{}\" has more than one sequence.".format(self.query_path))

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
        """Make a db, run blastn, parse results"""
        self.makedb()
        self.blastn()
        self.parse_results()

    def blastn(self):
        """Run the blastn using the current configuration"""
        self.run_cmd("blastn", **self.create_config())
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

        self.results = parse_results(self.raw_results, delim=delim)
        if save_as_json:
            path = os.path.join(self.results_out_path + ".json")
            self.dump_to_json(path)
        return self.results

    def dump_to_json(self, path):
        with open(path, 'w') as out:
            json.dump(self.results, out)

    def __str__(self):
        return "{}".format(self.create_config())


class Aligner(Blast):
    """
    A Blast object that stores the database files in a hidden temporary directory. Use entry points
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
        out = tempfile.mktemp(dir=db_output_directory)

        # seq_db
        seq_db = PySeqDB()
        seq_db.add_from_directory(subject_path)
        subject_path = os.path.join(db_output_directory, f"{db_name}.fsa")
        seq_db.concatenate_and_save(subject_path)


        self.seq_db = seq_db
        super(Aligner, self).__init__(db_name, subject_path, query_path, db_output_directory, out, **config)
        # add query_path
        self._query = self.seq_db.add(query_path)[0]

        self.fields = self.fields + (
            'subject_name',
            'subject_filename',
            'subject_circular',
            'query_name',
            'query_filename',
            'query_circular')

    def find_perfect_matches(self):
        """Finding perfect matches using python (i.e. not BLAST)"""

        # get the query sequence as a string
        query = PySequence.open(self.query_path)[0]
        for subj in self.seq_db.sequences:
            subj_seq = str(subj.seq).upper()
            q_seq = str(query.seq).upper()
            q_rc_seq = str(query.reverse_complement().seq).upper()

            for match in re.finditer(subj_seq, q_seq):
                print(match)
            for rc_match in re.finditer(subj_seq, q_rc_seq):
                print(rc_match)


    @classmethod
    def use_test_data(cls):
        """Create a Blast instance using predefined data located in tests"""
        dir_path = os.path.dirname(os.path.realpath(__file__))

        return cls('db',
                   os.path.join(dir_path, '..', 'tests/data/test_data/templates'),
                   os.path.join(dir_path, '..', 'tests/data/test_data/designs/pmodkan-ho-pact1-z4-er-vpr.gb'))

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
        raw = re.sub('Query_1', self._query.id, self.raw_results)
        self.results = parse_results(raw, delim, context={"db": self.seq_db.db})
        if save_as_json:
            path = os.path.join(self.results_out_path + ".json")
            self.dump_to_json(path)
        return self.results

