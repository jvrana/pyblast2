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
from .seqio import split_path, concat_seqs
from .utils import run_cmd, str_to_f_to_i


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

    def __init__(self, db_name, subj_in_dir, query_path, db_output_directory, results_out_path,
                 output_formatter=None, **config):
        """
        A Blast initializer for running blast searches.

        :param db_name: Name for database file structure. This name will be appended to all db
        files that blast creates.
        :param subj_in_dir: Input directory containing a list of subjects to align against the query
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

        # path to the input directory to build the blast database
        self.path_to_input_dir = os.path.abspath(subj_in_dir)

        # path to the query sequence
        self.path_to_query = os.path.abspath(query_path)

        # build the configuration
        if output_formatter is None:
            self.outfmt = Blast.outfmt
        self.config = {"outfmt": "\"{0}\"".format(' '.join(self.outfmt))}
        if config is not None:
            self.config.update(config)

        # path to the directory holding the blast database
        self.path_to_output_dir = os.path.abspath(db_output_directory)

        # path to the blast database
        self.db = os.path.join(db_output_directory, db_name)

        # path to the concatenated sequence file
        self.path_to_input_seq_file = None

        # dictionary of metadata collected from the blast database
        # this includes original filename and topology
        self.db_input_metadata = None

        # the raw results from running blast from the command line
        self.raw_results = None

        # the parsed results as a dictionary
        self.results = None  #

        # list of input sequences that made the blast database
        self.input_sequences = []

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
            raise Exception(error_message + "\n" + "*" * 50 + "\n" + help_message)

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
        """

        def _is_file(myfile):
            return os.path.isfile(os.path.abspath(myfile))

        def _is_dir(mydir):
            return os.path.isdir(os.path.abspath(mydir))

        outdir = split_path(self.results_out_path)[0]
        errors = []
        for file_ in [self.path_to_query]:
            if not _is_file(file_):
                errors.append("File not found: {}".format(file_))
        for dir_ in [outdir, self.path_to_output_dir, self.path_to_input_dir]:
            if not _is_dir(dir_):
                errors.append("Directory not found {}".format(dir_))
        if len(errors) > 0:
            raise ValueError("\n".join(errors))

    def create_config(self):
        """Create a configuration dictionary"""
        config_dict = {
            "db": self.db,
            "out": self.results_out_path,
            "query": self.path_to_query,
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

    def concat_templates(self):
        """
        Gathers all of the sequences in the input dir and concatenates them into a single fasta file

        :return: [concatenated fasta file, sequences used to make fassta file, metadata for the sequences
        :rtype: list
        """
        out = self.db + '.fsa'
        fasta, seqs, metadata = concat_seqs(self.path_to_input_dir, out, savemeta=True)
        self.db_input_metadata = metadata

        self.input_sequences = seqs
        return out, seqs, metadata

    def get_is_circular(self, seqid):
        """Whether the sequence given the sequence id has circular topology"""
        return self.db_input_metadata[seqid]['circular']

    def get_filename(self, seqid):
        """Get the filename from the sequence id"""
        return self.db_input_metadata[seqid]['filename']

    def makedb(self):
        """Creates a blastdb from sequences grabbed from the input directory"""
        out, seqs, metadata = self.concat_templates()
        return self.fasta_to_db(out)

    def fasta_to_db(self, fasta):
        """Create a blastdb from a concatenated fasta file"""
        self.run_cmd("makeblastdb", dbtype="nucl", title=self.name, out=self.db, **{"in": fasta})
        self.path_to_input_seq_file = fasta
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

        def cleanup_fields(match_fields, replacements=None):
            ''' Cleanup field names using replacements '''
            if replacements is None:
                replacements = {
                    ('.', ''),
                    (' ', '_'),
                    ('%', 'perc'),
                }
            match_fields = [x.strip() for x in match_fields]
            for i, f in enumerate(match_fields):
                for r in replacements:
                    match_fields[i] = match_fields[i].replace(r[0], r[1])
            return match_fields

        def extract_metadata(r, delim=','):
            g = re.search(
                '#\s*(?P<blast_ver>.+)\n' +
                '# Query:\s*(?P<query>.*)\n' +
                '# Database:\s*(?P<database>.+)\n' +
                '# Fields:\s*(?P<fields>.+)',
                r)
            metadata = g.groupdict()
            # clean up fields
            metadata['fields'] = re.split('\s*{}\s*'.format(delim), metadata['fields'])
            metadata['fields'] = cleanup_fields(metadata['fields'])
            return metadata

        def extract_raw_matches(r):
            return re.findall('\n([^#].*)', r)

        def validate_matches(raw_matches, fields):
            match_dicts = []
            for m in raw_matches:
                values = [str_to_f_to_i(v) for v in m.split('\t')]
                match_dicts.append(dict(list(zip(fields, values))))
            return match_dicts

        # print(self.results)
        results = self.raw_results
        if results.strip() == '':
            return {}
        meta = extract_metadata(results, delim)
        fields = meta['fields']
        raw_matches = extract_raw_matches(results)
        match_dicts = validate_matches(raw_matches, fields)

        if save_as_json:
            dir, filename, basename, ext = split_path(self.results_out_path)
            f = os.path.join(dir, basename + ".json")
            self.dump_to_json(f)
        self.results = match_dicts
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

    def __init__(self, db_name, subj_in_dir, query_path, **config):
        """
        Aligner: A Blast object that stores the database files in a hidden temporary directory.
        :param db_name: Name for database file structure. This name will be appended to all db files that blast creates.
        :param subj_in_dir: Input directory containing a list of subjects to align against the query
        :param query_path: Location of the fasta or genbank file containing the query
        :param config: Additional configurations to run for the blast search (see
        https://www.ncbi.nlm.nih.gov/books/NBK279682/)
        """
        db_output_directory = tempfile.mkdtemp()
        out = tempfile.mktemp(dir=db_output_directory)
        super(Aligner, self).__init__(db_name, subj_in_dir, query_path, db_output_directory, out, **config)

    @classmethod
    def use_test_data(cls):
        """Create a Blast instance using predefined data located in tests"""
        dir_path = os.path.dirname(os.path.realpath(__file__))

        return cls('db',
                   os.path.join(dir_path, '..', 'tests/data/test_data/templates'),
                   os.path.join(dir_path, '..', 'tests/data/test_data/designs/pmodkan-ho-pact1-z4-er-vpr.gb'))
