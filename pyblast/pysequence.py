from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from uuid import uuid4
import re
import os
from pyblast import utils


class PySeqDB(object):
    """A dictionary like class for storing sequences and prepping them for BLAST"""

    def __init__(self):
        self._db = {}

    @property
    def db(self):
        return dict(self._db)

    @property
    def sequences(self):
        return list(self.db.values())

    @property
    def ids(self):
        return list(self.db.keys())

    def add(self, path):
        pyseqs = PySequence.parse(path)
        for seq in pyseqs:
            seq.id = str(uuid4())
            self._db[seq.id] = seq
        return pyseqs

    def add_from_directory(self, path):
        sequences = []
        for filename in os.listdir(path):
            seq_path = os.path.join(path, filename)
            self.add(seq_path)
        return sequences

    def concatenate_and_save(self, out, fmt="fasta"):
        """Concatenate sequences and save as a single file"""
        with open(out, 'w') as handle:
            SeqIO.write(self.sequences, handle, fmt)
        return out


class PySequence(SeqRecord):
    """Extension of the :class:`SeqRecord` model. Has additional features such as:

        * filename
        * circular
    """

    def __init__(self, seq, id, name, description, dbxrefs, features, annotations, letter_annotations, filename=None, circular=None):
        super(PySequence, self).__init__(seq, id, name, description, dbxrefs, features, annotations, letter_annotations)
        self.filename = filename
        self.alias = id
        if circular is None:
            circular = False
        self.circular = circular

    @classmethod
    def from_seqrecord(cls, seq):
        """Create a PySequence from a SeqRecord object"""
        filename = None
        if hasattr(seq, 'filename'):
            filename = seq.filename

        circular = None
        if hasattr(seq, 'circular'):
            circular = seq.circular

        return cls(seq.seq, seq.id, seq.name, seq.description, seq.dbxrefs, seq.features, seq.annotations,
            seq.letter_annotations, filename, circular)

    @staticmethod
    def which_format(path):
        """Determine which format to parse the sequence from the path"""
        extension_options = {
            'genbank': ['.gb', 'ape'],
            'fasta': ['.fasta', '.fa', '.fsa', '.seq']
        }

        dir, filename = os.path.split(path)
        basename, ext = os.path.splitext(filename)
        for format, exts in extension_options.items():
            if ext in exts:
                return format

    @staticmethod
    def dna_at_path_is_circular(path):
        """
        Whether a genbank file at "path" is circular.

        :param path:
        :type path:
        :return:
        :rtype:
        """
        with open(path, 'rU') as myfile:
            first_line = myfile.readlines()[0]
            match = re.search("circular", first_line, re.IGNORECASE)
            return match is not None

    @classmethod
    def open(cls, path):
        return cls.parse(path)

    @classmethod
    def parse(cls, path):
        """Parse a sequence file to a PySequence"""
        with open(path, 'rU') as handle:
            seqs = list(SeqIO.parse(handle, cls.which_format(path)))

            circular = False

            if len(seqs) == 1 and cls.dna_at_path_is_circular(path):
                circular = True

            for seq in seqs:
                seq.circular = circular
                seq.filename = path

            return [PySequence.from_seqrecord(seq) for seq in seqs]

    @staticmethod
    def save_sequences(path, sequences):
        """
        Save a list of sequences to a file.

        :param path: output path
        :type path: basestring
        :param sequences: list of sequences
        :type sequences: list
        :return:
        :rtype:
        """
        with open(path, 'w') as handle:
            SeqIO.write(sequences, handle, PySequence.which_format(path))
        return path

    @staticmethod
    def sanitize_filename(filename, replacements=None):
        """
        Santitized filename according to replacements list.

        :param filename:
        :type filename:
        :param replacements:
        :type replacements:
        :return:
        :rtype:
        """
        if replacements is None:
            replacements = [(' ', '_')]
        new_filename = ""
        for repl in replacements:
            new_filename = re.sub(repl[0], repl[1], filename)
        return new_filename

    @staticmethod
    def sanitize_filenames(dir, replacements=None, odir=None):
        """
        Sanitizes all filenames according to replacements list

        :param dir: input directory
        :type dir: str
        :param replacements: list of tuples indicating replacements
        :type replacements: list
        :param odir: output directory to save new files (optional)
        :type odir: str
        :return: None
        :rtype: None
        """
        if odir is None:
            odir = dir
        for filename in os.listdir(dir):
            if os.path.isfile(os.path.join(dir, filename)):
                newfilename = PySequence.sanitize_filename(filename, replacements=replacements)
                os.rename(os.path.join(dir, filename), os.path.join(odir, newfilename))

    def save(self, path):
        """
        Save PySequence to a file

        :param path: output path
        :type path: basestring
        :return:
        :rtype:
        """
        return PySequence.save_sequences(path, [self])
