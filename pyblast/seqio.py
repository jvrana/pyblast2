"""Sequence parser"""

import json
import os
import re
import uuid

from Bio import SeqIO


# import csv
# from Bio.Alphabet.IUPAC import ambiguous_dna
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
# from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
# from Bio.SeqFeature import CompoundLocation

def split_path(path):
    """Split filename into dir, filename, basename, and ext"""
    dir, filename = os.path.split(path)
    basename, ext = os.path.splitext(filename)
    return [dir, filename, basename, ext]


def file_format(ext, **kwargs):
    """Return format from an extension"""
    extension_options = {
        "genbank": ['.gb', '.ape'],
        "fasta": ['.fasta', '.fa', '.fsa', '.seq']
    }
    extension_options.update(kwargs)
    for format, exts in extension_options.items():
        if ext in exts:
            return format
    raise ValueError("File format \"{}\" not recognized".format(ext))


def format_decorator(f):
    """
    Implies a SeqIO format based on filename suffix
    :param f:
    :return: wrapped function
    """

    def wrapped(*args, **kwargs):
        args = list(args)
        ext = split_path(args[0])[-1]
        kwargs['format'] = file_format(ext, **kwargs)
        return f(*args, **kwargs)

    return wrapped


# TODO: locus_ID gets truncated, fix this...
@format_decorator
def open_sequence(path, format=None):
    """Open a sequence from a path"""
    seqs = []
    with open(path, 'rU') as handle:
        s = list(SeqIO.parse(handle, format))
        seqs += s
    return seqs


@format_decorator
def save_sequence(path, sequences, format=None):
    """
    Save a sequence

    :param path:
    :type path:
    :param sequences:
    :type sequences:
    :param format:
    :type format:
    :param fmt:
    :type fmt:
    :return:
    :rtype:
    """
    with open(path, 'w') as handle:
        SeqIO.write(sequences, handle, format)
    return path


@format_decorator
def determine_format(path, format=None):
    """
    Determine the format of a sequence at a path

    :param path:
    :type path:
    :param format:
    :type format:
    :return:
    :rtype:
    """
    return format


def dna_at_path_is_circular(path):
    """
    Whether a genbank file at "path" is circular.

    :param path:
    :type path:
    :return:
    :rtype:
    """
    with open(path) as myfile:
        first_line = myfile.readlines()[0]
        match = re.search("circular", first_line, re.IGNORECASE)
        return match is not None


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
            newfilename = sanitize_filename(filename, replacements=replacements)
            os.rename(os.path.join(dir, filename), os.path.join(odir, newfilename))


def concat_seqs(idir, out, savemeta=False):
    """
    Concatenates a directory of sequences into a single fasta file

    :param idir: input directory
    :type idir: str
    :param out: output path
    :type out: str
    :param savemeta: save metadata associated with each file in separate meta json file? (e.g {filename: ...,
    circular: ...})
    :type savemeta: bool
    :return: ( output path (str), sequences (list of SeqIO), metadata (dict) )
    :rtype: tuple
    """
    sequences = []
    metadata = {}
    for filename in os.listdir(idir):
        seq_path = os.path.join(idir, filename)
        seqs = open_sequence(seq_path)
        sequences += seqs

        for seq in seqs:
            seq.id = str(uuid.uuid4())
            seq.filename = filename
            seq.circular = dna_at_path_is_circular(seq_path)
            metadata[seq.id] = {'circular': seq.circular, 'filename': seq.filename}

    with open(out, "w") as handle:
        SeqIO.write(sequences, handle, "fasta")

    if savemeta:
        with open(out.split('.')[0] + '.json', 'w') as handle:
            json.dump(metadata, handle)

    return out, sequences, metadata

    # for filename in glob(os.path.join(dir, '*')):
    #     if os.path.isfile(filename):
    #         print(filename)
    #         for r in replacements:
    #             new_filename = re.sub(r[0], r[1], filename)
    #             os.rename(filename, new_filename)
