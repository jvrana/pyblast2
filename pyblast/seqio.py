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


# TODO: locus_ID gets truncated, fix this...
@format_decorator
def open_sequence(path, format=None):
    """Open a sequence from a path"""
    seqs = []
    with open(path, 'rU') as handle:
        seqs = list(SeqIO.parse(handle, format))

    for seq in seqs:
        seq.filename = path

    if len(seqs) == 1:
        if dna_at_path_is_circular(path):
            seqs[0].circular = True
    else:
        for s in seqs:
            s.circular = False
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
