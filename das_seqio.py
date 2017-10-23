'''
Project: jdna
File: das_seqio
Author: Justin
Date: 2/6/17

Description: 

'''

import shlex
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
import os
import json
from glob import glob
import re
import coral

def run_cmd(cmd_str, **kwargs):
    """
    Runs a subprocess command with kwargs arguments
    :param cmd_str: formattable command line string (e.g. 'ls {directory_location}')
    :param kwargs: dictionary of arguments for command
    :return: None
    :param cmd_str:
    :param kwargs:
    :return:
    """
    cmd_str = cmd_str.format(**kwargs)
    args = shlex.split(cmd_str)
    print("CMD: {}".format(cmd_str))
    output = subprocess.Popen(args)
    output.wait()
    return output


def format_decorator(f):
    """
    Implies a SeqIO format based on filename suffix
    :param f:
    :return: wrapped function
    """

    def wrapped(*args, **kwargs):
        args = list(args)
        g = re.search('^(.+)\.(\w+)$', args[0])
        prefix = g.group(1)
        suffix = g.group(2)
        formats = {"gb": "genbank", "fsa": "fasta", "fasta": "fasta"}
        formats.update(kwargs)
        kwargs['format'] = formats[suffix]
        return f(*args, **kwargs)

    return wrapped

@format_decorator
def determine_format(filename, format=None):
    return format

def sanitize_filenames(dir):
    replacements = [(' ', '_')]
    for filename in glob(os.path.join(dir, '*')):
        if os.path.isfile(filename):
            print(filename)
            for r in replacements:
                new_filename = re.sub(r[0], r[1], filename)
                os.rename(filename, new_filename)

# TODO: locus_ID gets truncated, fix this...
@format_decorator
def open_sequence(filepath, format=None, **fmt):
    """
    Open sequences from path as Bio.Seq formats
    :param filepath:
    :param format:
    :param fmt:
    :return:
    """
    g = re.search('^(.+)\.(\w+)$', filepath)
    prefix, suffix = g.group(1), g.group(2)
    seqs = []
    with open(filepath, 'rU') as handle:
        s = list(SeqIO.parse(handle, format))
        # if len(s) == 1:
        #     s[0].id = prefix
        seqs += s
    return seqs


@format_decorator
def save_sequence(filename, sequences, format=None, **fmt):
    with open(filename, 'w') as handle:
        SeqIO.write(sequences, handle, format)
    return filename


def concat_seqs(idir, out, savemeta=False):
    """
    Concatenates a directory of sequences into a single fasta file
    :param idir: input directory path
    :param out: output path
    :return: full output path
    """
    sequences = []
    metadata = {}
    filenames = glob(os.path.join(idir, '*.*'))
    for filename in filenames:
        seqs = open_sequence(filename)
        sequences += seqs

        # TODO: this is really hacky, recode this
        # TODO: this requires Coral, do we want another dependency?
        for s in seqs:
            metadata[s.id] = {'filename': filename, 'circular': False, 'seqrecord': s}
        if len(seqs) == 1:
            try:
                c = coral.seqio.read_dna(filename)
                metadata[seqs[0].id]['circular'] = c.circular
                metadata[seqs[0].id]['coral'] = c
            except:
                pass

    with open(out, "w") as handle:
        SeqIO.write(sequences, handle, "fasta")

    if savemeta:
        with open(out.split('.')[0] + '.json', 'w') as handle:
            json.dump(metadata, handle)

    return out, sequences, metadata


def gb_to_fsa(input_path, output_path):
    sequences = open_sequence(input_path)
    with open(output_path, 'w') as handle:
        SeqIO.write(sequences, handle, 'fasta')
    return output_path


def seq_is_circular(path):
    d = coral.seqio.read_dna(path)
    return d.circular


def primers_json_to_fasta(json_file, output_path):
    seqs = []
    pjs = json.load(json_file)
    for primer_json in pjs:
        seq = primer_json['Overhang Sequence'] + primer_json['Anneal Sequence']
        seqs.append(Seq())

