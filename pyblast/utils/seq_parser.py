"""seq_parser.py"""

import logging
import os
import tempfile
from glob import glob

from marshmallow import ValidationError

from pyblast.schema import SequenceSchema
from pyblast.utils.dna_bases import rc_dict


def json_to_fasta_tempfile(jsondata, prefix="", id="name"):
    """Writes serialized sequence model data to a temporary fasta file"""
    fd, temp_path = tempfile.mkstemp(
        prefix="query_{}__".format(prefix), suffix=".fasta"
    )
    with open(temp_path, "w") as out:
        out.write(json_to_fasta_data(jsondata, id=id))
        out.close()
    return temp_path


def json_to_fasta_data(jsondata, id="name"):
    """Converts a serialized sequence model to fasta format"""
    loaded = load_sequence_jsons(jsondata)
    seq = dump_sequence_jsons(loaded)

    def convert(data):
        if "sequence" not in data:
            pass
        return ">{id}\n{bases}\n".format(id=data[id], bases=data["bases"].upper())

    if type(jsondata) is list:
        return "\n".join([convert(x).strip() for x in seq])
    else:
        return convert(jsondata)


def fasta_to_json(fasta, id="name"):
    """
    Dumps a fasta file to serialize sequence

    :param fasta:
    :type fasta:
    :param id:
    :type id:
    :return:
    :rtype:
    """
    data = []
    sequences = fasta.split(">")[1:]
    for seq in sequences:
        tokens = seq.split("\n")
        header = tokens[0]
        seqs = tokens[1:]
        cols = header.split("|")
        data.append(
            {
                "{}".format(id): cols[0],
                "bases": "".join(seqs).strip(),
                "circular": False,
            }
        )
    return dump_sequence_jsons(data)


def concat_fasta_to_tempfile(dir):
    # concatenate files if subject_path is directory
    seqs = []
    fasta_files = glob(os.path.join(dir, "*.fsa"))
    fasta_files += glob(os.path.join(dir, "*.fasta"))
    for fsa in fasta_files:
        with open(fsa, "r") as f:
            seqs += fasta_to_json(f.read())
    return json_to_fasta_tempfile(seqs)


def complement(seq):
    """Complement a DNA sequence"""
    return "".join([rc_dict[x] for x in seq])


def reverse_complement(seq):
    """Reverse complement a sequence"""
    return complement(seq)[::-1]


def load_sequence_jsons(preloaded_data):
    many = issubclass(preloaded_data.__class__, list)
    schema = SequenceSchema(many=many)
    return schema.load(preloaded_data)


def dump_sequence_jsons(postloaded_data):
    """Forces json data into the SequenceSchema """
    many = issubclass(postloaded_data.__class__, list)
    schema = SequenceSchema(many=many)
    try:
        dumped = schema.dump(postloaded_data)
    except ValidationError as err:
        logging.error(err.messages)
        if many:
            dumped = err.valid_data
            dumped = [
                s[1] for s in enumerate(dumped) if s[0] not in err.messages.keys()
            ]
        else:
            raise err
    return dumped
