import logging
import os
import tempfile
from glob import glob

from marshmallow import ValidationError

from pyblast.schema import SequenceSchema
from .dna_bases import rc_dict


def json_to_fasta_tempfile(jsondata, prefix="", id="name"):
    """Writes JSON data to a temporary fasta file"""
    fd, temp_path = tempfile.mkstemp(prefix="query_{}__".format(prefix), suffix=".fasta")
    with open(temp_path, 'w') as out:
        out.write(json_to_fasta_data(jsondata, id=id))
        out.close()
    return temp_path

def json_to_fasta_data(jsondata, id="name"):
    """Converts json to fasta format"""

    seq = parse_sequence_jsons(jsondata)

    def convert(data):
        if "sequence" not in data:
            pass
        return ">{id}\n{sequence}\n".format(id=data[id], sequence=data["sequence"].upper())

    if type(jsondata) is list:
        return '\n'.join([convert(x).strip() for x in seq])
    else:
        return convert(jsondata)


def fasta_to_json(fasta, id="name"):
    data = []
    sequences = fasta.split('>')[1:]
    for seq in sequences:
        tokens = seq.split('\n')
        header = tokens[0]
        seqs = tokens[1:]
        cols = header.split('|')
        data.append({
            "{}".format(id): cols[0],
            "sequence": ''.join(seqs).strip(),
            "circular": False
        })
    return parse_sequence_jsons(data)


def concat_fasta_to_tempfile(dir):
    # concatenate files if subject_path is directory
    seqs = []
    fasta_files = glob(os.path.join(dir, "*.fsa"))
    fasta_files += glob(os.path.join(dir, "*.fasta"))
    for fsa in fasta_files:
        with open(fsa, 'r') as f:
            seqs += fasta_to_json(f.read())
    return json_to_fasta_tempfile(seqs)


def reverse_complement(seq):
    return ''.join([rc_dict[x] for x in seq])


def parse_sequence_jsons(data):
    many = type(data) is list
    schema = SequenceSchema(many=many)
    try:
        result = schema.load(data)
    except ValidationError as err:
        logging.error(err.messages)
        result = err.valid_data
        result = [s[1] for s in enumerate(result) if s[0] not in err.messages.keys()]
    return result
