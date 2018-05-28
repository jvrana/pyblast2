import os
import tempfile
from glob import glob
from pyblast.schema import SequenceSchema


def json_to_fasta_tempfile(jsondata, prefix="", id="name"):
    """Writes JSON data to a temporary fasta file"""
    fd, temp_path = tempfile.mkstemp(prefix="query_{}__".format(prefix), suffix=".fasta")
    with open(temp_path, 'w') as out:
        out.write(json_to_fasta_data(jsondata, id=id))
        out.close()
    return temp_path


def json_to_fasta_data(jsondata, id="name"):
    """Converts json to fasta format"""

    schema = SequenceSchema()
    seq = schema.load(jsondata, many=type(jsondata) is list)

    def convert(data):
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
    schema = SequenceSchema(many=True)
    return schema.load(data)

def concat_fasta_to_tempfile(dir):
    # concatenate files if subject_path is directory
    seqs = []
    fasta_files = glob(os.path.join(dir, "*.fsa"))
    fasta_files += glob(os.path.join(dir, "*.fasta"))
    for fsa in fasta_files:
        with open(fsa, 'r') as f:
            seqs += fasta_to_json(f.read())
    return json_to_fasta_tempfile(seqs)
