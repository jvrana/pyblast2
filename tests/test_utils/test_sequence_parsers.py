from pyblast.utils import fasta_to_json, json_to_fasta_data, json_to_fasta_tempfile, concat_fasta_to_tempfile
import os
import json

def test_fast_to_json(here):
    data = """
>16079__W17-pJZC-gRNATarSeqSwap-R .
A
GGG
TTT


>1867__BFL1_pETCON_R .



GCTTTTGTTCGGATCCGCCCCCCTCGAGACCTGATTTCGGTTCAAATTTTTTC



"""
    seqs = fasta_to_json(data)
    assert len(seqs) == 2

    assert seqs[0] == {
            "name": "16079__W17-pJZC-gRNATarSeqSwap-R .",
            "sequence": "AGGGTTT",
            "circular": False,
            "description": "",
            "features": [],
            "size": 7,
            "id": seqs[0]["id"],
            "notes": {}
        }

    assert seqs[1] == {
            "name": "1867__BFL1_pETCON_R .",
            "sequence": "GCTTTTGTTCGGATCCGCCCCCCTCGAGACCTGATTTCGGTTCAAATTTTTTC",
            "circular": False,
            "description": "",
            "features": [],
            "size": 53,
            "id": seqs[1]["id"],
            "notes": {}
        }


def test_json_to_fasta_data():
    here = os.path.dirname(os.path.abspath(__file__))

    with open(os.path.join(here, "dna.json")) as f:
        data = json.load(f)
        fasta = json_to_fasta_data(data)
    expected = """>pMODKan-HO-pACT1-Z4-
AGTCAGATAC
"""
    assert fasta == expected


def test_json_to_fasta_tempfile():
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "dna.json")) as f:
        data = json.load(f)
        path = json_to_fasta_tempfile(data)
    expected = """>pMODKan-HO-pACT1-Z4-
AGTCAGATAC
"""
    with open(path, 'r') as f:
        assert f.read() == expected

def test_concat_fasta_to_tempfile():
    here = os.path.dirname(os.path.abspath(__file__))
    concat = concat_fasta_to_tempfile(os.path.join(here, "files"))
    expected = """>myseq
AGCAGGAG
>myseq2
AGGCGGGAGAGG
>myseq3
AGGAGGGGGGG
>myseq4
TTTTTTTT"""
    with open(concat, 'r') as f:
        assert expected == f.read()
