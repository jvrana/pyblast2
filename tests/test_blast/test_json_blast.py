"""testing expected bases for query and subject for JSON Blast."""
import json
import os

import pytest

from pyblast.blast import JSONBlast
from pyblast.exceptions import PyBlastException
from pyblast.utils import reverse_complement


def test_basic_run():
    junk1 = "atgctatgctgatgctgctgtgctgatgctgatgtgtattgctgtatcgcgcgagttagc"
    junk2 = "g" * 30
    frag = "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta"

    query = {"id": "myquery", "bases": frag, "name": "myseq", "circular": False}

    subject = {
        "id": "mysubject",
        "bases": junk1 + frag + junk2,
        "name": "myseq",
        "circular": False,
    }

    # print(type(query))
    # print(type(subject))
    blaster = JSONBlast([subject], query)
    blaster.quick_blastn()
    assert len(blaster.results)


class TestJSONBlast:
    @pytest.fixture
    def aligner(self, here):
        with open(os.path.join(here, "data/test_data/templates.json"), "r") as f:
            subjects = json.load(f)
        with open(os.path.join(here, "data/test_data/query.json"), "r") as f:
            query = json.load(f)
        j = JSONBlast(query, subjects)
        j.quick_blastn()
        return j

    def test_example(self):
        j = JSONBlast.use_test_data()
        j.quick_blastn()
        pass

    def test_perfect_matches(self, here):
        with open(os.path.join(here, "data/test_data/primer.json"), "r") as f:
            subjects = json.load(f)
        with open(os.path.join(here, "data/test_data/query.json"), "r") as f:
            query = json.load(f)

        j = JSONBlast(subjects, query)
        j.quick_blastn_short()
        alignments = j.get_perfect()
        print(len(alignments))

    def test_make_from_json(self):
        j = JSONBlast(
            [
                {
                    "bases": "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta",
                    "name": "myseq",
                    "circular": False,
                }
            ],
            {
                "bases": "tggaagggctaattcactcccaaagaagacaagatatccttgatctgtggatctaccacacacaaggctacttccctgattagcagaactacacaccagggccaggggtcagatatccactgacctttggatggtgctacaagctagtaccagttgagccagataaggtagaagaggccaataaaggagagaacaccagcttgttacaccctgtgagcctgcatgggatggatgacccggagagagaagtgttagagtggaggtttgacagccgcctagcatttcatcacgtggcccgagagctgcatccggagtacttcaagaactgctgatatcgagcttgctacaagggactttccgctggggactttccagggaggcgtggcctgggcgggactggggagtggcgagccctcagatcctgcatataagcagctgctttttgcctgtactgggtctctctggttagaccagatctgagcctgggagctctctggctaactagggaacccactgcttaagcctcaataaagcttgccttgagtgcttcaagtagtgtgtgcccgtctgttgtgtgactctggtaactagagatccctcagacccttttagtcagtgtggaaaatctctagcagtggcgcccgaacagggacttgaaagcgaaagggaaaccagaggagctctctcgacgcaggactcggcttgctgaagcgcgcacggcaagaggcgaggggcggcgactggtgagtacgccaaaaattttgactagcggaggctagaaggagagagatgggtgcgagagcgtcagtattaagcgggggagaattagatcgcgatgggaaaaaattcggttaaggccagggggaaagaaaaaatataaattaaaacatatagtatgggcaagcagggagctagaacgattcgcagttaatcctggcctgttagaaacatcagaaggctgtagacaaatactgggacagctacaaccatcccttcagacaggatcagaagaacttagatcattatataatacagtagcaaccctctattgtgtgcatcaaaggatagagataaaagacaccaaggaagctttagacaagatagaggaagagcaaaacaaaagtaagaccaccgcacagcaagcggccggccgctgatcttcagacctggaggaggagatatgagggacaattggagaagtgaattatataaatataaagtagtaaaaattgaaccattaggagtagcacccaccaaggcaaagagaagagtggtgcagagagaaaaaagagcagtgggaataggagctttgttccttgggttcttgggagcagcaggaagcactatgggcgcagcgtcaatgacgctgacggtacaggccagacaattattgtctggtatagtgcagcagcagaacaatttgctgagggctattgaggcgcaacagcatctgttgcaactcacagtctggggcatcaagcagctccaggcaagaatcctggctgtgaaagatacctaaaggatcaacagctcctggggatttggggttgctctggaaaactcatttgcaccactgctgtgccttggaatgctagttggagtaataaatctctggaacagatttggaatcacacgacctggatggagtgggacagagaaattaacaattacacaagcttaatacactccttaattgaagaatcgcaaaaccagcaagaaaagaatgaacaagaattattggaattagataaatgggcaagtttgtggaattggtttaacataacaaattggctgtggtatataaaattattcataatgatagtaggaggcttggtaggtttaagaatagtttttgctgtactttctatagtgaatagagttaggcagggatattcaccattatcgtttcagacccacctcccaaccccgaggggacccgacaggcccgaaggaatagaagaagaaggtggagagagagacagagacagatccattcgattagtgaacggatctcgacggtatcgccaaatggcagtattcatccacaattttaaaagaaaaggggggattggggggtacagtgcaggggaaagaatagtagacataatagcaacagacatacaaactaaagaattacaaaaacaaattacaaaaattcaaaattttcgggtttattacagggacagcagagatccagtttggatcgataagcttgatatcgaattcctgcagccccgataaaataaaagattttatttagtctccagaaaaaggggggaatgaaagaccccacctgtaggtttggcaagctagctgcagtaacgccattttgcaaggcatggaaaaataccaaaccaagaatagagaagttcagatcaagggcgggtacatgaaaatagctaacgttgggccaaacaggatatctgcggtgagcagtttcggccccggcccggggccaagaacagatggtcaccgcagtttcggccccggcccgaggccaagaacagatggtccccagatatggcccaaccctcagcagtttcttaagacccatcagatgtttccaggctcccccaaggacctgaaatgaccctgcgccttatttgaattaaccaatcagcctgcttctcgcttctgttcgcgcgcttctgcttcccgagctctataaaagagctcacaacccctcactcggcgcgccagtcctccgacagactgagtcgcccgggggggatctggagctctcgagaattctcacgcgtctgcaggatatcaagcttgcggtaccgcgggcccggccaccatggacaagaagtacagcatcggcctggccatcggcaccaactctgtgggctgggccgtgatcaccgacgagtacaaggtgcccagcaagaaattcaaggtgctgggcaacaccgaccggcacagcatcaagaagaacctgatcggcgccctgctgttcgacagcggagaaacagccgaggccacccggctgaagagaaccgccagaagaagatacaccagacggaagaaccggatctgctatctgcaagagatcttcagcaacgagatggccaaggtggacgacagcttcttccacagactggaagagtccttcctggtggaagaggataagaagcacgagcggcaccccatcttcggcaacatcgtggacgaggtggcctaccacgagaagtaccccaccatctaccacctgagaaagaaactggtggacagcaccgacaaggccgacctgcggctgatctatctggccctggcccacatgatcaagttccggggccacttcctgatcgagggcgacctgaaccccgacaacagcgacgtggacaagctgttcatccagctggtgcagacctacaaccagctgttcgaggaaaaccccatcaacgccagcggcgtggacgccaaggccatcctgtctgccagactgagcaagagcagacggctggaaaatctgatcgcccagctgcccggcgagaagaagaatggcctgttcggcaacctgattgccctgagcctgggcctgacccccaacttcaagagcaacttcgacctggccgaggatgccaaactgcagctgagcaaggacacctacgacgacgacctggacaacctgctggcccagatcggcgaccagtacgccgacctgtttctggccgccaagaacctgtccgacgccatcctgctgagcgacatcctgagagtgaacaccgagatcaccaaggcccccctgagcgcctctatgatcaagagatacgacgagcaccaccaggacctgaccctgctgaaagctctcgtgcggcagcagctgcctgagaagtacaaagagattttcttcgaccagagcaagaacggctacgccggctacatcgatggcggagccagccaggaagagttctacaagttcatcaagcccatcctggaaaagatggacggcaccgaggaactgctcgtgaagctgaacagagaggacctgctgcggaagcagcggaccttcgacaacggcagcatcccccaccagatccacctgggagagctgcacgccattctgcggcggcaggaagatttttacccattcctgaaggacaaccgggaaaagatcgagaagatcctgaccttccgcatcccctactacgtgggccctctggccaggggaaacagcagattcgcctggatgaccagaaagagcgaggaaaccatcaccccctggaacttcgaggaagtggtggacaagggcgccagcgcccagagcttcatcgagcggatgaccaacttcgataagaacctgcccaacgagaaggtgctgcccaagcacagcctgctgtacgagtacttcaccgtgtacaacgagctgaccaaagtgaaatacgtgaccgagggaatgagaaagcccgccttcctgagcggcgagcagaaaaaagccatcgtggacctgctgttcaagaccaaccggaaagtgaccgtgaagcagctgaaagaggactacttcaagaaaatcgagtgcttcgactccgtggaaatctccggcgtggaagatcggttcaacgcctccctgggcacataccacgatctgctgaaaattatcaaggacaaggacttcctggacaatgaggaaaacgaggacattctggaagatatcgtgctgaccctgacactgtttgaggacagagagatgatcgaggaacggctgaaaacctatgcccacctgttcgacgacaaagtgatgaagcagctgaagcggcggagatacaccggctggggcaggctgagccggaagctgatcaacggcatccgggacaagcagtccggcaagacaatcctggatttcctgaagtccgacggcttcgccaacagaaacttcatgcagctgatccacgacgacagcctgacctttaaagaggacatccagaaagcccaggtgtccggccagggcgatagcctgcacgagcacattgccaatctggccggcagccccgccattaagaagggcatcctgcagacagtgaaggtggtggacgagctcgtgaaagtgatgggccggcacaagcccgagaacatcgtgatcgaaatggccagagagaaccagaccacccagaagggacagaagaacagccgcgagagaatgaagcggatcgaagagggcatcaaagagctgggcagccagatcctgaaagaacaccccgtggaaaacacccagctgcagaacgagaagctgtacctgtactacctgcagaatgggcgggatatgtacgtggaccaggaactggacatcaaccggctgtccgactacgatgtggacgctatcgtgcctcagagctttctgaaggacgactccatcgataacaaagtgctgactcggagcgacaagaaccggggcaagagcgacaacgtgccctccgaagaggtcgtgaagaagatgaagaactactggcgccagctgctgaatgccaagctgattacccagaggaagttcgacaatctgaccaaggccgagagaggcggcctgagcgaactggataaggccggcttcatcaagagacagctggtggaaacccggcagatcacaaagcacgtggcacagatcctggactcccggatgaacactaagtacgacgagaacgacaaactgatccgggaagtgaaagtgatcaccctgaagtccaagctggtgtccgatttccggaaggatttccagttttacaaagtgcgcgagatcaacaactaccaccacgcccacgacgcctacctgaacgccgtcgtgggaaccgccctgatcaaaaagtaccctaagctggaaagcgagttcgtgtacggcgactacaaggtgtacgacgtgcggaagatgatcgccaagagcgagcaggaaatcggcaaggctaccgccaagtacttcttctacagcaacatcatgaactttttcaagaccgagattaccctggccaacggcgagatccggaagcggcctctgatcgagacaaacggcgaaacaggcgagatcgtgtgggataagggccgggactttgccaccgtgcggaaagtgctgtctatgccccaagtgaatatcgtgaaaaagaccgaggtgcagacaggcggcttcagcaaagagtctatcctgcccaagaggaacagcgacaagctgatcgccagaaagaaggactgggaccctaagaagtacggcggcttcgacagccccaccgtggcctattctgtgctggtggtggccaaagtggaaaagggcaagtccaagaaactgaagagtgtgaaagagctgctggggatcaccatcatggaaagaagcagcttcgagaagaatcccatcgactttctggaagccaagggctacaaagaagtgaaaaaggacctgatcatcaagctgcctaagtactccctgttcgagctggaaaacggccggaagagaatgctggcctctgccggcgaactgcagaagggaaacgaactggccctgccctccaaatatgtgaacttcctgtacctggccagccactatgagaagctgaagggctcccccgaggataatgagcagaaacagctgtttgtggaacagcacaaacactacctggacgagatcatcgagcagatcagcgagttctccaagagagtgatcctggccgacgctaatctggacaaggtgctgagcgcctacaacaagcacagagacaagcctatcagagagcaggccgagaatatcatccacctgtttaccctgaccaatctgggagcccctgccgccttcaagtactttgacaccaccatcgaccggaagaggtacaccagcaccaaagaggtgctggacgccaccctgatccaccagagcatcaccggcctgtacgagacacggatcgacctgtctcagctgggaggcgacgcctatccctatgacgtgcccgattatgccagcctgggcagcggctcccccaagaaaaaacgcaaggtggaagatcctaagaaaaagcggaaagtggacggcattggtagtgggagcaacggcagcagcggatccagcgagctgattaaggagaacatgcacatgaagctgtacatggagggcaccgtggacaaccatcacttcaagtgcacatccgagggcgaaggcaagccctacgagggcacccagaccatgagaatcaaggtggtcgagggcggccctctccccttcgccttcgacatcctggctactagcttcctctacggcagcaagaccttcatcaaccacacccagggcatccccgacttcttcaagcagtccttccctgagggcttcacatgggagagagtcaccacatacgaagacgggggcgtgctgaccgctacccaggacaccagcctccaggacggctgcctcatctacaacgtcaagatcagaggggtgaacttcacatccaacggccctgtgatgcagaagaaaacactcggctgggaggccttcaccgagacgctgtaccccgctgacggcggcctggaaggcagaaacgacatggccctgaagctcgtgggcgggagccatctgatcgcaaacatcaagaccacatatagatccaagaaacccgctaagaacctcaagatgcctggcgtctactatgtggactacagactggaaagaatcaaggaggccaacaacgagacctacgtcgagcagcacgaggtggcagtggccagatactgcgacctccctagcaaactggggcacaagcttaattagtaaggccgcgactctagagtcgacctgcaggcatgcaagcttgatatcaagcttatcgataatcaacctctggattacaaaatttgtgaaagattgactggtattcttaactatgttgctccttttacgctatgtggatacgctgctttaatgcctttgtatcatgctattgcttcccgtatggctttcattttctcctccttgtataaatcctggttgctgtctctttatgaggagttgtggcccgttgtcaggcaacgtggcgtggtgtgcactgtgtttgctgacgcaacccccactggttggggcattgccaccacctgtcagctcctttccgggactttcgctttccccctccctattgccacggcggaactcatcgccgcctgccttgcccgctgctggacaggggctcggctgttgggcactgacaattccgtggtgttgtcggggaaatcatcgtcctttccttggctgctcgcctgtgttgccacctggattctgcgcgggacgtccttctgctacgtcccttcggccctcaatccagcggaccttccttcccgcggcctgctgccggctctgcggcctcttccgcgtcttcgccttcgccctcagacgagtcggatctccctttgggccgcctccccgcatcgataccgtcgacctcgagggaattaattcgagctcggtacctttaagaccaatgacttacaaggcagctgtagatcttagccactttttaaaagaaaaggggggactggaagggctaattcactcccaacgaagacaagatctgctttttgcttgtactgggtctctctggttagaccagatctgagcctgggagctctctggctaactagggaacccactgcttaagcctcaataaagcttgccttgagtgcttcaagtagtgtgtgcccgtctgttgtgtgactctggtaactagagatccctcagacccttttagtcagtgtggaaaatctctagcagcatctagaattaattccgtgtattctatagtgtcacctaaatcgtatgtgtatgatacataaggttatgtattaattgtagccgcgttctaacgacaatatgtacaagcctaattgtgtagcatctggcttactgaagcagaccctatcatctctctcgtaaactgccgtcagagtcggtttggttggacgaaccttctgagtttctggtaacgccgtcccgcacccggaaatggtcagcgaaccaatcagcagggtcatcgctagccagatcctctacgccggacgcatcgtggccggcatcaccggcgccacaggtgcggttgctggcgcctatatcgccgacatcaccgatggggaagatcgggctcgccacttcgggctcatgagcgcttgtttcggcgtgggtatggtggcaggccccgtggccgggggactgttgggcgccatctccttgcatgcaccattccttgcggcggcggtgctcaacggcctcaacctactactgggctgcttcctaatgcaggagtcgcataagggagagcgtcgaatggtgcactctcagtacaatctgctctgatgccgcatagttaagccagccccgacacccgccaacacccgctgacgcgccctgacgggcttgtctgctcccggcatccgcttacagacaagctgtgaccgtctccgggagctgcatgtgtcagaggttttcaccgtcatcaccgaaacgcgcgagacgaaagggcctcgtgatacgcctatttttataggttaatgtcatgataataatggtttcttagacgtcaggtggcacttttcggggaaatgtgcgcggaacccctatttgtttatttttctaaatacattcaaatatgtatccgctcatgagacaataaccctgataaatgcttcaataatattgaaaaaggaagagtatgagtattcaacatttccgtgtcgcccttattcccttttttgcggcattttgccttcctgtttttgctcacccagaaacgctggtgaaagtaaaagatgctgaagatcagttgggtgcacgagtgggttacatcgaactggatctcaacagcggtaagatccttgagagttttcgccccgaagaacgttttccaatgatgagcacttttaaagttctgctatgtggcgcggtattatcccgtattgacgccgggcaagagcaactcggtcgccgcatacactattctcagaatgacttggttgagtactcaccagtcacagaaaagcatcttacggatggcatgacagtaagagaattatgcagtgctgccataaccatgagtgataacactgcggccaacttacttctgacaacgatcggaggaccgaaggagctaaccgcttttttgcacaacatgggggatcatgtaactcgccttgatcgttgggaaccggagctgaatgaagccataccaaacgacgagcgtgacaccacgatgcctgtagcaatggcaacaacgttgcgcaaactattaactggcgaactacttactctagcttcccggcaacaattaatagactggatggaggcggataaagttgcaggaccacttctgcgctcggcccttccggctggctggtttattgctgataaatctggagccggtgagcgtgggtctcgcggtatcattgcagcactggggccagatggtaagccctcccgtatcgtagttatctacacgacggggagtcaggcaactatggatgaacgaaatagacagatcgctgagataggtgcctcactgattaagcattggtaactgtcagaccaagtttactcatatatactttagattgatttaaaacttcatttttaatttaaaaggatctaggtgaagatcctttttgataatctcatgaccaaaatcccttaacgtgagttttcgttccactgagcgtcagaccccgtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaacaaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttccgaaggtaactggcttcagcagagcgcagataccaaatactgtctttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggactcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaactgagatacctacagcgtgagctatgagaaagcgccacgcttcccgaagggagaaaggcggacaggtatccggtaagcggcagggtcggaacaggagagcgcacgagggagcttccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccacctctgacttgagcgtcgatttttgtgatgctcgtcaggggggcggagcctatggaaaaacgccagcaacgcggcctttttacggttcctggccttttgctggccttttgctcacatgttctttcctgcgttatcccctgattctgtggataaccgtattaccgcctttgagtgagctgataccgctcgccgcagccgaacgaccgagcgcagcgagtcagtgagcgaggaagcggaagagcgcccaatacgcaaaccgcctctccccgcgcgttggccgattcattaatgcagctgtggaatgtgtgtcagttagggtgtggaaagtccccaggctccccagcaggcagaagtatgcaaagcatgcatctcaattagtcagcaaccaggtgtggaaagtccccaggctccccagcaggcagaagtatgcaaagcatgcatctcaattagtcagcaaccatagtcccgcccctaactccgcccatcccgcccctaactccgcccagttccgcccattctccgccccatggctgactaattttttttatttatgcagaggccgaggccgcctcggcctctgagctattccagaagtagtgaggaggcttttttggaggcctaggcttttgcaaaaagcttggacacaagacaggcttgcgagatatgtttgagaataccactttatcccgcgtcagggagaggcagtgcgtaaaaagacgcggactcatgtgaaatactggtttttagtgcgccagatctctataatctcgcgcaacctattttcccctcgaacactttttaagccgtagataaacaggctgggacacttcacatgagcgaaaaatacatcgtcacctgggacatgttgcagatccatgcacgtaaactcgcaagccgactgatgccttctgaacaatggaaaggcattattgccgtaagccgtggcggtctgtaccgggtgcgttactggcgcgtgaactgggtattcgtcatgtcgataccgtttgtatttccagctacgatcacgacaaccagcgcgagcttaaagtgctgaaacgcgcagaaggcgatggcgaaggcttcatcgttattgatgacctggtggataccggtggtactgcggttgcgattcgtgaaatgtatccaaaagcgcactttgtcaccatcttcgcaaaaccggctggtcgtccgctggttgatgactatgttgttgatatcccgcaagatacctggattgaacagccgtgggatatgggcgtcgtattcgtcccgccaatctccggtcgctaatcttttcaacgcctggcactgccgggcgttgttctttttaacttcaggcgggttacaatagtttccagtaagtattctggaggctgcatccatgacacaggcaaacctgagcgaaaccctgttcaaaccccgctttaaacatcctgaaacctcgacgctagtccgccgctttaatcacggcgcacaaccgcctgtgcagtcggcccttgatggtaaaaccatccctcactggtatcgcatgattaaccgtctgatgtggatctggcgcggcattgacccacgcgaaatcctcgacgtccaggcacgtattgtgatgagcgatgccgaacgtaccgacgatgatttatacgatacggtgattggctaccgtggcggcaactggatttatgagtgggccccggatctttgtgaaggaaccttacttctgtggtgtgacataattggacaaactacctacagagatttaaagctctaaggtaaatataaaatttttaagtgtataatgtgttaaactactgattctaattgtttgtgtattttagattccaacctatggaactgatgaatgggagcagtggtggaatgcctttaatgaggaaaacctgttttgctcagaagaaatgccatctagtgatgatgaggctactgctgactctcaacattctactcctccaaaaaagaagagaaaggtagaagaccccaaggactttccttcagaattgctaagttttttgagtcatgctgtgtttagtaatagaactcttgcttgctttgctatttacaccacaaaggaaaaagctgcactgctatacaagaaaattatggaaaaatattctgtaacctttataagtaggcataacagttataatcataacatactgttttttcttactccacacaggcatagagtgtctgctattaataactatgctcaaaaattgtgtacctttagctttttaatttgtaaaggggttaataaggaatatttgatgtatagtgccttgactagagatcataatcagccataccacatttgtagaggttttacttgctttaaaaaacctcccacacctccccctgaacctgaaacataaaatgaatgcaattgttgttgttaacttgtttattgcagcttataatggttacaaataaagcaatagcatcacaaatttcacaaataaagcatttttttcactgcattctagttgtggtttgtccaaactcatcaatgtatcttatcatgtctggatcaactggataactcaagctaaccaaaatcatcccaaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcattttaaagaaattgtatttgttaaatatgtactacaaacttagtagt",
                "name": "myseq2",
                "circular": False,
            },
        )
        j.quick_blastn()
        results = j.results
        print(results)

    def test_make_from_json_with_preloaded_id(self):
        """We expect the alignment results to contain the very same id and
        information in the query and subject alignments as the sequences that
        went into the JSONBlast alignment.

        Since this is a perfect alignment, we also expect the sequences
        to be spit back out.
        """
        j = JSONBlast(
            [
                {
                    "id": "ABCDEFG",
                    "bases": "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta",
                    "name": "myseq",
                    "circular": False,
                }
            ],
            {
                "id": "1234",
                "bases": "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta",
                "name": "myseq2",
                "circular": False,
            },
        )
        j.quick_blastn()
        results = j.results

        expected_query = {
            "bases": "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta".upper(),
            "name": "myseq2",
            "circular": False,
            "length": 93,
            "start": 1,
            "end": 93,
            "raw_end": 93,
            "origin_record_id": "1234",
            "origin_sequence_length": 93,
            "strand": 1,
        }

        expected_subject = {
            "bases": "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta".upper(),
            "name": "myseq",
            "circular": False,
            "length": 93,
            "start": 1,
            "end": 93,
            "raw_end": 93,
            "strand": 1,
            "origin_record_id": "ABCDEFG",
            "origin_sequence_length": 93,
        }
        del results[0]["query"]["sequence_id"]
        del results[0]["subject"]["sequence_id"]
        del results[0]["query"]["origin_key"]
        del results[0]["subject"]["origin_key"]

        assert results[0]["query"] == expected_query
        assert results[0]["subject"] == expected_subject

    def test_json_blast_with_orderdict(self):
        """JSONBlast should handle the subclasses of dict as inputs as well as
        dict."""
        from collections import OrderedDict

        j = JSONBlast(
            [
                OrderedDict(
                    {
                        "id": "ABCDEFG",
                        "bases": "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta",
                        "name": "myseq",
                        "circular": False,
                    }
                )
            ],
            OrderedDict(
                {
                    "id": "1234",
                    "bases": "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta",
                    "name": "myseq2",
                    "circular": False,
                }
            ),
        )
        j.quick_blastn()
        results = j.results

        expected_query = {
            "bases": "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta".upper(),
            "name": "myseq2",
            "circular": False,
            "length": 93,
            "start": 1,
            "end": 93,
            "raw_end": 93,
            "origin_record_id": "1234",
            "origin_sequence_length": 93,
            "strand": 1,
        }

        expected_subject = {
            "bases": "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta".upper(),
            "name": "myseq",
            "circular": False,
            "length": 93,
            "start": 1,
            "end": 93,
            "raw_end": 93,
            "strand": 1,
            "origin_record_id": "ABCDEFG",
            "origin_sequence_length": 93,
        }

        del results[0]["query"]["sequence_id"]
        del results[0]["subject"]["sequence_id"]
        del results[0]["query"]["origin_key"]
        del results[0]["subject"]["origin_key"]

        assert results[0]["query"] == expected_query
        assert results[0]["subject"] == expected_subject

    def test_json_blast_with_no_hits(self):
        """JSONBlast should handle the subclasses of dict as inputs as well as
        dict."""
        from collections import OrderedDict

        j = JSONBlast(
            [
                OrderedDict(
                    {
                        "id": "ABCDEFG",
                        "bases": "aggagagagagaggagggagagaagaggagagagagaga",
                        "name": "myseq",
                        "circular": False,
                    }
                )
            ],
            OrderedDict(
                {
                    "id": "1234",
                    "bases": "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta",
                    "name": "myseq2",
                    "circular": False,
                }
            ),
        )
        j.quick_blastn()
        results = j.results
        print(results)

    def test_json_blast_with_repeats(self):
        """We expect no hits from repeats."""
        from collections import OrderedDict

        j = JSONBlast(
            [
                OrderedDict(
                    {
                        "id": "ABCDEFG",
                        "bases": "agtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtag",
                        "name": "myseq",
                        "circular": False,
                    }
                )
            ],
            OrderedDict(
                {
                    "id": "1234",
                    "bases": "agtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtag",
                    "name": "myseq2",
                    "circular": False,
                }
            ),
        )
        j.quick_blastn()
        results = j.results
        assert results == []

    def test_json_blast_with_objects_raises_pyblasterror(self):
        class Sequence:
            def __init__(self, **kwargs):
                self.__dict__.update(**kwargs)

        query = Sequence(
            **{
                "id": "ABCDEFG",
                "bases": "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta",
                "name": "myseq",
                "circular": False,
            }
        )

        from copy import copy

        subject = copy(query)

        with pytest.raises(PyBlastException):
            j = JSONBlast([subject], [query])
            j.quick_blastn()
            assert len(j.results) == 1


class TestJSONBlastForExpectedSequences:
    @pytest.fixture
    def frag(self):
        return "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta"

    @pytest.fixture
    def seqs(self, frag):
        """Alignments that wraps around the query more than onces."""

        # class Sequence(object):
        #
        #     def __init__(self, **kwargs):
        #         self.__dict__.update(**kwargs)

        junk1 = "atgctatgctgatgctgctgtgctgatgctgatgtgtattgctgtatcgcgcgagttagc"
        junk2 = "g" * 30

        query = {"id": "myquery", "bases": frag, "name": "myseq", "circular": False}

        subject = {
            "id": "mysubject",
            "bases": junk1 + frag + junk2,
            "name": "myseq",
            "circular": False,
        }

        return [subject], query

    @pytest.fixture
    def alignments(self, seqs):
        subjects, query = seqs
        j = JSONBlast(subjects, query, span_origin=False)
        j.quick_blastn()
        results = j.get_perfect()
        return results

    @pytest.fixture
    def alignments_rc_subject(self, seqs):
        subjects, query = seqs
        subjects[0]["bases"] = reverse_complement(subjects[0]["bases"])
        j = JSONBlast(subjects, query, span_origin=False)
        j.quick_blastn()
        results = j.get_perfect()
        return results

    def test_expected_query_sequence(self, alignments, frag):
        assert len(alignments) == 1
        alignment = alignments[0]
        assert alignment["query"]["strand"] == 1
        assert alignment["query"]["bases"].upper() == frag.upper()

    def test_expected_subject_sequence(self, alignments, frag):
        assert len(alignments) == 1
        alignment = alignments[0]
        assert alignment["subject"]["strand"] == 1
        assert alignment["subject"]["bases"].upper() == frag.upper()

    def test_expected_length(self, alignments):
        for align in alignments:
            assert align["meta"]["alignment length"] == len(align["query"]["bases"])
            assert align["meta"]["alignment length"] == len(align["subject"]["bases"])
            assert align["query"]["end"] - align["query"]["start"] + 1 == len(
                align["query"]["bases"]
            )
            assert align["subject"]["end"] - align["subject"]["start"] + 1 == len(
                align["subject"]["bases"]
            )

    def test_expected_query_sequence_rc_subject(self, alignments_rc_subject, frag):
        assert len(alignments_rc_subject) == 1
        alignment = alignments_rc_subject[0]
        assert alignment["query"]["strand"] == 1
        assert alignment["query"]["bases"].upper() == frag.upper()

    def test_expected_subject_sequence_rc_subject(self, alignments_rc_subject, frag):
        assert len(alignments_rc_subject) == 1
        alignment = alignments_rc_subject[0]
        assert alignment["subject"]["strand"] == -1
        assert alignment["subject"]["bases"].upper() == frag.upper()
