"""
testing expected bases for query and subject for JSON Blast
"""

import pytest
import os
from pyblast.blast import JSONBlast
from pyblast.exceptions import PyBlastException
from pyblast.utils import reverse_complement
import json


class TestJSONBlast:

    @pytest.fixture
    def aligner(self, here):
        with open(os.path.join(here, "data/test_data/templates.json"), 'r') as f:
            subjects = json.load(f)
        with open(os.path.join(here, "data/test_data/query.json"), 'r') as f:
            query = json.load(f)
        j = JSONBlast(query, subjects)
        j.quick_blastn()
        return j

    def test_example(self):
        j = JSONBlast.use_test_data()
        j.quick_blastn()
        pass

    def test_perfect_matches(self, here):
        with open(os.path.join(here, "data/test_data/primer.json"), 'r') as f:
            subjects = json.load(f)
        with open(os.path.join(here, "data/test_data/query.json"), 'r') as f:
            query = json.load(f)

        j = JSONBlast(subjects, query)
        j.quick_blastn_short()
        alignments = j.results.get_perfect().alignments
        print(len(j.results.alignments))
        print(len(alignments))

    def test_make_from_json(self):
        j = JSONBlast([{
                           "sequence": "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta",
                           "name": "myseq", "circular": False}], {
                          "sequence": "tggaagggctaattcactcccaaagaagacaagatatccttgatctgtggatctaccacacacaaggctacttccctgattagcagaactacacaccagggccaggggtcagatatccactgacctttggatggtgctacaagctagtaccagttgagccagataaggtagaagaggccaataaaggagagaacaccagcttgttacaccctgtgagcctgcatgggatggatgacccggagagagaagtgttagagtggaggtttgacagccgcctagcatttcatcacgtggcccgagagctgcatccggagtacttcaagaactgctgatatcgagcttgctacaagggactttccgctggggactttccagggaggcgtggcctgggcgggactggggagtggcgagccctcagatcctgcatataagcagctgctttttgcctgtactgggtctctctggttagaccagatctgagcctgggagctctctggctaactagggaacccactgcttaagcctcaataaagcttgccttgagtgcttcaagtagtgtgtgcccgtctgttgtgtgactctggtaactagagatccctcagacccttttagtcagtgtggaaaatctctagcagtggcgcccgaacagggacttgaaagcgaaagggaaaccagaggagctctctcgacgcaggactcggcttgctgaagcgcgcacggcaagaggcgaggggcggcgactggtgagtacgccaaaaattttgactagcggaggctagaaggagagagatgggtgcgagagcgtcagtattaagcgggggagaattagatcgcgatgggaaaaaattcggttaaggccagggggaaagaaaaaatataaattaaaacatatagtatgggcaagcagggagctagaacgattcgcagttaatcctggcctgttagaaacatcagaaggctgtagacaaatactgggacagctacaaccatcccttcagacaggatcagaagaacttagatcattatataatacagtagcaaccctctattgtgtgcatcaaaggatagagataaaagacaccaaggaagctttagacaagatagaggaagagcaaaacaaaagtaagaccaccgcacagcaagcggccggccgctgatcttcagacctggaggaggagatatgagggacaattggagaagtgaattatataaatataaagtagtaaaaattgaaccattaggagtagcacccaccaaggcaaagagaagagtggtgcagagagaaaaaagagcagtgggaataggagctttgttccttgggttcttgggagcagcaggaagcactatgggcgcagcgtcaatgacgctgacggtacaggccagacaattattgtctggtatagtgcagcagcagaacaatttgctgagggctattgaggcgcaacagcatctgttgcaactcacagtctggggcatcaagcagctccaggcaagaatcctggctgtgaaagatacctaaaggatcaacagctcctggggatttggggttgctctggaaaactcatttgcaccactgctgtgccttggaatgctagttggagtaataaatctctggaacagatttggaatcacacgacctggatggagtgggacagagaaattaacaattacacaagcttaatacactccttaattgaagaatcgcaaaaccagcaagaaaagaatgaacaagaattattggaattagataaatgggcaagtttgtggaattggtttaacataacaaattggctgtggtatataaaattattcataatgatagtaggaggcttggtaggtttaagaatagtttttgctgtactttctatagtgaatagagttaggcagggatattcaccattatcgtttcagacccacctcccaaccccgaggggacccgacaggcccgaaggaatagaagaagaaggtggagagagagacagagacagatccattcgattagtgaacggatctcgacggtatcgccaaatggcagtattcatccacaattttaaaagaaaaggggggattggggggtacagtgcaggggaaagaatagtagacataatagcaacagacatacaaactaaagaattacaaaaacaaattacaaaaattcaaaattttcgggtttattacagggacagcagagatccagtttggatcgataagcttgatatcgaattcctgcagccccgataaaataaaagattttatttagtctccagaaaaaggggggaatgaaagaccccacctgtaggtttggcaagctagctgcagtaacgccattttgcaaggcatggaaaaataccaaaccaagaatagagaagttcagatcaagggcgggtacatgaaaatagctaacgttgggccaaacaggatatctgcggtgagcagtttcggccccggcccggggccaagaacagatggtcaccgcagtttcggccccggcccgaggccaagaacagatggtccccagatatggcccaaccctcagcagtttcttaagacccatcagatgtttccaggctcccccaaggacctgaaatgaccctgcgccttatttgaattaaccaatcagcctgcttctcgcttctgttcgcgcgcttctgcttcccgagctctataaaagagctcacaacccctcactcggcgcgccagtcctccgacagactgagtcgcccgggggggatctggagctctcgagaattctcacgcgtctgcaggatatcaagcttgcggtaccgcgggcccggccaccatggacaagaagtacagcatcggcctggccatcggcaccaactctgtgggctgggccgtgatcaccgacgagtacaaggtgcccagcaagaaattcaaggtgctgggcaacaccgaccggcacagcatcaagaagaacctgatcggcgccctgctgttcgacagcggagaaacagccgaggccacccggctgaagagaaccgccagaagaagatacaccagacggaagaaccggatctgctatctgcaagagatcttcagcaacgagatggccaaggtggacgacagcttcttccacagactggaagagtccttcctggtggaagaggataagaagcacgagcggcaccccatcttcggcaacatcgtggacgaggtggcctaccacgagaagtaccccaccatctaccacctgagaaagaaactggtggacagcaccgacaaggccgacctgcggctgatctatctggccctggcccacatgatcaagttccggggccacttcctgatcgagggcgacctgaaccccgacaacagcgacgtggacaagctgttcatccagctggtgcagacctacaaccagctgttcgaggaaaaccccatcaacgccagcggcgtggacgccaaggccatcctgtctgccagactgagcaagagcagacggctggaaaatctgatcgcccagctgcccggcgagaagaagaatggcctgttcggcaacctgattgccctgagcctgggcctgacccccaacttcaagagcaacttcgacctggccgaggatgccaaactgcagctgagcaaggacacctacgacgacgacctggacaacctgctggcccagatcggcgaccagtacgccgacctgtttctggccgccaagaacctgtccgacgccatcctgctgagcgacatcctgagagtgaacaccgagatcaccaaggcccccctgagcgcctctatgatcaagagatacgacgagcaccaccaggacctgaccctgctgaaagctctcgtgcggcagcagctgcctgagaagtacaaagagattttcttcgaccagagcaagaacggctacgccggctacatcgatggcggagccagccaggaagagttctacaagttcatcaagcccatcctggaaaagatggacggcaccgaggaactgctcgtgaagctgaacagagaggacctgctgcggaagcagcggaccttcgacaacggcagcatcccccaccagatccacctgggagagctgcacgccattctgcggcggcaggaagatttttacccattcctgaaggacaaccgggaaaagatcgagaagatcctgaccttccgcatcccctactacgtgggccctctggccaggggaaacagcagattcgcctggatgaccagaaagagcgaggaaaccatcaccccctggaacttcgaggaagtggtggacaagggcgccagcgcccagagcttcatcgagcggatgaccaacttcgataagaacctgcccaacgagaaggtgctgcccaagcacagcctgctgtacgagtacttcaccgtgtacaacgagctgaccaaagtgaaatacgtgaccgagggaatgagaaagcccgccttcctgagcggcgagcagaaaaaagccatcgtggacctgctgttcaagaccaaccggaaagtgaccgtgaagcagctgaaagaggactacttcaagaaaatcgagtgcttcgactccgtggaaatctccggcgtggaagatcggttcaacgcctccctgggcacataccacgatctgctgaaaattatcaaggacaaggacttcctggacaatgaggaaaacgaggacattctggaagatatcgtgctgaccctgacactgtttgaggacagagagatgatcgaggaacggctgaaaacctatgcccacctgttcgacgacaaagtgatgaagcagctgaagcggcggagatacaccggctggggcaggctgagccggaagctgatcaacggcatccgggacaagcagtccggcaagacaatcctggatttcctgaagtccgacggcttcgccaacagaaacttcatgcagctgatccacgacgacagcctgacctttaaagaggacatccagaaagcccaggtgtccggccagggcgatagcctgcacgagcacattgccaatctggccggcagccccgccattaagaagggcatcctgcagacagtgaaggtggtggacgagctcgtgaaagtgatgggccggcacaagcccgagaacatcgtgatcgaaatggccagagagaaccagaccacccagaagggacagaagaacagccgcgagagaatgaagcggatcgaagagggcatcaaagagctgggcagccagatcctgaaagaacaccccgtggaaaacacccagctgcagaacgagaagctgtacctgtactacctgcagaatgggcgggatatgtacgtggaccaggaactggacatcaaccggctgtccgactacgatgtggacgctatcgtgcctcagagctttctgaaggacgactccatcgataacaaagtgctgactcggagcgacaagaaccggggcaagagcgacaacgtgccctccgaagaggtcgtgaagaagatgaagaactactggcgccagctgctgaatgccaagctgattacccagaggaagttcgacaatctgaccaaggccgagagaggcggcctgagcgaactggataaggccggcttcatcaagagacagctggtggaaacccggcagatcacaaagcacgtggcacagatcctggactcccggatgaacactaagtacgacgagaacgacaaactgatccgggaagtgaaagtgatcaccctgaagtccaagctggtgtccgatttccggaaggatttccagttttacaaagtgcgcgagatcaacaactaccaccacgcccacgacgcctacctgaacgccgtcgtgggaaccgccctgatcaaaaagtaccctaagctggaaagcgagttcgtgtacggcgactacaaggtgtacgacgtgcggaagatgatcgccaagagcgagcaggaaatcggcaaggctaccgccaagtacttcttctacagcaacatcatgaactttttcaagaccgagattaccctggccaacggcgagatccggaagcggcctctgatcgagacaaacggcgaaacaggcgagatcgtgtgggataagggccgggactttgccaccgtgcggaaagtgctgtctatgccccaagtgaatatcgtgaaaaagaccgaggtgcagacaggcggcttcagcaaagagtctatcctgcccaagaggaacagcgacaagctgatcgccagaaagaaggactgggaccctaagaagtacggcggcttcgacagccccaccgtggcctattctgtgctggtggtggccaaagtggaaaagggcaagtccaagaaactgaagagtgtgaaagagctgctggggatcaccatcatggaaagaagcagcttcgagaagaatcccatcgactttctggaagccaagggctacaaagaagtgaaaaaggacctgatcatcaagctgcctaagtactccctgttcgagctggaaaacggccggaagagaatgctggcctctgccggcgaactgcagaagggaaacgaactggccctgccctccaaatatgtgaacttcctgtacctggccagccactatgagaagctgaagggctcccccgaggataatgagcagaaacagctgtttgtggaacagcacaaacactacctggacgagatcatcgagcagatcagcgagttctccaagagagtgatcctggccgacgctaatctggacaaggtgctgagcgcctacaacaagcacagagacaagcctatcagagagcaggccgagaatatcatccacctgtttaccctgaccaatctgggagcccctgccgccttcaagtactttgacaccaccatcgaccggaagaggtacaccagcaccaaagaggtgctggacgccaccctgatccaccagagcatcaccggcctgtacgagacacggatcgacctgtctcagctgggaggcgacgcctatccctatgacgtgcccgattatgccagcctgggcagcggctcccccaagaaaaaacgcaaggtggaagatcctaagaaaaagcggaaagtggacggcattggtagtgggagcaacggcagcagcggatccagcgagctgattaaggagaacatgcacatgaagctgtacatggagggcaccgtggacaaccatcacttcaagtgcacatccgagggcgaaggcaagccctacgagggcacccagaccatgagaatcaaggtggtcgagggcggccctctccccttcgccttcgacatcctggctactagcttcctctacggcagcaagaccttcatcaaccacacccagggcatccccgacttcttcaagcagtccttccctgagggcttcacatgggagagagtcaccacatacgaagacgggggcgtgctgaccgctacccaggacaccagcctccaggacggctgcctcatctacaacgtcaagatcagaggggtgaacttcacatccaacggccctgtgatgcagaagaaaacactcggctgggaggccttcaccgagacgctgtaccccgctgacggcggcctggaaggcagaaacgacatggccctgaagctcgtgggcgggagccatctgatcgcaaacatcaagaccacatatagatccaagaaacccgctaagaacctcaagatgcctggcgtctactatgtggactacagactggaaagaatcaaggaggccaacaacgagacctacgtcgagcagcacgaggtggcagtggccagatactgcgacctccctagcaaactggggcacaagcttaattagtaaggccgcgactctagagtcgacctgcaggcatgcaagcttgatatcaagcttatcgataatcaacctctggattacaaaatttgtgaaagattgactggtattcttaactatgttgctccttttacgctatgtggatacgctgctttaatgcctttgtatcatgctattgcttcccgtatggctttcattttctcctccttgtataaatcctggttgctgtctctttatgaggagttgtggcccgttgtcaggcaacgtggcgtggtgtgcactgtgtttgctgacgcaacccccactggttggggcattgccaccacctgtcagctcctttccgggactttcgctttccccctccctattgccacggcggaactcatcgccgcctgccttgcccgctgctggacaggggctcggctgttgggcactgacaattccgtggtgttgtcggggaaatcatcgtcctttccttggctgctcgcctgtgttgccacctggattctgcgcgggacgtccttctgctacgtcccttcggccctcaatccagcggaccttccttcccgcggcctgctgccggctctgcggcctcttccgcgtcttcgccttcgccctcagacgagtcggatctccctttgggccgcctccccgcatcgataccgtcgacctcgagggaattaattcgagctcggtacctttaagaccaatgacttacaaggcagctgtagatcttagccactttttaaaagaaaaggggggactggaagggctaattcactcccaacgaagacaagatctgctttttgcttgtactgggtctctctggttagaccagatctgagcctgggagctctctggctaactagggaacccactgcttaagcctcaataaagcttgccttgagtgcttcaagtagtgtgtgcccgtctgttgtgtgactctggtaactagagatccctcagacccttttagtcagtgtggaaaatctctagcagcatctagaattaattccgtgtattctatagtgtcacctaaatcgtatgtgtatgatacataaggttatgtattaattgtagccgcgttctaacgacaatatgtacaagcctaattgtgtagcatctggcttactgaagcagaccctatcatctctctcgtaaactgccgtcagagtcggtttggttggacgaaccttctgagtttctggtaacgccgtcccgcacccggaaatggtcagcgaaccaatcagcagggtcatcgctagccagatcctctacgccggacgcatcgtggccggcatcaccggcgccacaggtgcggttgctggcgcctatatcgccgacatcaccgatggggaagatcgggctcgccacttcgggctcatgagcgcttgtttcggcgtgggtatggtggcaggccccgtggccgggggactgttgggcgccatctccttgcatgcaccattccttgcggcggcggtgctcaacggcctcaacctactactgggctgcttcctaatgcaggagtcgcataagggagagcgtcgaatggtgcactctcagtacaatctgctctgatgccgcatagttaagccagccccgacacccgccaacacccgctgacgcgccctgacgggcttgtctgctcccggcatccgcttacagacaagctgtgaccgtctccgggagctgcatgtgtcagaggttttcaccgtcatcaccgaaacgcgcgagacgaaagggcctcgtgatacgcctatttttataggttaatgtcatgataataatggtttcttagacgtcaggtggcacttttcggggaaatgtgcgcggaacccctatttgtttatttttctaaatacattcaaatatgtatccgctcatgagacaataaccctgataaatgcttcaataatattgaaaaaggaagagtatgagtattcaacatttccgtgtcgcccttattcccttttttgcggcattttgccttcctgtttttgctcacccagaaacgctggtgaaagtaaaagatgctgaagatcagttgggtgcacgagtgggttacatcgaactggatctcaacagcggtaagatccttgagagttttcgccccgaagaacgttttccaatgatgagcacttttaaagttctgctatgtggcgcggtattatcccgtattgacgccgggcaagagcaactcggtcgccgcatacactattctcagaatgacttggttgagtactcaccagtcacagaaaagcatcttacggatggcatgacagtaagagaattatgcagtgctgccataaccatgagtgataacactgcggccaacttacttctgacaacgatcggaggaccgaaggagctaaccgcttttttgcacaacatgggggatcatgtaactcgccttgatcgttgggaaccggagctgaatgaagccataccaaacgacgagcgtgacaccacgatgcctgtagcaatggcaacaacgttgcgcaaactattaactggcgaactacttactctagcttcccggcaacaattaatagactggatggaggcggataaagttgcaggaccacttctgcgctcggcccttccggctggctggtttattgctgataaatctggagccggtgagcgtgggtctcgcggtatcattgcagcactggggccagatggtaagccctcccgtatcgtagttatctacacgacggggagtcaggcaactatggatgaacgaaatagacagatcgctgagataggtgcctcactgattaagcattggtaactgtcagaccaagtttactcatatatactttagattgatttaaaacttcatttttaatttaaaaggatctaggtgaagatcctttttgataatctcatgaccaaaatcccttaacgtgagttttcgttccactgagcgtcagaccccgtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaacaaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttccgaaggtaactggcttcagcagagcgcagataccaaatactgtctttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggactcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaactgagatacctacagcgtgagctatgagaaagcgccacgcttcccgaagggagaaaggcggacaggtatccggtaagcggcagggtcggaacaggagagcgcacgagggagcttccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccacctctgacttgagcgtcgatttttgtgatgctcgtcaggggggcggagcctatggaaaaacgccagcaacgcggcctttttacggttcctggccttttgctggccttttgctcacatgttctttcctgcgttatcccctgattctgtggataaccgtattaccgcctttgagtgagctgataccgctcgccgcagccgaacgaccgagcgcagcgagtcagtgagcgaggaagcggaagagcgcccaatacgcaaaccgcctctccccgcgcgttggccgattcattaatgcagctgtggaatgtgtgtcagttagggtgtggaaagtccccaggctccccagcaggcagaagtatgcaaagcatgcatctcaattagtcagcaaccaggtgtggaaagtccccaggctccccagcaggcagaagtatgcaaagcatgcatctcaattagtcagcaaccatagtcccgcccctaactccgcccatcccgcccctaactccgcccagttccgcccattctccgccccatggctgactaattttttttatttatgcagaggccgaggccgcctcggcctctgagctattccagaagtagtgaggaggcttttttggaggcctaggcttttgcaaaaagcttggacacaagacaggcttgcgagatatgtttgagaataccactttatcccgcgtcagggagaggcagtgcgtaaaaagacgcggactcatgtgaaatactggtttttagtgcgccagatctctataatctcgcgcaacctattttcccctcgaacactttttaagccgtagataaacaggctgggacacttcacatgagcgaaaaatacatcgtcacctgggacatgttgcagatccatgcacgtaaactcgcaagccgactgatgccttctgaacaatggaaaggcattattgccgtaagccgtggcggtctgtaccgggtgcgttactggcgcgtgaactgggtattcgtcatgtcgataccgtttgtatttccagctacgatcacgacaaccagcgcgagcttaaagtgctgaaacgcgcagaaggcgatggcgaaggcttcatcgttattgatgacctggtggataccggtggtactgcggttgcgattcgtgaaatgtatccaaaagcgcactttgtcaccatcttcgcaaaaccggctggtcgtccgctggttgatgactatgttgttgatatcccgcaagatacctggattgaacagccgtgggatatgggcgtcgtattcgtcccgccaatctccggtcgctaatcttttcaacgcctggcactgccgggcgttgttctttttaacttcaggcgggttacaatagtttccagtaagtattctggaggctgcatccatgacacaggcaaacctgagcgaaaccctgttcaaaccccgctttaaacatcctgaaacctcgacgctagtccgccgctttaatcacggcgcacaaccgcctgtgcagtcggcccttgatggtaaaaccatccctcactggtatcgcatgattaaccgtctgatgtggatctggcgcggcattgacccacgcgaaatcctcgacgtccaggcacgtattgtgatgagcgatgccgaacgtaccgacgatgatttatacgatacggtgattggctaccgtggcggcaactggatttatgagtgggccccggatctttgtgaaggaaccttacttctgtggtgtgacataattggacaaactacctacagagatttaaagctctaaggtaaatataaaatttttaagtgtataatgtgttaaactactgattctaattgtttgtgtattttagattccaacctatggaactgatgaatgggagcagtggtggaatgcctttaatgaggaaaacctgttttgctcagaagaaatgccatctagtgatgatgaggctactgctgactctcaacattctactcctccaaaaaagaagagaaaggtagaagaccccaaggactttccttcagaattgctaagttttttgagtcatgctgtgtttagtaatagaactcttgcttgctttgctatttacaccacaaaggaaaaagctgcactgctatacaagaaaattatggaaaaatattctgtaacctttataagtaggcataacagttataatcataacatactgttttttcttactccacacaggcatagagtgtctgctattaataactatgctcaaaaattgtgtacctttagctttttaatttgtaaaggggttaataaggaatatttgatgtatagtgccttgactagagatcataatcagccataccacatttgtagaggttttacttgctttaaaaaacctcccacacctccccctgaacctgaaacataaaatgaatgcaattgttgttgttaacttgtttattgcagcttataatggttacaaataaagcaatagcatcacaaatttcacaaataaagcatttttttcactgcattctagttgtggtttgtccaaactcatcaatgtatcttatcatgtctggatcaactggataactcaagctaaccaaaatcatcccaaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcattttaaagaaattgtatttgttaaatatgtactacaaacttagtagt",
                          "name": "myseq2", "circular": False})
        j.quick_blastn()
        results = j.results.alignments
        print(results)

    def test_make_from_json_with_preloaded_id(self):
        """We expect the alignment results to contain the very same id and information in the query and subject
        alignments as the sequences that went into the JSONBlast alignment. Since this is a perfect alignment,
        we also expect the sequences to be spit back out."""
        j = JSONBlast([
            {"id": "ABCDEFG",
             "sequence": "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta",
             "name": "myseq",
             "circular": False}
        ], {"id": "1234",
            "sequence": "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta",
            "name": "myseq2",
            "circular": False})
        j.quick_blastn()
        results = j.results.alignments

        expected_query = {
            "bases": "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta".upper(),
            "name": "myseq2",
            "circular": False,
            "length": 93,
            "start": 1,
            "end": 93,
            "sequence_id": "1234",
            "strand": "plus"
        }

        expected_subject = {
            "bases": "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta".upper(),
            "name": "myseq",
            "circular": False,
            "length": 93,
            "start": 1,
            "end": 93,
            "strand": "plus",
            "sequence_id": "ABCDEFG"
        }

        assert results[0]['query'] == expected_query
        assert results[0]['subject'] == expected_subject

    def test_json_blast_with_orderdict(self):
        """JSONBlast should handle the subclasses of dict as inputs as well as dict"""
        from collections import OrderedDict
        j = JSONBlast([
            OrderedDict({"id": "ABCDEFG",
                         "sequence": "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta",
                         "name": "myseq",
                         "circular": False})
        ], OrderedDict({"id": "1234",
                        "sequence": "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta",
                        "name": "myseq2",
                        "circular": False}))
        j.quick_blastn()
        results = j.results.alignments

        expected_query = {
            "bases": "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta".upper(),
            "name": "myseq2",
            "circular": False,
            "length": 93,
            "start": 1,
            "end": 93,
            "sequence_id": "1234",
            "strand": "plus"
        }

        expected_subject = {
            "bases": "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta".upper(),
            "name": "myseq",
            "circular": False,
            "length": 93,
            "start": 1,
            "end": 93,
            "strand": "plus",
            "sequence_id": "ABCDEFG"
        }

        assert results[0]['query'] == expected_query
        assert results[0]['subject'] == expected_subject

    def test_json_blast_with_no_hits(self):
        """JSONBlast should handle the subclasses of dict as inputs as well as dict"""
        from collections import OrderedDict
        j = JSONBlast([
            OrderedDict({"id": "ABCDEFG",
                         "sequence": "aggagagagagaggagggagagaagaggagagagagaga",
                         "name": "myseq",
                         "circular": False})
        ], OrderedDict({"id": "1234",
                        "sequence": "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta",
                        "name": "myseq2",
                        "circular": False}))
        j.quick_blastn()
        results = j.results.alignments
        print(results)

    def test_json_blast_with_repeats(self):
        """We expect no hits from repeats"""
        from collections import OrderedDict
        j = JSONBlast([
            OrderedDict({"id": "ABCDEFG",
                         "sequence": "agtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtag",
                         "name": "myseq",
                         "circular": False})
        ], OrderedDict({"id": "1234",
                        "sequence": "agtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtag",
                        "name": "myseq2",
                        "circular": False}))
        j.quick_blastn()
        results = j.results.alignments
        assert results == ()

    def test_json_blast_with_preloaded_objects(self):
        class Sequence(object):

            def __init__(self, **kwargs):
                self.__dict__.update(**kwargs)

        query = Sequence(
            **{"id": "ABCDEFG",
             "bases": "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta",
             "name": "myseq",
             "circular": False}
        )

        from copy import copy
        subject = copy(query)

        j = JSONBlast([subject], query, preloaded=True)
        j.quick_blastn()
        results = j.results.alignments
        assert len(results) > 0

    def test_json_blast_with_objects_raises_pyblasterror(self):
        class Sequence(object):

            def __init__(self, **kwargs):
                self.__dict__.update(**kwargs)

        query = Sequence(
            **{"id": "ABCDEFG",
             "bases": "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta",
             "name": "myseq",
             "circular": False}
        )

        from copy import copy
        subject = copy(query)

        with pytest.raises(PyBlastException):
            j = JSONBlast([subject], query, preloaded=False)
            j.quick_blastn()
            assert len(j.results.alignments) == 1

class TestJSONBlastForExpectedSequences:

    @pytest.fixture
    def frag(self):
        return "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta"

    @pytest.fixture
    def seqs(self, frag):
        """Alignments that wraps around the query more than onces"""

        class Sequence(object):

            def __init__(self, **kwargs):
                self.__dict__.update(**kwargs)

        junk1 = "atgctatgctgatgctgctgtgctgatgctgatgtgtattgctgtatcgcgcgagttagc"
        junk2 = "g" * 30

        query = Sequence(
            **{"id": "myquery",
               "bases": frag,
               "name": "myseq",
               "circular": False}
        )

        subject = Sequence(
            **{"id": "mysubject",
               "bases": junk1 + frag + junk2,
               "name": "myseq",
               "circular": False}
        )

        return [subject], query

    @pytest.fixture
    def alignments(self, seqs):
        subjects, query = seqs
        j = JSONBlast(subjects, query, preloaded=True, span_origin=False)
        j.quick_blastn()
        results = j.results.get_perfect()
        return results.alignments

    @pytest.fixture
    def alignments_rc_subject(self, seqs):
        subjects, query = seqs
        subjects[0].bases = reverse_complement(subjects[0].bases)
        j = JSONBlast(subjects, query, preloaded=True, span_origin=False)
        j.quick_blastn()
        results = j.results.get_perfect()
        return results.alignments

    def test_expected_query_sequence(self, alignments, frag):
        assert len(alignments) == 1
        alignment = alignments[0]
        assert alignment['query']['strand'] == 'plus'
        assert alignment['query']['bases'].upper() == frag.upper()

    def test_expected_subject_sequence(self, alignments, frag):
        assert len(alignments) == 1
        alignment = alignments[0]
        assert alignment['subject']['strand'] == 'plus'
        assert alignment['subject']['bases'].upper() == frag.upper()

    def test_expected_length(self, alignments):
        for align in alignments:
            assert align['meta']['alignment_length'] == len(align['query']['bases'])
            assert align['meta']['alignment_length'] == len(align['subject']['bases'])
            assert align['query']['end'] - align['query']['start'] + 1 == len(align['query']['bases'])
            assert align['subject']['end'] - align['subject']['start'] + 1 == len(align['subject']['bases'])

    def test_expected_query_sequence_rc_subject(self, alignments_rc_subject, frag):
        assert len(alignments_rc_subject) == 1
        alignment = alignments_rc_subject[0]
        assert alignment['query']['strand'] == 'plus'
        assert alignment['query']['bases'].upper() == frag.upper()

    def test_expected_subject_sequence_rc_subject(self, alignments_rc_subject, frag):
        assert len(alignments_rc_subject) == 1
        alignment = alignments_rc_subject[0]
        assert alignment['subject']['strand'] == 'minus'
        assert alignment['subject']['bases'].upper() == frag.upper()

