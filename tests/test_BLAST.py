import os
import json
import pytest

from pyblast import Blast, Aligner, JSONBlast
from pyblast.schema import SequenceSchema
from pyblast.utils.seq_parser import reverse_complement
from pyblast.exceptions import PyBlastException


def pytest_namespace(here):
    return {'b':
                Blast('db',
                      os.path.join(here, 'data/test_data/db.fsa'),
                      os.path.join(here, 'data/test_data/query.fsa'),
                      os.path.join(here, 'data/blast_results'),
                      os.path.join(here, 'data/blast_results/results.out')
                      )
            }


@pytest.fixture
def b(here):
    return Blast('db',
                 os.path.join(here, 'data/test_data/db.fsa'),
                 os.path.join(here, 'data/test_data/query.fsa'),
                 os.path.join(here, 'data/blast_results'),
                 os.path.join(here, 'data/blast_results/results.out')
                 )


def test_makedb(b):
    b.makedb()

def test_blasn(b):
    b.makedb()
    b.blastn()

def test_parse_results(b):
    b.makedb()
    b.blastn()
    b.parse_results()

    results = b.results.alignments
    res = results[0]
    assert 'query' in res
    assert 'subject' in res
    assert 'meta' in res

    query = res['query']
    assert 'sequence_id' in query
    assert 'start' in query
    assert 'end' in query
    assert 'length' in query
    assert query['name'] is None
    assert query['circular'] is None

    subject = res['subject']
    assert 'sequence_id' in subject
    assert 'start' in subject
    assert 'end' in subject
    assert 'length' in subject
    assert subject['name'] is None
    assert subject['circular'] is None
    assert subject['strand'] in ['plus', 'minus']

    meta = res['meta']
    assert 'score' in meta
    assert 'evalue' in meta
    assert 'bit_score' in meta
    assert 'identical' in meta
    assert 'gaps_open' in meta
    assert 'gaps' in meta
    assert 'alignment_length' in meta


def test_quick_blastn(b):
    b.quick_blastn()
    results = b.raw_results


class TestAligner:

    @pytest.fixture
    def aligner(self, here):
        template_dictionary = os.path.join(here, 'data/test_data/db.fsa')
        query_path = os.path.join(here, 'data/test_data/query.fsa')
        db_name = 'db'

        a = Aligner(db_name, template_dictionary, query_path)
        a.quick_blastn()
        return a

    def test_query(self, aligner):
        results = aligner.results

        for res in results.alignments:
            assert res['query']['circular'] is None
            assert res['query']['name'] is None

    def test_subject(self, aligner):
        results = aligner.results
        for res in results.alignments:
            assert res['subject']['circular'] is None
            assert res['subject']['name'] is None
            assert res['subject']['strand'] in ['plus', 'minus']

    def test_example(self):
        a = Aligner.use_test_data()
        a.quick_blastn()
    # def test_get_metadata():
    #     a = Aligner.use_test_data()
    #     a.quick_blastn()
    #     a.get_filename(a.input_sequences[0].id)
    #     a.get_is_circular(a.input_sequences[0].id)

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

        j = JSONBlast(subjects, query, task="blastn-short")
        j.quick_blastn()
        alignments = j.results.get_perfect().get_with_perfect_subjects().alignments
        print(len(j.results.alignments))
        print(len(alignments))

    def test_make_from_json(self):
        j = JSONBlast([{"sequence": "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta", "name": "myseq", "circular": False}], {"sequence": "tggaagggctaattcactcccaaagaagacaagatatccttgatctgtggatctaccacacacaaggctacttccctgattagcagaactacacaccagggccaggggtcagatatccactgacctttggatggtgctacaagctagtaccagttgagccagataaggtagaagaggccaataaaggagagaacaccagcttgttacaccctgtgagcctgcatgggatggatgacccggagagagaagtgttagagtggaggtttgacagccgcctagcatttcatcacgtggcccgagagctgcatccggagtacttcaagaactgctgatatcgagcttgctacaagggactttccgctggggactttccagggaggcgtggcctgggcgggactggggagtggcgagccctcagatcctgcatataagcagctgctttttgcctgtactgggtctctctggttagaccagatctgagcctgggagctctctggctaactagggaacccactgcttaagcctcaataaagcttgccttgagtgcttcaagtagtgtgtgcccgtctgttgtgtgactctggtaactagagatccctcagacccttttagtcagtgtggaaaatctctagcagtggcgcccgaacagggacttgaaagcgaaagggaaaccagaggagctctctcgacgcaggactcggcttgctgaagcgcgcacggcaagaggcgaggggcggcgactggtgagtacgccaaaaattttgactagcggaggctagaaggagagagatgggtgcgagagcgtcagtattaagcgggggagaattagatcgcgatgggaaaaaattcggttaaggccagggggaaagaaaaaatataaattaaaacatatagtatgggcaagcagggagctagaacgattcgcagttaatcctggcctgttagaaacatcagaaggctgtagacaaatactgggacagctacaaccatcccttcagacaggatcagaagaacttagatcattatataatacagtagcaaccctctattgtgtgcatcaaaggatagagataaaagacaccaaggaagctttagacaagatagaggaagagcaaaacaaaagtaagaccaccgcacagcaagcggccggccgctgatcttcagacctggaggaggagatatgagggacaattggagaagtgaattatataaatataaagtagtaaaaattgaaccattaggagtagcacccaccaaggcaaagagaagagtggtgcagagagaaaaaagagcagtgggaataggagctttgttccttgggttcttgggagcagcaggaagcactatgggcgcagcgtcaatgacgctgacggtacaggccagacaattattgtctggtatagtgcagcagcagaacaatttgctgagggctattgaggcgcaacagcatctgttgcaactcacagtctggggcatcaagcagctccaggcaagaatcctggctgtgaaagatacctaaaggatcaacagctcctggggatttggggttgctctggaaaactcatttgcaccactgctgtgccttggaatgctagttggagtaataaatctctggaacagatttggaatcacacgacctggatggagtgggacagagaaattaacaattacacaagcttaatacactccttaattgaagaatcgcaaaaccagcaagaaaagaatgaacaagaattattggaattagataaatgggcaagtttgtggaattggtttaacataacaaattggctgtggtatataaaattattcataatgatagtaggaggcttggtaggtttaagaatagtttttgctgtactttctatagtgaatagagttaggcagggatattcaccattatcgtttcagacccacctcccaaccccgaggggacccgacaggcccgaaggaatagaagaagaaggtggagagagagacagagacagatccattcgattagtgaacggatctcgacggtatcgccaaatggcagtattcatccacaattttaaaagaaaaggggggattggggggtacagtgcaggggaaagaatagtagacataatagcaacagacatacaaactaaagaattacaaaaacaaattacaaaaattcaaaattttcgggtttattacagggacagcagagatccagtttggatcgataagcttgatatcgaattcctgcagccccgataaaataaaagattttatttagtctccagaaaaaggggggaatgaaagaccccacctgtaggtttggcaagctagctgcagtaacgccattttgcaaggcatggaaaaataccaaaccaagaatagagaagttcagatcaagggcgggtacatgaaaatagctaacgttgggccaaacaggatatctgcggtgagcagtttcggccccggcccggggccaagaacagatggtcaccgcagtttcggccccggcccgaggccaagaacagatggtccccagatatggcccaaccctcagcagtttcttaagacccatcagatgtttccaggctcccccaaggacctgaaatgaccctgcgccttatttgaattaaccaatcagcctgcttctcgcttctgttcgcgcgcttctgcttcccgagctctataaaagagctcacaacccctcactcggcgcgccagtcctccgacagactgagtcgcccgggggggatctggagctctcgagaattctcacgcgtctgcaggatatcaagcttgcggtaccgcgggcccggccaccatggacaagaagtacagcatcggcctggccatcggcaccaactctgtgggctgggccgtgatcaccgacgagtacaaggtgcccagcaagaaattcaaggtgctgggcaacaccgaccggcacagcatcaagaagaacctgatcggcgccctgctgttcgacagcggagaaacagccgaggccacccggctgaagagaaccgccagaagaagatacaccagacggaagaaccggatctgctatctgcaagagatcttcagcaacgagatggccaaggtggacgacagcttcttccacagactggaagagtccttcctggtggaagaggataagaagcacgagcggcaccccatcttcggcaacatcgtggacgaggtggcctaccacgagaagtaccccaccatctaccacctgagaaagaaactggtggacagcaccgacaaggccgacctgcggctgatctatctggccctggcccacatgatcaagttccggggccacttcctgatcgagggcgacctgaaccccgacaacagcgacgtggacaagctgttcatccagctggtgcagacctacaaccagctgttcgaggaaaaccccatcaacgccagcggcgtggacgccaaggccatcctgtctgccagactgagcaagagcagacggctggaaaatctgatcgcccagctgcccggcgagaagaagaatggcctgttcggcaacctgattgccctgagcctgggcctgacccccaacttcaagagcaacttcgacctggccgaggatgccaaactgcagctgagcaaggacacctacgacgacgacctggacaacctgctggcccagatcggcgaccagtacgccgacctgtttctggccgccaagaacctgtccgacgccatcctgctgagcgacatcctgagagtgaacaccgagatcaccaaggcccccctgagcgcctctatgatcaagagatacgacgagcaccaccaggacctgaccctgctgaaagctctcgtgcggcagcagctgcctgagaagtacaaagagattttcttcgaccagagcaagaacggctacgccggctacatcgatggcggagccagccaggaagagttctacaagttcatcaagcccatcctggaaaagatggacggcaccgaggaactgctcgtgaagctgaacagagaggacctgctgcggaagcagcggaccttcgacaacggcagcatcccccaccagatccacctgggagagctgcacgccattctgcggcggcaggaagatttttacccattcctgaaggacaaccgggaaaagatcgagaagatcctgaccttccgcatcccctactacgtgggccctctggccaggggaaacagcagattcgcctggatgaccagaaagagcgaggaaaccatcaccccctggaacttcgaggaagtggtggacaagggcgccagcgcccagagcttcatcgagcggatgaccaacttcgataagaacctgcccaacgagaaggtgctgcccaagcacagcctgctgtacgagtacttcaccgtgtacaacgagctgaccaaagtgaaatacgtgaccgagggaatgagaaagcccgccttcctgagcggcgagcagaaaaaagccatcgtggacctgctgttcaagaccaaccggaaagtgaccgtgaagcagctgaaagaggactacttcaagaaaatcgagtgcttcgactccgtggaaatctccggcgtggaagatcggttcaacgcctccctgggcacataccacgatctgctgaaaattatcaaggacaaggacttcctggacaatgaggaaaacgaggacattctggaagatatcgtgctgaccctgacactgtttgaggacagagagatgatcgaggaacggctgaaaacctatgcccacctgttcgacgacaaagtgatgaagcagctgaagcggcggagatacaccggctggggcaggctgagccggaagctgatcaacggcatccgggacaagcagtccggcaagacaatcctggatttcctgaagtccgacggcttcgccaacagaaacttcatgcagctgatccacgacgacagcctgacctttaaagaggacatccagaaagcccaggtgtccggccagggcgatagcctgcacgagcacattgccaatctggccggcagccccgccattaagaagggcatcctgcagacagtgaaggtggtggacgagctcgtgaaagtgatgggccggcacaagcccgagaacatcgtgatcgaaatggccagagagaaccagaccacccagaagggacagaagaacagccgcgagagaatgaagcggatcgaagagggcatcaaagagctgggcagccagatcctgaaagaacaccccgtggaaaacacccagctgcagaacgagaagctgtacctgtactacctgcagaatgggcgggatatgtacgtggaccaggaactggacatcaaccggctgtccgactacgatgtggacgctatcgtgcctcagagctttctgaaggacgactccatcgataacaaagtgctgactcggagcgacaagaaccggggcaagagcgacaacgtgccctccgaagaggtcgtgaagaagatgaagaactactggcgccagctgctgaatgccaagctgattacccagaggaagttcgacaatctgaccaaggccgagagaggcggcctgagcgaactggataaggccggcttcatcaagagacagctggtggaaacccggcagatcacaaagcacgtggcacagatcctggactcccggatgaacactaagtacgacgagaacgacaaactgatccgggaagtgaaagtgatcaccctgaagtccaagctggtgtccgatttccggaaggatttccagttttacaaagtgcgcgagatcaacaactaccaccacgcccacgacgcctacctgaacgccgtcgtgggaaccgccctgatcaaaaagtaccctaagctggaaagcgagttcgtgtacggcgactacaaggtgtacgacgtgcggaagatgatcgccaagagcgagcaggaaatcggcaaggctaccgccaagtacttcttctacagcaacatcatgaactttttcaagaccgagattaccctggccaacggcgagatccggaagcggcctctgatcgagacaaacggcgaaacaggcgagatcgtgtgggataagggccgggactttgccaccgtgcggaaagtgctgtctatgccccaagtgaatatcgtgaaaaagaccgaggtgcagacaggcggcttcagcaaagagtctatcctgcccaagaggaacagcgacaagctgatcgccagaaagaaggactgggaccctaagaagtacggcggcttcgacagccccaccgtggcctattctgtgctggtggtggccaaagtggaaaagggcaagtccaagaaactgaagagtgtgaaagagctgctggggatcaccatcatggaaagaagcagcttcgagaagaatcccatcgactttctggaagccaagggctacaaagaagtgaaaaaggacctgatcatcaagctgcctaagtactccctgttcgagctggaaaacggccggaagagaatgctggcctctgccggcgaactgcagaagggaaacgaactggccctgccctccaaatatgtgaacttcctgtacctggccagccactatgagaagctgaagggctcccccgaggataatgagcagaaacagctgtttgtggaacagcacaaacactacctggacgagatcatcgagcagatcagcgagttctccaagagagtgatcctggccgacgctaatctggacaaggtgctgagcgcctacaacaagcacagagacaagcctatcagagagcaggccgagaatatcatccacctgtttaccctgaccaatctgggagcccctgccgccttcaagtactttgacaccaccatcgaccggaagaggtacaccagcaccaaagaggtgctggacgccaccctgatccaccagagcatcaccggcctgtacgagacacggatcgacctgtctcagctgggaggcgacgcctatccctatgacgtgcccgattatgccagcctgggcagcggctcccccaagaaaaaacgcaaggtggaagatcctaagaaaaagcggaaagtggacggcattggtagtgggagcaacggcagcagcggatccagcgagctgattaaggagaacatgcacatgaagctgtacatggagggcaccgtggacaaccatcacttcaagtgcacatccgagggcgaaggcaagccctacgagggcacccagaccatgagaatcaaggtggtcgagggcggccctctccccttcgccttcgacatcctggctactagcttcctctacggcagcaagaccttcatcaaccacacccagggcatccccgacttcttcaagcagtccttccctgagggcttcacatgggagagagtcaccacatacgaagacgggggcgtgctgaccgctacccaggacaccagcctccaggacggctgcctcatctacaacgtcaagatcagaggggtgaacttcacatccaacggccctgtgatgcagaagaaaacactcggctgggaggccttcaccgagacgctgtaccccgctgacggcggcctggaaggcagaaacgacatggccctgaagctcgtgggcgggagccatctgatcgcaaacatcaagaccacatatagatccaagaaacccgctaagaacctcaagatgcctggcgtctactatgtggactacagactggaaagaatcaaggaggccaacaacgagacctacgtcgagcagcacgaggtggcagtggccagatactgcgacctccctagcaaactggggcacaagcttaattagtaaggccgcgactctagagtcgacctgcaggcatgcaagcttgatatcaagcttatcgataatcaacctctggattacaaaatttgtgaaagattgactggtattcttaactatgttgctccttttacgctatgtggatacgctgctttaatgcctttgtatcatgctattgcttcccgtatggctttcattttctcctccttgtataaatcctggttgctgtctctttatgaggagttgtggcccgttgtcaggcaacgtggcgtggtgtgcactgtgtttgctgacgcaacccccactggttggggcattgccaccacctgtcagctcctttccgggactttcgctttccccctccctattgccacggcggaactcatcgccgcctgccttgcccgctgctggacaggggctcggctgttgggcactgacaattccgtggtgttgtcggggaaatcatcgtcctttccttggctgctcgcctgtgttgccacctggattctgcgcgggacgtccttctgctacgtcccttcggccctcaatccagcggaccttccttcccgcggcctgctgccggctctgcggcctcttccgcgtcttcgccttcgccctcagacgagtcggatctccctttgggccgcctccccgcatcgataccgtcgacctcgagggaattaattcgagctcggtacctttaagaccaatgacttacaaggcagctgtagatcttagccactttttaaaagaaaaggggggactggaagggctaattcactcccaacgaagacaagatctgctttttgcttgtactgggtctctctggttagaccagatctgagcctgggagctctctggctaactagggaacccactgcttaagcctcaataaagcttgccttgagtgcttcaagtagtgtgtgcccgtctgttgtgtgactctggtaactagagatccctcagacccttttagtcagtgtggaaaatctctagcagcatctagaattaattccgtgtattctatagtgtcacctaaatcgtatgtgtatgatacataaggttatgtattaattgtagccgcgttctaacgacaatatgtacaagcctaattgtgtagcatctggcttactgaagcagaccctatcatctctctcgtaaactgccgtcagagtcggtttggttggacgaaccttctgagtttctggtaacgccgtcccgcacccggaaatggtcagcgaaccaatcagcagggtcatcgctagccagatcctctacgccggacgcatcgtggccggcatcaccggcgccacaggtgcggttgctggcgcctatatcgccgacatcaccgatggggaagatcgggctcgccacttcgggctcatgagcgcttgtttcggcgtgggtatggtggcaggccccgtggccgggggactgttgggcgccatctccttgcatgcaccattccttgcggcggcggtgctcaacggcctcaacctactactgggctgcttcctaatgcaggagtcgcataagggagagcgtcgaatggtgcactctcagtacaatctgctctgatgccgcatagttaagccagccccgacacccgccaacacccgctgacgcgccctgacgggcttgtctgctcccggcatccgcttacagacaagctgtgaccgtctccgggagctgcatgtgtcagaggttttcaccgtcatcaccgaaacgcgcgagacgaaagggcctcgtgatacgcctatttttataggttaatgtcatgataataatggtttcttagacgtcaggtggcacttttcggggaaatgtgcgcggaacccctatttgtttatttttctaaatacattcaaatatgtatccgctcatgagacaataaccctgataaatgcttcaataatattgaaaaaggaagagtatgagtattcaacatttccgtgtcgcccttattcccttttttgcggcattttgccttcctgtttttgctcacccagaaacgctggtgaaagtaaaagatgctgaagatcagttgggtgcacgagtgggttacatcgaactggatctcaacagcggtaagatccttgagagttttcgccccgaagaacgttttccaatgatgagcacttttaaagttctgctatgtggcgcggtattatcccgtattgacgccgggcaagagcaactcggtcgccgcatacactattctcagaatgacttggttgagtactcaccagtcacagaaaagcatcttacggatggcatgacagtaagagaattatgcagtgctgccataaccatgagtgataacactgcggccaacttacttctgacaacgatcggaggaccgaaggagctaaccgcttttttgcacaacatgggggatcatgtaactcgccttgatcgttgggaaccggagctgaatgaagccataccaaacgacgagcgtgacaccacgatgcctgtagcaatggcaacaacgttgcgcaaactattaactggcgaactacttactctagcttcccggcaacaattaatagactggatggaggcggataaagttgcaggaccacttctgcgctcggcccttccggctggctggtttattgctgataaatctggagccggtgagcgtgggtctcgcggtatcattgcagcactggggccagatggtaagccctcccgtatcgtagttatctacacgacggggagtcaggcaactatggatgaacgaaatagacagatcgctgagataggtgcctcactgattaagcattggtaactgtcagaccaagtttactcatatatactttagattgatttaaaacttcatttttaatttaaaaggatctaggtgaagatcctttttgataatctcatgaccaaaatcccttaacgtgagttttcgttccactgagcgtcagaccccgtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaacaaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttccgaaggtaactggcttcagcagagcgcagataccaaatactgtctttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggactcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaactgagatacctacagcgtgagctatgagaaagcgccacgcttcccgaagggagaaaggcggacaggtatccggtaagcggcagggtcggaacaggagagcgcacgagggagcttccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccacctctgacttgagcgtcgatttttgtgatgctcgtcaggggggcggagcctatggaaaaacgccagcaacgcggcctttttacggttcctggccttttgctggccttttgctcacatgttctttcctgcgttatcccctgattctgtggataaccgtattaccgcctttgagtgagctgataccgctcgccgcagccgaacgaccgagcgcagcgagtcagtgagcgaggaagcggaagagcgcccaatacgcaaaccgcctctccccgcgcgttggccgattcattaatgcagctgtggaatgtgtgtcagttagggtgtggaaagtccccaggctccccagcaggcagaagtatgcaaagcatgcatctcaattagtcagcaaccaggtgtggaaagtccccaggctccccagcaggcagaagtatgcaaagcatgcatctcaattagtcagcaaccatagtcccgcccctaactccgcccatcccgcccctaactccgcccagttccgcccattctccgccccatggctgactaattttttttatttatgcagaggccgaggccgcctcggcctctgagctattccagaagtagtgaggaggcttttttggaggcctaggcttttgcaaaaagcttggacacaagacaggcttgcgagatatgtttgagaataccactttatcccgcgtcagggagaggcagtgcgtaaaaagacgcggactcatgtgaaatactggtttttagtgcgccagatctctataatctcgcgcaacctattttcccctcgaacactttttaagccgtagataaacaggctgggacacttcacatgagcgaaaaatacatcgtcacctgggacatgttgcagatccatgcacgtaaactcgcaagccgactgatgccttctgaacaatggaaaggcattattgccgtaagccgtggcggtctgtaccgggtgcgttactggcgcgtgaactgggtattcgtcatgtcgataccgtttgtatttccagctacgatcacgacaaccagcgcgagcttaaagtgctgaaacgcgcagaaggcgatggcgaaggcttcatcgttattgatgacctggtggataccggtggtactgcggttgcgattcgtgaaatgtatccaaaagcgcactttgtcaccatcttcgcaaaaccggctggtcgtccgctggttgatgactatgttgttgatatcccgcaagatacctggattgaacagccgtgggatatgggcgtcgtattcgtcccgccaatctccggtcgctaatcttttcaacgcctggcactgccgggcgttgttctttttaacttcaggcgggttacaatagtttccagtaagtattctggaggctgcatccatgacacaggcaaacctgagcgaaaccctgttcaaaccccgctttaaacatcctgaaacctcgacgctagtccgccgctttaatcacggcgcacaaccgcctgtgcagtcggcccttgatggtaaaaccatccctcactggtatcgcatgattaaccgtctgatgtggatctggcgcggcattgacccacgcgaaatcctcgacgtccaggcacgtattgtgatgagcgatgccgaacgtaccgacgatgatttatacgatacggtgattggctaccgtggcggcaactggatttatgagtgggccccggatctttgtgaaggaaccttacttctgtggtgtgacataattggacaaactacctacagagatttaaagctctaaggtaaatataaaatttttaagtgtataatgtgttaaactactgattctaattgtttgtgtattttagattccaacctatggaactgatgaatgggagcagtggtggaatgcctttaatgaggaaaacctgttttgctcagaagaaatgccatctagtgatgatgaggctactgctgactctcaacattctactcctccaaaaaagaagagaaaggtagaagaccccaaggactttccttcagaattgctaagttttttgagtcatgctgtgtttagtaatagaactcttgcttgctttgctatttacaccacaaaggaaaaagctgcactgctatacaagaaaattatggaaaaatattctgtaacctttataagtaggcataacagttataatcataacatactgttttttcttactccacacaggcatagagtgtctgctattaataactatgctcaaaaattgtgtacctttagctttttaatttgtaaaggggttaataaggaatatttgatgtatagtgccttgactagagatcataatcagccataccacatttgtagaggttttacttgctttaaaaaacctcccacacctccccctgaacctgaaacataaaatgaatgcaattgttgttgttaacttgtttattgcagcttataatggttacaaataaagcaatagcatcacaaatttcacaaataaagcatttttttcactgcattctagttgtggtttgtccaaactcatcaatgtatcttatcatgtctggatcaactggataactcaagctaaccaaaatcatcccaaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcattttaaagaaattgtatttgttaaatatgtactacaaacttagtagt", "name": "myseq2", "circular": False})
        j.quick_blastn()
        results = j.results.alignments
        print(results)

    def test_make_from_json_with_preloaded_id(self):
        """We expect the alignment results to contain the very same id and information in the query and subject
        alignments as the sequences that went into the JSONBlast alignment. Since this is a perfect alignment,
        we also expect the sequences to be spit back out."""
        j = JSONBlast(
            [
                {"id": "ABCDEFG",
                 "sequence": "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta",
                 "name": "myseq",
                 "circular": False}
            ],
            {"id": "1234",
             "sequence": "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta",
             "name": "myseq2",
             "circular": False}
        )
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
        j = JSONBlast(
            [
                OrderedDict({"id": "ABCDEFG",
                 "sequence": "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta",
                 "name": "myseq",
                 "circular": False})
            ],
            OrderedDict({"id": "1234",
             "sequence": "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta",
             "name": "myseq2",
             "circular": False})
        )
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
        j = JSONBlast(
            [
                OrderedDict({"id": "ABCDEFG",
                 "sequence": "aggagagagagaggagggagagaagaggagagagagaga",
                 "name": "myseq",
                 "circular": False})
            ],
            OrderedDict({"id": "1234",
             "sequence": "aaacttcccaccccataccctattaccactgccaattacctagtggtttcatttactctaaacctgtgattcctctgaattattttcatttta",
             "name": "myseq2",
             "circular": False})
        )
        j.quick_blastn()
        results = j.results.alignments
        print(results)

    def test_json_blast_with_repeats(self):
        """We expect no hits from repeats"""
        from collections import OrderedDict
        j = JSONBlast(
            [
                OrderedDict({"id": "ABCDEFG",
                 "sequence": "agtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtag",
                 "name": "myseq",
                 "circular": False})
            ],
            OrderedDict({"id": "1234",
             "sequence": "agtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtagagtag",
             "name": "myseq2",
             "circular": False})
        )
        j.quick_blastn()
        results = j.results.alignments
        print(results)

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

        r = SequenceSchema().dump(query)

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

        r = SequenceSchema().dump(query)

        from copy import copy
        subject = copy(query)

        from collections import OrderedDict
        with pytest.raises(PyBlastException):
            j = JSONBlast([subject], query, preloaded=False)