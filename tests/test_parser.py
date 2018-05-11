from pyblast.parser import _serialize_data
from pyblast.pysequence import PySequence


class TestSerializer:
    raw_data = {'query acc.': 'Query_1', 'subject acc.': '756c05c4-f3f2-4e3b-a344-a4ed75827529', 'score': 4219,
                'evalue': 0, 'bit score': 7792, 'alignment length': 4219, 'identical': 4219, 'gap opens': 0, 'gaps': 0,
                'query length': 10781, 'q. start': 1374, 'q. end': 5592, 'subject length': 7883, 's. start': 1,
                's. end': 4219, 'subject strand': 'plus',
                'query seq': 'AAAAA',
                'subject seq': 'TCGCG', }

    context = {'db': {
        '756c05c4-f3f2-4e3b-a344-a4ed75827529': PySequence('', '', 'myseq', '', None, [], {}, {}, filename='myfilename',
                                                           circular=True)
    }}

    serialized = _serialize_data([raw_data], context)
    expected = {
        'query': {
            'query acc.': 'Query_1',
            'query length': 10781,
            'q. start': 1374, 'q. end': 5592,
            'circular': None,
            'filename': None,
            'name': None,
            'query seq': 'AAAAA',
        },
        'subject': {
            'subject acc.': '756c05c4-f3f2-4e3b-a344-a4ed75827529',
            'subject length': 7883,
            's. start': 1,
            's. end': 4219,
            'subject strand': 'plus',
            'circular': True,
            'filename': 'myfilename',
            'name': 'myseq',
            'subject seq': 'TCGCG'
        },
        'meta': {
            'score': 4219, 'evalue': 0,
            'bit score': 7792, 'alignment length': 4219,
            'identical': 4219, 'gap opens': 0, 'gaps': 0
        }
    }

    def test_serializes_subject(self):
        assert TestSerializer.expected['subject'] == TestSerializer.serialized[0]['subject']

    def test_serializes_query(self):
        assert TestSerializer.expected['query'] == TestSerializer.serialized[0]['query']

    def test_serializes_meta(self):
        assert TestSerializer.expected['meta'] == TestSerializer.serialized[0]['meta']
