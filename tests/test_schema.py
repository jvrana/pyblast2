import pytest
from marshmallow import Schema, fields, ValidationError

from pyblast.schema import QuerySchema, SubjectSchema, AlignmentMetaSchema

def test_marshmallow_is_functional():
    class UserSchema(Schema):
        name = fields.String()
        email = fields.Email(data_key='emailAddress')
        d = fields.Dict(default={})
        method = fields.Method(deserialize='get_my_data')

        def get_my_data(self, obj):
            return 12345

    s = UserSchema()

    data = {
        'name': 'Mike',
        'email': 'foo@bar.com'
    }
    result = s.dump(data)
    assert result == {'name': u'Mike', 'emailAddress': 'foo@bar.com', 'd': {}}

    data = {
        'name': 'Mike',
        'emailAddress': 'foo@bar.com'
    }
    result = s.load(data)
    assert result == {'name': u'Mike', 'email': 'foo@bar.com'}


class TestSchema:
    preloaded_data = {
        "query": {
            "query acc.": 'Query_1',
            'q. start': 1374,
            'q. end': 5592,
            'query seq': 'AGTAT',
            'name': None,
            'filename': None,
            'circular': None,
            'query length': 10781
        },
        "subject": {
            "subject acc.": '756c05c4-f3f2-4e3b-a344-a4ed75827529',
            's. start': 1,
            's. end': 4219,
            'subject seq': 'aggagagag',
            'name': None,
            'filename': None,
            'circular': None,
            'subject length': 7883,
            'subject strand': 'plus',
        },
        "meta": {'score': 4219,
                 'evalue': 0,
                 'bit score':
                     7792,
                 'alignment length': 4219,
                 'identical': 4219,
                 'gap opens': 0,
                 'gaps': 0
                 }
    }

    class TestLoad:
        def test_query_load(self):
            schema = QuerySchema()
            loaded = schema.load(TestSchema.preloaded_data['query'])
            expected = {'sequence_id': 'Query_1',
                        'bases': 'AGTAT',
                        'start': 1374,
                        'end': 5592,
                        'name': None,
                        'circular': None,
                        'length': 10781,
                        "strand": "plus"
                        }
            assert loaded == expected

        def test_subject_load(self):
            schema = SubjectSchema()
            loaded = schema.load(TestSchema.preloaded_data['subject'])
            expected = {'sequence_id': '756c05c4-f3f2-4e3b-a344-a4ed75827529',
                        'bases': 'aggagagag',
                        'start': 1,
                        'end': 4219,
                        'name': None,
                        'circular': None,
                        'length': 7883,
                        'strand': 'plus', }
            assert loaded == expected

        def test_meta_load(self):
            schema = AlignmentMetaSchema()
            loaded = schema.load(TestSchema.preloaded_data['meta'])
            expected = {'score': 4219,
                        'evalue': 0,
                        'bit_score': 7792,
                        'alignment_length': 4219,
                        'identical': 4219,
                        'gaps_open': 0,
                        'gaps': 0
                        }
            assert loaded == expected

    class TestDump:

        @pytest.fixture
        def unmarshalled_data(self):
            return {'sequence_id': 'Query_1',
                    'start': 1374,
                    'end': 5592,
                    'length': 10781
                    }

        @pytest.fixture
        def seq(self):
            return {
                "name": "myseq",
                "circular": True,
            }

        def test_with_context(self, unmarshalled_data, seq):
            context = {'db': {'Query_1': seq}}
            schema_with_context = QuerySchema()
            schema_with_context.context = context
            dumped_with_context = schema_with_context.dump(unmarshalled_data)
            assert dumped_with_context['name'] == 'myseq'
            assert dumped_with_context['circular']

        def test_with_no_context(self, unmarshalled_data):
            # no context
            schema_no_context = QuerySchema()
            dumped_no_context = schema_no_context.dump(unmarshalled_data)
            assert dumped_no_context['name'] is None
            assert dumped_no_context['circular'] is None

        def test_seq_missing_from_db(self, unmarshalled_data, seq):
            context = {'db': {'Query_2': seq}}
            schema_with_context = QuerySchema()
            schema_with_context.context = context
            with pytest.raises(ValidationError):
                dumped_with_context = schema_with_context.dump(unmarshalled_data)

    class TestLoadValidationErrors:

        def test_strand_validation(self):
            schema = SubjectSchema()
            data = TestSchema.preloaded_data['subject']
            data['subject strand'] = '+'
            with pytest.raises(ValidationError):
                schema.load(data)
