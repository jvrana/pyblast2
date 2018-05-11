import pytest
from marshmallow import Schema, fields, ValidationError

from pyblast.models import QuerySchema, SubjectSchema, AlignmentMetaSchema


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
            'name': None,
            'filename': None,
            'circular': None,
            'query length': 10781
        },
        "subject": {
            "subject acc.": '756c05c4-f3f2-4e3b-a344-a4ed75827529',
            's. start': 1,
            's. end': 4219,
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
            expected = {"acc": 'Query_1',
                        'start': 1374,
                        'end': 5592,
                        'name': None,
                        'filename': None,
                        'circular': None,
                        'length': 10781
                        }
            assert loaded == expected

        def test_subject_load(self):
            schema = SubjectSchema()
            loaded = schema.load(TestSchema.preloaded_data['subject'])
            expected = {"acc": '756c05c4-f3f2-4e3b-a344-a4ed75827529',
                        'start': 1,
                        'end': 4219,
                        'name': None,
                        'filename': None,
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

    class TestLoadValidationErrors:

        def test_strand_validation(self):
            schema = SubjectSchema()
            data = TestSchema.preloaded_data['subject']
            data['subject strand'] = '+'
            with pytest.raises(ValidationError):
                schema.load(data)
