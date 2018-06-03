import pytest
from marshmallow import ValidationError

from pyblast.schema import SequenceSchema, FeatureSchema


class TestSchema:
    preloaded_data = {
        "sequence": "gagagagagagagagggagagagcctatttaatat",
        "circular": True,
        "name": "pBbS8c-RFP",
        "description": "",
        "features": [
            {
                "name": "anonymous feature",
                "type": "misc_feature",
                "id": "5590c1978979df000a4f02c7",
                "start": 1,
                "end": 3,
                "strand": 1,
                "notes": {},
            },
            {
                "name": "coding region 1",
                "id": "5590c1978979df000a4f02c7",
                "type": "CDS",
                "start": 12,
                "end": 9,
                "strand": -1,
                "notes": {},
            }
        ],
    }


class TestLoad:
    def test_sequence_load(self):
        schema = SequenceSchema()
        loaded = schema.load(TestSchema.preloaded_data)
        loaded.pop('id', None)
        expected = {"size": 35, "notes": {}}
        expected.update(TestSchema.preloaded_data)
        assert loaded == expected

    def test_validation_error_wrong_size(self):
        schema = SequenceSchema()
        test_data = TestSchema.preloaded_data
        test_data['size'] = 50
        with pytest.raises(ValidationError):
            schema.load(test_data)

    def test_feature_load(self):
        schema = FeatureSchema()
        loaded = schema.load(TestSchema.preloaded_data["features"][0])
        expected = {
            "name": "anonymous feature",
            "type": "misc_feature",
            "id": "5590c1978979df000a4f02c7",
            "start": 1,
            "end": 3,
            "strand": 1,
            "notes": {},
        }
        assert loaded == expected

    def test_feature_load_no_id(self):
        schema = FeatureSchema()
        test_data = TestSchema.preloaded_data["features"][0]
        old_id = test_data['id']
        test_data.pop('id', None)
        loaded = schema.load(test_data)
        loaded_id = loaded['id']
        expected = {
            "name": "anonymous feature",
            "type": "misc_feature",
            "id": loaded_id,
            "start": 1,
            "end": 3,
            "strand": 1,
            "notes": {},
        }
        assert not loaded_id == old_id
        assert len(loaded_id) > 10
        assert loaded == expected

    def test_feature_load_error_wrong_strand(self):
        schema = FeatureSchema()
        test_data = TestSchema.preloaded_data["features"][0]
        test_data['strand'] = 0
        with pytest.raises(ValidationError):
            schema.load(test_data)


# class TestDump:
#
#     @pytest.fixture
#     def unmarshalled_data(self):
#         return {"acc": 'Query_1',
#                 'start': 1374,
#                 'end': 5592,
#                 'length': 10781
#                 }
#
#     @pytest.fixture
#     def seq(self):
#         return PySequence('', '', 'myseq', '', None, [], {}, {},
#                           filename='myfilename',
#                           circular=True)
#
#     def test_with_context(self, unmarshalled_data, seq):
#         context = {'db': {'Query_1': seq}}
#         schema_with_context = QuerySchema()
#         schema_with_context.context = context
#         dumped_with_context = schema_with_context.dump(unmarshalled_data)
#         assert dumped_with_context['filename'] == 'myfilename'
#         assert dumped_with_context['name'] == 'myseq'
#         assert dumped_with_context['circular']
#
#     def test_with_no_context(self, unmarshalled_data):
#         # no context
#         schema_no_context = QuerySchema()
#         dumped_no_context = schema_no_context.dump(unmarshalled_data)
#         assert dumped_no_context['filename'] is None
#         assert dumped_no_context['name'] is None
#         assert dumped_no_context['circular'] is None
#
#     def test_seq_missing_from_db(self, unmarshalled_data, seq):
#         context = {'db': {'Query_2': seq}}
#         schema_with_context = QuerySchema()
#         schema_with_context.context = context
#         dumped_with_context = schema_with_context.dump(unmarshalled_data)
#         assert dumped_with_context['filename'] is None
#         assert dumped_with_context['name'] is None
#         assert dumped_with_context['circular'] is None
#
# class TestLoadValidationErrors:
#
#     def test_strand_validation(self):
#         schema = SubjectSchema()
#         data = TestSchema.preloaded_data['subject']
#         data['subject strand'] = '+'
#         with pytest.raises(ValidationError):
#             schema.load(data)
