from uuid import uuid4

from marshmallow import Schema, fields, validates, ValidationError, validates_schema, pre_load


class SequenceSchema(Schema):
    size = fields.Integer()
    sequence = fields.String(required=True)
    circular = fields.Boolean(required=True)
    name = fields.String(required=True)
    id = fields.String(missing=lambda: str(uuid4()))
    description = fields.String(missing="")
    features = fields.Nested("FeatureSchema", many=True, missing=list())
    notes = fields.Dict(missing=dict())

    @pre_load
    def make_default_size(self, data):
        data = {k: v for k, v in data.items() if v is not None}
        if 'size' not in data:
            data['size'] = len(data['sequence'])
        return data

    @validates_schema
    def validate_size(self, data):
        len1 = len(data['sequence'])
        if not data['size'] == len1:
            raise ValidationError("size ({}) must be equal to sequence length ({})".format(data['size'], len1))


class FeatureSchema(Schema):
    name = fields.String(required=True)
    type = fields.String(required=True)
    id = fields.String(missing=lambda: str(uuid4()))
    start = fields.Int(required=True)
    end = fields.Int(required=True)
    strand = fields.Int(required=True)
    notes = fields.Dict(missing=dict())

    @validates('strand')
    def validate_strand(self, value):
        if value not in [1, -1]:
            raise ValidationError("Strand must by 1 [FORWARD] or -1 [REVERSE].")


class SequenceSchemaMixIn:
    """Sequence Schema methods common to QuerySchema and SubjectSchema"""
    # seq = fields.Method(serialize="get_sequence", deserialize="return_value", allow_none=True)
    name = fields.Method(serialize="get_name", allow_none=True)
    circular = fields.Method(serialize="get_circular", allow_none=True)

    def get_sequence(self, obj):
        """Searches the context for the sequence using the accession id"""
        if 'db' not in self.context:
            return None
        db = self.context['db']
        acc = obj['acc']
        if acc not in db:
            # TODO: raise error is acc is not in context?
            return None
        return db[obj['acc']]

    def get_name(self, obj):
        return self.get_sequence(obj)["name"]

    def get_circular(self, obj):
        return self.get_sequence(obj)["circular"]


class QuerySchema(Schema, SequenceSchemaMixIn):
    acc = fields.String(data_key="query acc.", required=True)
    # filename = fields.String(data_key="query filename")
    # circular = fields.Boolean(data_key="query circular")
    start = fields.Integer(data_key="q. start", required=True)
    end = fields.Integer(data_key="q. end", required=True)
    length = fields.Integer(data_key='query length', required=True)
    sequence = fields.String(data_key='query seq')
    # alignment = fields.Nested("AlignmentSchema")


class SubjectSchema(Schema, SequenceSchemaMixIn):
    acc = fields.String(data_key="subject acc.", required=True)
    # filename = fields.String(data_key="subject filename")
    # circular = fields.Boolean(data_key="subject circular")
    start = fields.Integer(data_key="s. start", required=True)
    end = fields.Integer(data_key="s. end", required=True)
    length = fields.Integer(data_key='subject length', required=True)
    sequence = fields.String(data_key='subject seq')
    strand = fields.String(data_key='subject strand')

    @validates("strand")
    def strand_validation(self, value):
        expected_values = ['plus', 'minus']
        if value not in expected_values:
            raise ValidationError(
                "subject_strand value must be one of the following: {}".format(','.join(expected_values)))
    # alignment = fields.Nested("AlignmentSchema")


class AlignmentMetaSchema(Schema):
    """Schema for alignment score metadata for an alignment"""
    score = fields.Float(required=True)
    evalue = fields.Float(required=True)
    bit_score = fields.Integer(data_key='bit score', required=True)

    identical = fields.Integer(required=True)
    gaps_open = fields.Integer(data_key="gap opens", required=True)
    gaps = fields.Integer(required=True)
    alignment_length = fields.Integer(data_key='alignment length', required=True)
    alignment = fields.Nested("AlignmentSchema")


class AlignmentSchema(Schema):
    """Schema for a BLAST alignment"""
    query = fields.Nested('QuerySchema', )
    subject = fields.Nested('SubjectSchema', )
    meta = fields.Nested(AlignmentMetaSchema, exclude=("alignment",))