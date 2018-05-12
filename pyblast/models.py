from marshmallow import Schema, fields, validates, ValidationError


class SequenceSchema:
    """Sequence Schema methods common to QuerySchema and SubjectSchema"""
    # seq = fields.Method(serialize="get_sequence", deserialize="return_value", allow_none=True)
    filename = fields.Method(serialize="get_filename", deserialize="return_value", allow_none=True)
    name = fields.Method(serialize="get_name", deserialize="return_value", allow_none=True)
    circular = fields.Method(serialize="get_circular", deserialize="return_value", allow_none=True)

    #
    # name = fields.Function(lambda seq, context: context['db'][seq['acc']])
    # name = fields.Function(lambda seq, context: context['db'][seq['acc']])

    def return_value(self, value):
        return value

    def get_sequence(self, obj):
        if 'db' not in self.context:
            return None
        db = self.context['db']
        acc = obj['acc']
        if acc not in db:
            return None
        return db[obj['acc']]

    def get_seq_attribute(self, obj, attr):
        seq = self.get_sequence(obj)
        if seq is not None:
            return getattr(seq, attr)

    def get_name(self, obj):
        return self.get_seq_attribute(obj, "name")

    def get_filename(self, obj):
        return self.get_seq_attribute(obj, "filename")

    def get_circular(self, obj):
        return self.get_seq_attribute(obj, "circular")


class QuerySchema(Schema, SequenceSchema):
    acc = fields.String(data_key="query acc.", required=True)
    # filename = fields.String(data_key="query filename")
    # circular = fields.Boolean(data_key="query circular")
    start = fields.Integer(data_key="q. start", required=True)
    end = fields.Integer(data_key="q. end", required=True)
    length = fields.Integer(data_key='query length', required=True)
    sequence = fields.String(data_key='query seq')
    # alignment = fields.Nested("AlignmentSchema")


class SubjectSchema(Schema, SequenceSchema):
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
    meta = fields.Nested("AlignmentMetaSchema", exclude=("alignment",))
