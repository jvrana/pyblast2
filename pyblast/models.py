from marshmallow import Schema, fields, validates, ValidationError


# class Alignment(object):
#     def __init__(self, query_acc, subject_acc, query_seq, subject_seq,
#                  score, evalue, bit_score,
#                  identical, gaps_open, gaps, alignment_length,
#                  query_length, subject_length, q_start, q_end,
#                  s_start, s_end, subject_strand):
#         self.query_acc = query_acc
#         self.subject_acc = subject_acc
#         self.query_seq = query_seq
#         self.subject_seq = subject_seq
#         self.score = score
#         self.evalue = evalue
#         self.bit_score = bit_score
#         self.identical = identical
#         self.gaps_open = gaps_open
#         self.gaps = gaps
#         self.alignment_length = alignment_length
#         self.query_length = query_length
#         self.subject_length = subject_length
#         self.q_start = q_start
#         self.q_end = q_end
#         self.s_start = s_start
#         self.s_end = s_end
#         self.subject_strand = subject_strand

class Sequence:

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



class QuerySchema(Schema, Sequence):

    acc = fields.String(data_key="query acc.", required=True)
    # filename = fields.String(data_key="query filename")
    # circular = fields.Boolean(data_key="query circular")
    start = fields.Integer(data_key="q. start", required=True)
    end = fields.Integer(data_key="q. end", required=True)
    length = fields.Integer(data_key='query length', required=True)
    sequence = fields.String(data_key='query seq')
    # alignment = fields.Nested("AlignmentSchema")


class SubjectSchema(Schema, Sequence):
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
            raise ValidationError("subject_strand value must be one of the following: {}".format(','.join(expected_values)))
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
    subject = fields.Nested('SubjectSchema',)
    meta = fields.Nested("AlignmentMetaSchema", exclude=("alignment",))



class AlignmentSchemaOld(Schema):

    # query and subject
    query_acc = fields.String(data_key='query acc.', required=True)
    subject_acc = fields.String(data_key='subject acc.', required=True)
    query_seq = fields.String(data_key='query seq', required=True)
    subject_seq = fields.String(data_key='subject seq', required=True)

    # alignment scores
    score = fields.Float(required=True)
    evalue = fields.Float(required=True)
    bit_score = fields.Integer(data_key='bit score', required=True)

    identical = fields.Integer(required=True)
    gaps_open = fields.Integer(data_key="gap opens", required=True)
    gaps = fields.Integer(required=True)
    alignment_length = fields.Integer(data_key='alignment length', required=True)

    query_length = fields.Integer(data_key='query length', required=True)
    subject_length = fields.Integer(data_key='subject length', required=True)
    q_start = fields.Integer(data_key='q. start', required=True)
    q_end = fields.Integer(data_key='q. end', required=True)
    s_start = fields.Integer(data_key='s. start', required=True)
    s_end = fields.Integer(data_key='s. end', required=True)
    subject_strand = fields.String(data_key='subject strand', required=True)



    # subject = fields.Method(serialize="get_subject", dump_only=True)
    # subject = fields.Method(serialize="get_subject", dump_only=True)
    # seq_db = fields.Dict(default={}, load_only=True)
    #
    # def get_subject(self, obj):
    #     seq_db = obj['seq_db']
    #     acc = obj['subject_acc']
    #     if acc in seq_db:
    #         return seq_db[acc]
    #
    # def get_query(self, obj):
    #     seq_db = obj['seq_db']
    #     acc = obj['query_acc']
    #     if acc in seq_db:
    #         return seq_db[acc]
    #     return 'none'
    #
    #
    @validates("subject_strand")
    def strand_validation(self, value):
        expected_values = ['plus', 'minus']
        if value not in expected_values:
            raise ValidationError("subject_strand value must be one of the following: {}".format(','.join(expected_values)))


# class AlignmentMetaData(Schema):
#
#     # subject_name = fields.String(required=True)
#     # subject_filename = fields.String(required=True)
#     # subject_circular = fields.Boolean(required=True)
#     #
#     # query_name = fields.String(required=True)
#     # query_filename = fields.String(required=True)
#     # query_circular = fields.Boolean(required=True)
#
#     subject = fields.Method(serialize="get_subject", dump_only=True)
#
#     seq_db = fields.Dict(load_only=True, default={})



