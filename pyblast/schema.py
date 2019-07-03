"""schema.py"""

from uuid import uuid4

from marshmallow import Schema, fields, validates, validates_schema, pre_load, EXCLUDE
from marshmallow import ValidationError

from pyblast.utils.dna_bases import rc_dict


class SequenceSchema(Schema):
    """
    Sequence JSON schema.

    .. code-block::

        size = fields.Method("get_size")
        bases = fields.String(required=True, data_key='sequence')
        circular = fields.Boolean(required=True)
        name = fields.String(required=True)
        id = fields.String(default=lambda: str(uuid4()))
        description = fields.String(default="", allow_none=True)
        features = fields.Nested("FeatureSchema", many=True, default=list())
    """

    class Meta:
        unknown = EXCLUDE

    size = fields.Method("get_size")
    bases = fields.String(required=True)
    circular = fields.Boolean(required=True)
    name = fields.String(required=True)
    id = fields.String(default=lambda: str(uuid4()))
    description = fields.String(default="", allow_none=True)
    features = fields.Nested("FeatureSchema", many=True, default=list())

    def get_size(self, obj):
        if hasattr(obj, "bases"):
            return len(obj.bases)
        else:
            return len(obj["bases"])

    def clean_data(self, data):
        if issubclass(data.__class__, dict):
            return {k: v for k, v in data.items() if v is not None}
        return data

    # @validates_schema
    # def validate_size(self, data):
    #     len1 = len(data['sequence'])
    #     if not data['size'] == len1:
    #         raise ValidationError("size ({}) must be equal to sequence length ({})".format(data['size'], len1))

    @validates("bases")
    def validate_sequence(self, bases):
        unknown_bases = []
        for b in bases:
            if b not in rc_dict and b not in unknown_bases:
                unknown_bases.append(b)
        if len(unknown_bases) > 0:
            raise ValidationError(
                "Found unknown base(s): {}".format(", ".join(unknown_bases))
            )
        if bases.strip() == "":
            raise ValidationError("Sequence cannot be empty")


class FeatureSchema(Schema):
    """
    .. code-block::

        name = fields.String(required=True)
        type = fields.String(required=True)
        id = fields.String(missing=lambda: str(uuid4()))
        start = fields.Int(required=True)
        end = fields.Int(required=True)
        strand = fields.Int(required=True)

    """

    name = fields.String(required=True)
    type = fields.String(required=True)
    id = fields.String(missing=lambda: str(uuid4()))
    start = fields.Int(required=True)
    end = fields.Int(required=True)
    strand = fields.Int(required=True)

    @validates("strand")
    def validate_strand(self, value):
        if value not in [1, -1]:
            raise ValidationError("Strand must by 1 [FORWARD] or -1 [REVERSE].")

    @pre_load
    def clean_data(self, data):
        return {k: v for k, v in data.items() if v is not None}


class SequenceSchemaMixIn:
    """Sequence Schema methods common to QuerySchema and SubjectSchema.

    .. code-block::

        name = fields.Method(serialize="get_name", deserialize="get_val", allow_none=True, )
        circular = fields.Method(serialize="get_circular", deserialize="get_val", allow_none=True, )

    """

    # seq = fields.Method(serialize="get_sequence", deserialize="return_value", allow_none=True)
    name = fields.Method(serialize="get_name", deserialize="get_val", allow_none=True)
    circular = fields.Method(
        serialize="get_circular", deserialize="get_val", allow_none=True
    )

    def get_val(self, x):
        return x

    def get_sequence(self, obj):
        """Searches the context for the sequence using the accession id"""
        # if 'db' not in self.context:
        #     return {
        #         "name": None,
        #         "circular": None
        #     }
        db = self.context["db"]
        seq_id = obj["sequence_id"]
        if seq_id not in db:
            raise ValidationError(
                'No sequence was found with accession ID "{}". Sequence must be found '
                "in schema context with key 'db' or must be run with no context".format(
                    seq_id
                ),
                "acc",
            )
        return db[obj["sequence_id"]]

    def get_name(self, obj):
        if "db" not in self.context:
            return None
        seq = self.get_sequence(obj)
        return seq["name"]

    def get_circular(self, obj):
        if "db" not in self.context:
            return None
        seq = self.get_sequence(obj)
        return seq["circular"]

    def load_sequence_id(self, value):
        return str(value)

    def get_sequence_id(self, obj):
        return obj["sequence_id"]


class QuerySchema(Schema, SequenceSchemaMixIn):
    """
    .. code-block::

        sequence_id = fields.Method("get_sequence_id", deserialize="load_sequence_id", data_key='query acc.', required=True)
        start = fields.Integer(data_key="q. start", required=True)
        end = fields.Integer(data_key="q. end", required=True)
        length = fields.Integer(data_key='query length', required=True)
        bases = fields.String(data_key='query seq', required=True)
        strand = fields.String(missing="plus", default="plus")
    """

    # acc = fields.String(data_key="query acc.", required=True)
    class Meta:
        unknown = EXCLUDE

    sequence_id = fields.Method(
        "get_sequence_id",
        deserialize="load_sequence_id",
        data_key="query acc.",
        required=True,
    )
    start = fields.Integer(data_key="q. start", required=True)
    end = fields.Integer(data_key="q. end", required=True)
    length = fields.Integer(data_key="query length", required=True)
    bases = fields.String(data_key="query seq", required=True)
    strand = fields.String(missing="plus", default="plus")


class SubjectSchema(Schema, SequenceSchemaMixIn):
    """
    .. code-block::

        sequence_id = fields.Method("get_sequence_id", deserialize="load_sequence_id", data_key='subject acc.', required=True)
        # acc = fields.String(data_key="subject acc.", required=True)
        start = fields.Integer(data_key="s. start", required=True)
        end = fields.Integer(data_key="s. end", required=True)
        length = fields.Integer(data_key='subject length', required=True)
        bases = fields.String(data_key='subject seq', required=True)
        strand = fields.String(data_key='subject strand', required=True)
    """

    class Meta:
        unknown = EXCLUDE

    sequence_id = fields.Method(
        "get_sequence_id",
        deserialize="load_sequence_id",
        data_key="subject acc.",
        required=True,
    )
    # acc = fields.String(data_key="subject acc.", required=True)
    start = fields.Integer(data_key="s. start", required=True)
    end = fields.Integer(data_key="s. end", required=True)
    length = fields.Integer(data_key="subject length", required=True)
    bases = fields.String(data_key="subject seq", required=True)
    strand = fields.String(data_key="subject strand", required=True)

    @validates("strand")
    def strand_validation(self, value):
        expected_values = ["plus", "minus"]
        if value not in expected_values:
            raise ValidationError(
                "subject_strand value must be one of the following: {}".format(
                    ",".join(expected_values)
                ),
                "strand",
            )


class AlignmentMetaSchema(Schema):
    """
    Schema for alignment score metadata for an alignment

    .. code-block::

        score = fields.Float(required=True)
        evalue = fields.Float(required=True)
        bit_score = fields.Integer(data_key='bit score', required=True)

        identical = fields.Integer(required=True)
        gaps_open = fields.Integer(data_key="gap opens", required=True)
        gaps = fields.Integer(required=True)
        alignment_length = fields.Integer(data_key='alignment length', required=True)
        alignment = fields.Nested("AlignmentSchema")
        span_origin = fields.Boolean(default=False)

    """

    class Meta:
        unknown = EXCLUDE

    score = fields.Float(required=True)
    evalue = fields.Float(required=True)
    bit_score = fields.Integer(data_key="bit score", required=True)

    identical = fields.Integer(required=True)
    gaps_open = fields.Integer(data_key="gap opens", required=True)
    gaps = fields.Integer(required=True)
    alignment_length = fields.Integer(data_key="alignment length", required=True)
    alignment = fields.Nested("AlignmentSchema")
    span_origin = fields.Boolean(default=False)


class AlignmentSchema(Schema):
    """Schema for a BLAST alignment

    .. code-block::
        query = fields.Nested('QuerySchema', )
        subject = fields.Nested('SubjectSchema', )
        meta = fields.Nested(AlignmentMetaSchema, exclude=("alignment",))
    """

    query = fields.Nested("QuerySchema")
    subject = fields.Nested("SubjectSchema")
    meta = fields.Nested(AlignmentMetaSchema, exclude=("alignment",))

    @validates_schema
    def validate_alignment_length(self, data):
        query_bases = data["query"]["bases"]
        subject_bases = data["subject"]["bases"]


# class BlastnParametersSchema(Schema):
#
#     word_size = fields.Integer(description="Length of initial exact match.")
#     gapopen = fields.Integer(description="DefaultCost to open a gap. See appendix “BLASTN reward/penalty values”.")
