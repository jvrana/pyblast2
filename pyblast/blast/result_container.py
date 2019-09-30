from pyblast.utils import Span


class ResultEntry(object):

    def __init__(self, start, raw_end, strand, circular, record_key, record_id, record, origin_key, origin_record_id, origin_record):
        self.start = start
        self.raw_end = raw_end
        self.strand = strand

        self.record_key = record_key
        self.record_id = record_id
        self.record = record

        self.name = record.name

        self.origin_key = origin_key
        self.origin_record_id = origin_record_id
        self.origin_record = origin_record

        self.origin_sequence_length = len(origin_record)
        self.circular = circular

        span = self.to_span()
        if
        self.start = span.a
        self.end = span.b - 1
        self.raw_end = span.c - 1

        if output_index is None:
            output_index = input_index
        data['raw_end'] = data['end']
        s, e = data["start"], data["end"]
        reverse = data["strand"] != 1
        if reverse:
            s, e = e, s
        if e - s < 0:
            raise ValueError("End cannot be less than start")
        elif data['circular']:
            span = self.parse_result_to_span(data=data,
                                             inclusive=inclusive,
                                             input_index=input_index,
                                             output_index=output_index)
            data["start"] = span.a
            data["end"] = span.b - 1
            data["raw_end"] = span.c - 1
            if reverse:
                data["start"], data["end"], data['raw_end'] = data["end"], data["start"], data['start']

    def to_span(self, is_inclusive=True, input_index=1, output_index=None):
        s, e, l, c = self.start, self.raw_end, self.origin_sequence_length, self.circular
        if self.strand == -1:
            s, e = e, s
        if is_inclusive:
            e += 1
        span = Span(s, e, l, cyclic=c, allow_wrap=True, index=input_index)
        if input_index != output_index:
            span = span.reindex(output_index)
        return span

    def json(self):
        return {
            'start': self.start,
            'end': self.end,
            'raw_end': self.raw_end,
            'name': self.name,
            'record': record,
            'record_key'
        }