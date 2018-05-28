"""Data parser for BLAST"""
import re

from pyblast.schema import QuerySchema, SubjectSchema, AlignmentSchema, AlignmentMetaSchema, SequenceSchema


def str_to_f_to_i(v):
    try:
        v = float(v)
    except ValueError:
        pass
    try:
        v = int(v)
    except ValueError:
        pass
    return v


def _extract_metadata(r, delim):
    """Extracts information from the raw text file BLAST produces"""
    g = re.search(
        '#\s*(?P<blast_ver>.+)\n' +
        '# Query:\s*(?P<query>.*)\n' +
        '# Database:\s*(?P<database>.+)\n' +
        '# Fields:\s*(?P<fields>.+)',
        r)
    metadata = g.groupdict()
    fields_array = re.split('\s*{}\s*'.format(delim), metadata['fields'])
    metadata['fields'] = fields_array
    return metadata


def _get_alignment_rows(r):
    """Split text into alignment rows"""
    return re.findall('\n([^#].*)', r)


def _validate_matches(raw_matches, fields):
    """Create a dictionary from the fields and rows"""
    match_dicts = []
    for m in raw_matches:
        values = [str_to_f_to_i(v) for v in m.split('\t')]
        match_dicts.append(dict(list(zip(fields, values))))
    return match_dicts


def _cleanup_raw_results(raw_text, delim):
    """
    Converts raw BLAST text into a flatten dictionary

    :param raw_text: raw text from BLAST results
    :type raw_text: basestring
    :param delim: delimiter for parsing
    :type delim: basestring
    :return: flattened dictionary
    :rtype: dict
    """

    # print(results)
    if raw_text.strip() == '':
        return {}
    meta = _extract_metadata(raw_text, delim)
    fields = tuple(meta['fields'])
    alignment_rows = _get_alignment_rows(raw_text)
    match_dicts = _validate_matches(alignment_rows, fields)
    return match_dicts


def _serialize_data(data, context):
    """
    Serializes a flat dictionary of data into Query, Subject, AlginmentMeta, and Alignment models

    :param data: flatten dictionary
    :type data: dict
    :return: serialized data
    :rtype: dict
    """

    query_schema = QuerySchema(many=True)
    subject_schema = SubjectSchema(many=True)
    meta_schema = AlignmentMetaSchema(many=True)
    alignment_schema = AlignmentSchema(many=True)
    alignment_schema.context = context

    # load marshall alignment
    queries = query_schema.load(data)
    subjects = subject_schema.load(data)
    metas = meta_schema.load(data)

    serialized = alignment_schema.dump([
        {"query": q, "subject": s, "meta": m} for
        q, s, m in zip(queries, subjects, metas)
    ])
    return serialized


def parse_results(data, delim, context=None):
    """Convert raw BLAST results into Alignment model"""
    # first clean up the data and get the flattened dictionary
    cleaned = _cleanup_raw_results(data, delim)

    # force data to deserialized according to the schemas
    # if context is None:
    #     context = {'db': {}}
    serialized = _serialize_data(cleaned, context)

    # force data back into models
    alignment_schema = AlignmentSchema(many=True)
    alignment_schema.context = context
    deserialized = alignment_schema.load(serialized)

    return deserialized

def parse_sequence_jsons(data):
    many = type(data) is list
    schema = SequenceSchema(many=many)
    return schema.load(data)
