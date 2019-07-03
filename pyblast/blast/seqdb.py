from .constants import Constants as C
from uuid import uuid4
from pyblast.exceptions import SeqRecordValidationError


class SeqRecordDB(object):
    def __init__(self):
        self.records = {}
        self.mapping = {}

    def incr(self):
        return str(uuid4())

    @classmethod
    def validate_records(cls, records):
        for r in records:
            cls.validate_circular(r)

    @staticmethod
    def validate_circular(r):
        if C.LINEAR not in r.annotations and C.CIRCULAR not in r.annotations:
            raise SeqRecordValidationError(
                "SeqRecord {} is missing a '{linear}' or '{circular}' annotation. This must be provided.".format(
                    r, circular=C.CIRCULAR, linear=C.LINEAR
                )
            )
        if C.LINEAR in r.annotations and C.CIRCULAR in r.annotations:
            if r.annotations[C.LINEAR] is r.annotations[C.CIRCULAR]:
                raise SeqRecordValidationError(
                    "SeqRecord {} has conflicting topology definitions for '{linear}' and '{circular}'".format(
                        r, circular=C.CIRCULAR, linear=C.LINEAR
                    )
                )

    @staticmethod
    def is_circular(r):
        return (
            r.annotations.get(C.CIRCULAR, False) is True
            or r.annotations.get(C.LINEAR, True) is False
        )

    def key(self, record):
        k = self.mapping.get(id(record), None)
        return k

    def post_transform_hook(self, key, *args, **kwargs):
        record = self.get(key)
        record.annotations[C.OLD_KEY] = record.id
        record.id = key
        return record

    def transform(self, key, transform, transform_label):
        record = self.get(key)
        new_record = transform(record)
        parent_key = self.add_one(record)
        new_key = self.add_one(new_record)
        new_record.annotations[C.PARENT] = parent_key
        new_record.annotations[C.TRANSFORMATION] = transform_label
        self.post_transform_hook(new_key)
        return new_key

    def add_with_transformation(self, record, transform, transform_label):
        key = self.add_one(record)
        return self.transform(key, transform, transform_label)

    def add_many_with_transformations(self, records, transform, transform_label):
        return [
            self.add_with_transformation(r, transform, transform_label) for r in records
        ]

    def get_many(self, keys, default=C.NULL):
        records = []
        for k in keys:
            records.append(self.get(k, default=default))
        return records

    def get(self, key, default=C.NULL):
        if default is not C.NULL:
            return self.records.get(key, default)
        return self.records[key]

    def get_origin(self, key, blacklist=None):
        r = self.get(key)
        if r:
            if C.PARENT in r.annotations:
                if blacklist and r.annotations[C.TRANSFORMATION] in blacklist:
                    return r
                else:
                    return self.get_origin(r.annotations[C.PARENT])
        return r

    def add_one(self, record, validate=True):
        if self.key(record):
            return self.key(record)
        else:
            if validate:
                self.validate_records([record])
            key = self.incr()
            self.records[key] = record
            self.mapping[id(record)] = key
            return key

    def add_many(self, records):
        self.validate_records(records)
        keys = []
        for r in records:
            keys.append(self.add_one(r, validate=False))
        return keys
