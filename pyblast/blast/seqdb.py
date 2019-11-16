from collections.abc import Sized
from uuid import uuid4

import networkx as nx

from pyblast.constants import Constants as C
from pyblast.exceptions import SeqRecordValidationError
from pyblast.utils import is_circular


class GraphDB(Sized):
    def __init__(self):
        self.graph = nx.DiGraph()
        self.mapping = {}

    @staticmethod
    def new_key():
        """Generate a new unique key."""
        return str(uuid4())

    def key(self, obj):
        """Model id to key.

        To ensure uniqueness of the models being added.
        """
        k = self.mapping.get(id(obj), None)
        return k

    def get_many(self, keys, default=C.NULL):
        """Get many objects from list of keys."""
        records = []
        for k in keys:
            records.append(self.get(k, default=default))
        return records

    def get(self, key, default=C.NULL):
        """Get object from a key."""
        if default is not C.NULL:
            if key not in self.graph:
                return default
            return self.graph.nodes[key]["object"]
        return self.graph.nodes[key]["object"]

    def add(self, obj):
        """Add an object."""
        if self.key(obj):
            return self.key(obj)
        else:
            key = self.new_key()
            self.graph.add_node(key, object=obj)
            self.mapping[id(obj)] = key
            return key

    def add_many(self, objects):
        """Add many objects."""
        keys = []
        for o in objects:
            keys.append(self.add(o))
        return keys

    def add_edge(self, k1, k2, **kwargs):
        """Add an edge."""
        return self.graph.add_edge(k1, k2, **kwargs)

    def to_dict(self):
        """Return key to object dictionary."""
        d = {}
        for key in self.graph:
            d[key] = self.get(key)
        return d

    def __len__(self):
        return len(self.to_dict())


class SeqRecordDB(GraphDB):
    def __init__(self):
        super().__init__()
        self._transformations = {}

    @staticmethod
    def is_circular(records):
        return is_circular(records)

    @property
    def records(self):
        """Return all records."""
        return self.to_dict()

    @classmethod
    def validate_records(cls, records):
        for r in records:
            cls.validate_circular(r)

    @staticmethod
    def validate_circular(r):
        """Ensure the record has a 'topology' of the expected type ['circular',
        'linear']"""

        valid_keys = [C.CIRCULAR, C.LINEAR, C.TOPOLOGY]
        if all([k not in r.annotations for k in valid_keys]):
            raise SeqRecordValidationError(
                "SeqRecord {} is missing a any of the valid keys in its annotation: "
                "{keys}. This must be provided.\n"
                "`Use pyblast.utils.make_circular` or `pyblast.utils.make_linear` to "
                "annotate SeqRecords".format(r, keys=valid_keys)
            )
        if C.LINEAR in r.annotations and C.CIRCULAR in r.annotations:
            if r.annotations[C.LINEAR] is r.annotations[C.CIRCULAR]:
                raise SeqRecordValidationError(
                    "SeqRecord {} has conflicting topology definitions for '{linear}' "
                    "and '{circular}'".format(r, circular=C.CIRCULAR, linear=C.LINEAR)
                )

    def transform(self, parent_key, transform, transform_label):
        """Transform a record.

        :param parent_key:
        :type parent_key:
        :param transform:
        :type transform:
        :param transform_label:
        :type transform_label:
        :return:
        :rtype:
        """
        record = self.get(parent_key)
        new_record = transform(record)
        new_key = self.add(new_record)
        self.add_edge(
            parent_key,
            new_key,
            **{C.TRANSFORMATION: transform_label, C.PARENT: parent_key},
        )
        new_record._transform = transform_label
        new_record.id = new_key
        return new_key

    def add_with_transformation(self, record, transform, transform_label):
        """Add a record and its transformation to the database. If the original
        record exists in the database, then that record is used. If the
        original record has had the same type of transformation applied, no new
        records are added and the previously transformed key is returned. Else,
        the `transform` function is applied and the new record is added, with a
        trace back to the original record.

        :param record: the record to apply a transformation
        :type record: SeqRecord
        :param transform: the transformation function that takes in a single argument,
            the record
        :type transform: callable
        :param transform_label: the transformation label of the transform function
        :type transform_label: str
        :return: the new key that refers to the transformed record
        :rtype: str
        """
        key = self.add(record)
        transform_key = (key, transform_label)
        if transform_key in self._transformations:
            return self._transformations[transform_key]
        else:
            new_key = self.transform(key, transform, transform_label)
            self._transformations[transform_key] = new_key
            return new_key

    def add_many_with_transformations(self, records, transform, transform_label):
        """Transform many records.

        See `add_with_transformation` for usage.
        """
        return [
            self.add_with_transformation(r, transform, transform_label) for r in records
        ]

    def get_origin_key(self, key):
        """Trace the record at 'key' back to its origin through its
        transformation trace and return the origin's key.

        :param key: the key of the record
        :type key: str
        :return: the original record key
        :rtype: str
        """
        tree = nx.bfs_tree(self.graph, key, reverse=True)
        roots = [n for n, d in dict(tree.out_degree()).items() if d == 0]
        return roots[0]

    def get_origin(self, key):
        """Trace the record at 'key' back to its origin through its
        transformation trace and return the origin record.

        :param key: the key of the record
        :type key: str
        :return: the original record
        :rtype: SeqRecord
        """
        return self.get(self.get_origin_key(key))

    def add(self, record, validate=True):
        """Add a record."""
        if self.key(record):
            return self.key(record)
        else:
            if validate:
                self.validate_records([record])
            return super().add(record)

    def add_many(self, records):
        """Add many records."""
        self.validate_records(records)
        return super().add_many(records)
