from .constants import Constants as C
from uuid import uuid4
from pyblast.exceptions import SeqRecordValidationError
import networkx as nx
from .utils import is_circular


class GraphDB(object):
    def __init__(self):
        self.graph = nx.DiGraph()
        self.mapping = {}

    def new_key(self):
        """Generate a new unique key"""
        return str(uuid4())

    def key(self, object):
        """Model id to key. To ensure uniqueness of the models being added."""
        k = self.mapping.get(id(object), None)
        return k

    def get_many(self, keys, default=C.NULL):
        """Get many objects from list of keys."""
        records = []
        for k in keys:
            records.append(self.get(k, default=default))
        return records

    def get(self, key, default=C.NULL):
        """Get object from a key"""
        if default is not C.NULL:
            if key not in self.graph:
                return default
            return self.graph.node[key]["object"]
        return self.graph.node[key]["object"]

    def add(self, object):
        """Add an object"""
        if self.key(object):
            return self.key(object)
        else:
            key = self.new_key()
            self.graph.add_node(key, object=object)
            self.mapping[id(object)] = key
            return key

    def add_many(self, objects):
        """Add many objects"""
        keys = []
        for o in objects:
            keys.append(self.add(o))
        return keys

    def add_edge(self, k1, k2, **kwargs):
        """Add an edge"""
        return self.graph.add_edge(k1, k2, **kwargs)

    def to_dict(self):
        """Return key to object dictionary"""
        d = {}
        for key in self.graph:
            d[key] = self.get(key)
        return d


class SeqRecordDB(GraphDB):
    @staticmethod
    def is_circular(records):
        return is_circular(records)

    @property
    def records(self):
        return self.to_dict()

    @classmethod
    def validate_records(cls, records):
        for r in records:
            cls.validate_circular(r)

    @staticmethod
    def validate_circular(r):
        valid_keys = [C.CIRCULAR, C.LINEAR, C.TOPOLOGY]
        if all([k not in r.annotations for k in valid_keys]):
            raise SeqRecordValidationError(
                "SeqRecord {} is missing a any of the valid keys in its annotation: {keys}. This must be provided.".format(
                    r, keys=valid_keys
                )
            )
        if C.LINEAR in r.annotations and C.CIRCULAR in r.annotations:
            if r.annotations[C.LINEAR] is r.annotations[C.CIRCULAR]:
                raise SeqRecordValidationError(
                    "SeqRecord {} has conflicting topology definitions for '{linear}' and '{circular}'".format(
                        r, circular=C.CIRCULAR, linear=C.LINEAR
                    )
                )

    def transform(self, parent_key, transform, transform_label):
        record = self.get(parent_key)
        new_record = transform(record)
        new_key = self.add(new_record)
        self.add_edge(
            parent_key,
            new_key,
            **{C.TRANSFORMATION: transform_label, C.PARENT: parent_key}
        )
        new_record.id = new_key
        return new_key

    def add_with_transformation(self, record, transform, transform_label):
        key = self.add(record)
        return self.transform(key, transform, transform_label)

    def add_many_with_transformations(self, records, transform, transform_label):
        return [
            self.add_with_transformation(r, transform, transform_label) for r in records
        ]

    def get_origin(self, key, blacklist=None):
        tree = nx.bfs_tree(self.graph, key, reverse=True)
        roots = [n for n, d in dict(tree.out_degree()).items() if d == 0]
        record = self.get(roots[0])
        return record

    def add(self, record, validate=True):
        if self.key(record):
            return self.key(record)
        else:
            if validate:
                self.validate_records([record])
            return super().add(record)

    def add_many(self, records):
        self.validate_records(records)
        return super().add_many(records)
