from uuid import uuid4

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from pyblast.exceptions import PyBlastException
from pyblast.utils import new_feature_location


class JSONParser:
    @staticmethod
    def clean_data(data):
        return {k: v for k, v in data.items() if v is not None}

    @classmethod
    def JSON_to_SeqFeature(cls, data, length=None):
        """Convert JSON to a Bio.SeqFeature object."""
        start = data["start"]
        end = data["end"]
        strand = data["strand"]
        return SeqFeature(
            location=new_feature_location(start, end, strand=strand, length=length),
            type=data["type"],
            id=data["name"],
        )

    @classmethod
    def __JSON_to_SeqRecord(cls, data):
        if not issubclass(type(data), dict) and not issubclass(type(data), list):
            raise PyBlastException(
                "Data must be a json object but found a '{}'".format(type(data))
            )
        data = cls.clean_data(data)
        return SeqRecord(
            seq=Seq(data["bases"]),
            annotations={"circular": data.get("circular", False)},
            features=[
                cls.JSON_to_SeqFeature(f, len(data["bases"]))
                for f in data.get("features", [])
            ],
            description=data.get("description", ""),
            name=data["name"],
            id=data.get("id", str(uuid4())),
        )

    @classmethod
    def JSON_to_SeqRecords(cls, datalist):
        """Convert JSON to SeqRecords."""
        records = []
        if isinstance(datalist, list):
            for data in datalist:
                records.append(cls.__JSON_to_SeqRecord(data))
        else:
            records.append(cls.__JSON_to_SeqRecord(datalist))
        return records
