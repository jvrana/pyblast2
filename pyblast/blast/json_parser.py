from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from uuid import uuid4
from Bio.Seq import Seq


class JSONParser(object):
    @staticmethod
    def clean_data(data):
        return {k: v for k, v in data.items() if v is not None}

    @staticmethod
    def JSON_to_SeqFeature(data, length=None):
        start = data["start"]
        end = data["end"]
        if start > end:
            if length is None:
                raise ValueError(
                    "A length must be provided to create a feature with start > end."
                )
            pass
            f1 = FeatureLocation(start, length, data["strand"])
            f2 = FeatureLocation(0, end, data["strand"])
            if data["strand"] == -1:
                location = CompoundLocation([f2, f1])
            else:
                location = CompoundLocation([f1, f2])
        else:
            location = FeatureLocation(
                data["start"], data["end"], strand=data["strand"]
            )
        return SeqFeature(location=location, type=data["type"], id=data["name"])

    @classmethod
    def __JSON_to_SeqRecord(cls, data):
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
        records = []
        if isinstance(datalist, list):
            for data in datalist:
                records.append(cls.__JSON_to_SeqRecord(data))
        else:
            records.append(cls.__JSON_to_SeqRecord(datalist))
        return records
