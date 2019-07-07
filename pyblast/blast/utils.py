from .constants import Constants as C


def is_circular(r):
    return (
        r.annotations.get(C.CIRCULAR, False) is True
        or r.annotations.get(C.LINEAR, True) is False
        or r.annotations.get(C.TOPOLOGY, "linear") == "circular"
    )
