"""dna_bases.py"""

s1 = "ATUGCYRSWKMBDHVN"
s2 = "TAACGRYSWMKVHDBN"

rc_dict = dict(zip(s1, s2))
rc_dict.update(dict(zip(s1.lower(), s2.lower())))

ambiguous = {
    "Y": ["C", "T"],
    "R": ["A", "G"],
    "S": ["G", "C"],
    "W": ["A", "T"],
    "K": ["T", "G"],
    "M": ["A", "C"],
    "B": ["C", "G", "T"],
    "D": ["A", "G", "T"],
    "H": ["A", "C", "T"],
    "V": ["A", "C", "G"],
    "N": ["A", "C", "G", "T"],
}
