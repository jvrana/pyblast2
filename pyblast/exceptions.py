class PyBlastWarning(Warning):
    """A generic warning for pyblastbio"""


class PyBlastException(Exception):
    """A generic exception for pyblastbio"""


class SeqRecordValidationError(Exception):
    """An error if a Bio.SeqRecord instance is invalid."""
