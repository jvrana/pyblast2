""".. module:: pyblast.

Submodules
==========

.. autosummary::
    :toctree: _autosummary

    BioBlast
    JSONBlast
    BioBlastFactory
    utils.Span
    utils
"""
from .__version__ import __authors__
from .__version__ import __homepage__
from .__version__ import __repo__
from .__version__ import __title__
from .__version__ import __version__
from .blast import BioBlast
from .blast import BioBlastFactory
from .blast import JSONBlast
from .cli import entrypoint
