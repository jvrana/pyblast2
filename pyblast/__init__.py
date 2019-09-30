"""

.. module:: pyblast

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

from .blast import BioBlast, JSONBlast, BioBlastFactory
from .__version__ import __version__, __title__, __authors__, __homepage__, __repo__
from .cli import entrypoint
