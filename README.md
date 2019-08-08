
[![travis build](https://img.shields.io/travis/jvrana/pyblast.svg)](https://travis-ci.org/jvrana/pyblast)
[![Coverage Status](https://coveralls.io/repos/github/jvrana/pyblast/badge.svg?branch=master)](https://coveralls.io/github/jvrana/pyblast?branch=master)
[![PyPI version](https://badge.fury.io/py/REPO.svg)](https://badge.fury.io/py/REPO)

![module_icon](images/module_icon.png?raw=true)

#### Build/Coverage Status
Branch | Build | Coverage
:---: | :---: | :---:
**master** | [![travis build](https://img.shields.io/travis/jvrana/pyblast/master.svg)](https://travis-ci.org/jvrana/pyblast/master) | [![Coverage Status](https://coveralls.io/repos/github/jvrana/pyblast/badge.svg?branch=master)](https://coveralls.io/github/jvrana/pyblast?branch=master)
**development** | [![travis build](https://img.shields.io/travis/jvrana/pyblast/development.svg)](https://travis-ci.org/jvrana/pyblast/development) | [![Coverage Status](https://coveralls.io/repos/github/jvrana/pyblast/badge.svg?branch=development)](https://coveralls.io/github/jvrana/pyblast?branch=development)

**this repo is not longer active**

# pyblast

This is a wrapper for other applications to run blast searches on SeqRecord objects and JSON objects.
Features include:
* Blast self installation
* Alignment to circular queries, using either linear or circular subjects

# Installation

You can install BLAST to the pyblast directory using the following command:

```
pyblast install
```

This will install it to pyblast/blast_bin. If you want BLAST installed somewhere else, move the *ncbi-blast-X.X.X+* folder
to your desired location and add *path/to/ncbi-blast-X.X.X+/bin* to you $PATH. **PyBlast** will prefer to use the blast stored
in your executable path. If it cannot find a blast executable there, it looks for it in that paths in the pyblast/blast_bin/_paths.txt.
file. _paths.txt is automatically updated when you run install_blast.py so theres no need to manage the paths manually.

After installing and verifying the `blastn` command works from the cmd line,

```
pip install pyblastbio
```

## Usage

This package is a python wrapper for the BLAST command line, intended to be run along with a microservice (e.g. Flask).
This package also includes a basic python-based installation script which is used in unit-testing.

### Input options

```python

from pyblast import BioBlast
from Bio.SeqRecord
