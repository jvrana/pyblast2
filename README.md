![](https://github.com/jvrana/pyblast/workflows/Python%package/badge.svg)
[![Coverage Status](https://coveralls.io/repos/github/jvrana/pyblast/badge.svg?branch=master)](https://coveralls.io/github/jvrana/pyblast?branch=master)
[![PyPI version](https://badge.fury.io/py/pyblastbio.svg)](https://badge.fury.io/py/pyblastbio)

# pyblast

This is a wrapper for other applications to run blast searches on SeqRecord objects and JSON objects. Intended to
be used in small python applications.

Features include:
* Automatic BLAST parsing to JSON
* Alignment to circular queries, using either linear or circular subjects
* Blast self installation

# Installation

You can install BLAST to the pyblast directory using the following command:

```
pyblast install
```

This will install it to pyblast/blast_bin in your python install location. If you want BLAST installed somewhere else, move the *ncbi-blast-X.X.X+* folder
to your desired location and add *path/to/ncbi-blast-X.X.X+/bin* to you $PATH. **PyBlast** will prefer to use the blast stored
in your executable path. If it cannot find a blast executable there, it looks for it in that paths in the pyblast/blast_bin/_paths.txt.
file. _paths.txt is automatically updated when you run install_blast.py so theres no need to manage the paths manually.

After installing and verifying the `blastn` command works from the cmd line,

```
pip install pyblastbio
```

## Usage

This package is a python wrapper for the BLAST command line, intended to be run along with a microservice (e.g. Flask) or for a quick alignment in a jupyter notebook or small python script/app.

This package also includes a basic python-based installation script which is used in unit-testing.

### Running a blast query on a Bio.SeqRecord object

We can do a quick alignment to some sequences using the following, which gives us a nice dictionary of the results:

```python
from pyblast import BioBlast
from pyblast.utils import make_linear, make_circular
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

queries = [
  SeqRecord(Seq("ACGTGATTCGTCGTGTAGTTGAGTGTTACGTTGCATGTCGTACGTGTGTAGTGTCGTGTAGTGCTGATGCTACGTGATCG"))
]
subjects = [
  SeqRecord(Seq("ACGTGATTCGTCGTGTAGTTGAGTGTTACGTTGCATGTCGTTACGTGATCG"))
]

# pyblast requires a 'topology' annotation on the SeqRecords.
# we can make records circular or linear using `make_linear` or `make_circular` methods
subjects = make_linear(subjects)
queries = make_linear(queries)

blast = BioBlast(subjects, queries)
results = blast.quick_blastn()
print(results)
```

```json
[
  {
    "query": {
      "start": 1,
      "end": 46,
      "bases": "ACGTGATTCGTCGTGTAGTTGAGTGTTACGTTGCATGTCGT-ACGTG",
      "strand": 1,
      "length": 80,
      "sequence_id": "11e17df2-579f-4234-a1e6-f4e3fadfe277",
      "circular": false,
      "name": "<unknown name>",
      "origin_key": "bbadd55c-9413-4394-a23c-0da983630b98",
      "origin_record_id": "<unknown id>",
      "origin_sequence_length": 80
    },
    "subject": {
      "start": 1,
      "end": 47,
      "bases": "ACGTGATTCGTCGTGTAGTTGAGTGTTACGTTGCATGTCGTTACGTG",
      "strand": 1,
      "length": 51,
      "sequence_id": "69248d23-1044-4a75-80c9-53b999796d48",
      "circular": false,
      "name": "<unknown name>",
      "origin_key": "1f627d51-93df-458b-ba36-9b5a7b483a4d",
      "origin_record_id": "<unknown id>",
      "origin_sequence_length": 51
    },
    "meta": {
      "query acc.": "11e17df2-579f-4234-a1e6-f4e3fadfe277",
      "subject acc.": "69248d23-1044-4a75-80c9-53b999796d48",
      "score": 43,
      "evalue": 0,
      "bit score": 80,
      "alignment length": 47,
      "identical": 46,
      "gap opens": 1,
      "gaps": 1,
      "query length": 80,
      "q. start": 1,
      "q. end": 46,
      "subject length": 51,
      "s. start": 1,
      "s. end": 47,
      "subject strand": "plus",
      "query seq": "ACGTGATTCGTCGTGTAGTTGAGTGTTACGTTGCATGTCGT-ACGTG",
      "subject seq": "ACGTGATTCGTCGTGTAGTTGAGTGTTACGTTGCATGTCGTTACGTG",
      "span_origin": true
    }
  }
]
```

### Running blast on circular subjects and queries

Pyblast handles alignments to circular subjects and queries as well. As you can see below, we get a complete alignment of the subject (1 to 50) to the circular query (82 over origin to 30). Circular subjects and circular queries can be mixed together, as well as multiple queries.

```
seq = "ACGTTGTAGTGTAGTTGATGATGATGTCTGTGTCGTGTGATGTGCTGTAGTGTTTAGGGGCGGCGCGGAGTATGCTG"
queries = [
	SeqRecord(Seq(seq))
]

subjects = [
	SeqRecord(Seq(seq[-20:] + seq[:30]))
]

# pyblast requires a 'topology' annotation on the SeqRecords.
# we can make records circular or linear using `make_linear` or `make_circular` methods
subjects = make_circular(subjects)
queries = make_circular(queries)

blast = BioBlast(subjects, queries)
results = blast.quick_blastn()
print(results)
```

```json
[
  {
    "query": {
      "start": 82,
      "end": 30,
      "strand": 1,
      "...": "..."
    },
    "subject": {
      "start": 1,
      "end": 50,
      "strand": 1,
      "...": "..."
    },
    "meta": {
    	"...": "..."
    }
]
```

### BioBlastFactory

In some cases, we will want to share the same sequences for different types of alignments. For example, we may want to align a set of primers and a set of templates to the same query records. In these types of cases, we can use the **BioBlastFactory**:

```python
from pyblast import BioBlastFactory

# initialize a new factory
factory = BioBlastFactory()

# add records accessible by keyword
factory.add_records(records1, "primers")
factory.add_records(records2, "templates")
factory.add_records(records3, "queries")

# we spawn new BioBlast alignmers from the keywords above
primer_alignment = factory("primers", "queries")
template_alignment = factory("templates", "queries")

# we can then run alignments, ensuring the queries in both results
# refer to the exact same query
primer_results = primer_alignment.quick_short_blastn()
template_results = template_alignment.quick_blastn()
```

### Utilities for reading files

**pyblast** includes utilities for reading in *fasta* and *genbank* files.

```python
from pyblast.utils import load_glob, load_genbank_glob, load_fasta_glob

# load many genbank files into a list of SeqRecords
# 'topology' is automatically detected here
# we enforce all record_ids to be unique (a requirement for pyblast)
records1 = load_genbank_glob("~/mydesigns/*.gb", force_unique_ids=True)

# load many fasta files into a list of SeqRecords
# 'topology' is NOT detected
# we enforce all record_ids to be unique (a requirement for pyblast)
records2 = make_linear(load_fasta_glob("~/mydesigns/*.fasta"), force_unique_ids=True)

```
