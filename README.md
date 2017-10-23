
[![travis build](https://img.shields.io/travis/jvrana/pyblast.svg)](https://travis-ci.org/jvrana/pyblast)
[![Coverage Status](https://coveralls.io/repos/github/jvrana/pyblast/badge.svg?branch=master)](https://coveralls.io/github/jvrana/pyblast?branch=master)
[![PyPI version](https://badge.fury.io/py/REPO.svg)](https://badge.fury.io/py/REPO)

![module_icon](images/module_icon.png?raw=true)

#### Build/Coverage Status
Branch | Build | Coverage
:---: | :---: | :---:
**master** | [![travis build](https://img.shields.io/travis/jvrana/pyblast/master.svg)](https://travis-ci.org/jvrana/pyblast/master) | [![Coverage Status](https://coveralls.io/repos/github/jvrana/pyblast/badge.svg?branch=master)](https://coveralls.io/github/jvrana/pyblast?branch=master)
**development** | [![travis build](https://img.shields.io/travis/jvrana/pyblast/development.svg)](https://travis-ci.org/jvrana/pyblast/development) | [![Coverage Status](https://coveralls.io/repos/github/jvrana/pyblast/badge.svg?branch=development)](https://coveralls.io/github/jvrana/pyblast?branch=development)


# Blast Setup

You can install BLAST to the pyblast directory using the following command:

```
python blast_bin/install_blast your_email
```



## Status

20170908 - This repo is currently in beta...

## Usage

To be used as a dropin for BioPython blast command lines, which I found
annoying to use...

```python
# files and directories
db_name = 'db'
templates = 'tests/data/test_data/templates'
db_out_dir = 'tests/data/blast_results'
results_out = 'tests/data/blast_results/results.out'

# define blast configuration
b = BLAST(db_name, templates, db_out_dir, results_out)

# concatenate the sequences in templates directory and create a database
b.makedb()

# run the search
b.run()

# parse the matches into JSON and save
results = b.parse_results()

# send your JSON to your other apps
"""
results = \
[
  {
    "query_acc": "Query_1",
    "subject_acc": "pMODKan-HO-pACT1-ZEV4",
    "score": 4219,
    "evalue": 0,
    "identical": 4219,
    "gap_opens": 0,
    "gaps": 0,
    "query_length": 10781,
    "q_start": 1374,
    "s_start": 1,
    "s_end": 4219,
    "subject_strand": "plus",
    "query_seq": "TCGCGCGTTTCGGTGATGACGGTGAAAACCTCTGACACATGCAGCTCCCGGAGACGGTCACAGCTTGTCTGTAAGCGGATGCCGGGAGCAGACAAGCCCGTCAGGGCGCGTCAGCGGGTGTTGGCGGGTGTCGGGGCTGGCTTAACTATGCGGCGTTTAAACTTAGCAGATGCGCGCACCTGCGTTGTTACCACAACTCTTATGAGGCCCGCGGACAGCATCAAACTGTAAGATTCCGCCACATTTTATACACTCTGGTCCTTTAACTGGCAAACCTTCGGGCGTAATGCCCAATTTTTCGCCTTTGTCTTTTGCCTTTTTCACTTCACGTGCTTCTGGTACATACTTGCAATTTATACAGTGATGACCGCTGAATTTGTATCTTCCATAGCATCTAGCACATACTCGATTTTTACCACTCCAATCTTTATAAAAATACTTGATTCCCTTTCTGGGACAAGCAACACAGTGTTTTAGATTCTTTTTTTGTGATATTTTAAGCTGTTCTCCCACACAGCAGCCTCGACATGATTTCACTTCTATTTTGTTGCCAAGCAAGAAATTTTTATGGCCTTCTATCGTAAGCCCATATACAGTACTCTCACCCTGGAAATCATCCGTGAAGCTGAAATATACGGGTTCCCTTTTTATAATTGGCGGAACTTCTCTTGTTTTGTGACCACTTCGACAATATGACAAAACATTCTGTGAAGTTGTTCCCCCAGCATCAGAGCAGATTGTACTGAGAGTGCACCGGCGCGCCAGATCTGTTTAGCTTGCCTCGTCCCCGCCGGGTCACCCGGCCAGCGACATGGGGGCCCAGAATACCCTCCTTGACAGTCTTGACGTGCGCAGCTCAGGGGCATGATGTGACTGTCGCCCGTACATTTAGCCCATACATCCCCATGTATAATCATTTGCATCCATACATTTTGATGGCCGCACGGCGCGAAGCAAAAATTACGGCTCCTCGCTGCAGACCTGCGAGCAGGGAAACGCTCCCCTCACAGACGCGTTGAATTGTCCCCACGCCGCGCCCCTGTAGAGAAATATAAAAGGTTAGGATTTGCCACTGAGGTTCTTCTTTCATATACTTCCTTTTAAAATCTTGCTAGGATACAGTTCTCACATCACATCCGAACATAAACAACCATGGGTAAGGAAAAGACTCACGTTTCGAGGCCGCGATTAAATTCCAACATGGATGCTGATTTATATGGGTATAAATGGGCTCGCGATAATGTCGGGCAATCAGGTGCGACAATCTATCGATTGTATGGGAAGCCCGATGCGCCAGAGTTGTTTCTGAAACATGGCAAAGGTAGCGTTGCCAATGATGTTACAGATGAGATGGTCAGACTAAACTGGCTGACGGAATTTATGCCTCTTCCGACCATCAAGCATTTTATCCGTACTCCTGATGATGCATGGTTACTCACCACTGCGATCCCCGGCAAAACAGCATTCCAGATATTAGAAGAATATCCTGATTCAGGTGAAAATATTGTTGATGCGCTGGCAGTGTTCCTGCGCCGGTTGCATTCGATTCCTGTTTGTAATTGTCCTTTTAACAGCGATCGCGTATTTCGTCTCGCTCAGGCGCAATCACGAATGAATAACGGTTTGGTTGATGCGAGTGATTTTGATGACGAGCGTAATGGCTGGCCTGTTGAACAAGTCTGGAAAGAAATGCATAAGCTTTTGCCATTCTCACCGGATTCAGTCGTCACTCATGGTGATTTCTCACTTGATAACCTTATTTTTGACGAGGGGAAATTAATAGGTTGTATTGATGTTGGACGAGTCGGAATCGCAGACCGATACCAGGATCTTGCCATCCTATGGAACTGCCTCGGTGAGTTTTCTCCTTCATTACAGAAACGGCTTTTTCAAAAATATGGTATTGATAATCCTGATATGAATAAATTGCAGTTTCATTTGATGCTCGATGAGTTTTTCTAATCAGTACTGACAATAAAAAGATTCTTGTTTTCAAGAACTTGTCATTTGTATAGTTTTTTTATATTGTAGTTGTTCTATTTTAATCAAATGTTAGCGTGATTTATATTTTTTTTCGCCTCGACATCATCTGCCCAGATGCGAAGTTAAGTGCGCAGAAAGTAATATCATGCGTCAATCGTATGTGAATGCTGGTCGCTATACTGCTGTCGATTCGATACTAACGCCGCCATCCAGTGTCCGCCAGGGTTTTCCCAGTCACGACGCCTCTACCTTGCAGACCCATATAATATAATAACTAAATAAGTAAATAAGACACACGCGAGAACATATATACACAATTACAGTAACAATAACAAGAGGACAGATACTACCAAAATGTGTGGGGAAGCGGGTAAGCTGCCACAGCAATTAATGCACAACATTTAACCTACATTCTTCCTTATCGGATCCTCAAAACCCTTAAAAACATATGCCTCACCCTAACATATTTTCCAATTAACCCTCAATATTTCTCTGTCACCCGGCCTCTATTTTCCATTTTCTTCTTTACCCGCCACGCGTTTTTTTCTTTCAAATTTTTTTCTTCCTTCTTCTTTTTCTTCCACGTCCTCTTGCATAAATAAATAAACCGTTTTGAAACCAAACTCGCCTCTCTCTCTCCTTTTTGAAATATTTTTGGGTTTGTTTGATCCTTTCCTTCCCAATCTCTCTTGTTTAATATATATTCATTTATATCACGCTCTCTTTTTATCTTCCTTTTTTTCCTCTCTCTTGTATTCTTCCTTCCCCTTTCTACTCAAACCAAGAAGAAAAAGAAAAGGTCAATCTTTGTTAAAGAATAGGATCTTCTACTACATCAGCTTTTAGATTTTTCACGCTTACTGCTTTTTTCTTCCCAAGATCGAAAATTTACTGAATTAACAGGGCCCCCCCTCGAGGTCGACGGTATCGATAAGCTTGAAGCAAGCCTCCTGAAAGATGGGTACCCGCCCATATGCTTGCCCTGTCGAGTCCTGCGATCGCCGCTTTTCTCGCCACGCCAATCTTACCCGCCATATCCGCATCCATACCGGTCAGAAGCCCTTCCAGTGTCGAATCTGCATGCGTAACTTCAGTCGTAATGCGAACCTTGTGCGCCACATCCGCACCCACACAGGATCCCAAAAGCCGTTCCAATGTCGGATCTGTATGCGGAACTTTAGTCGAAAGGCCGACCTGAGGCGTCACATTCGCACGCACACCGGCGAGAAGCCTTTTGCCTGTGACATTTGTGGGAGGAAGTTTGCCAGGAAGGGCGACCTCAAGAGGCATACCAAAATCCATACAGGTGGCGGAGGCACACCTGCAGCTGCGTCGACTCTAGAGGATCCATCTGCTGGAGACATGAGAGCTGCCAACCTTTGGCCAAGCCCGCTCATGATCAAACGCTCTAAGAAGAACAGCCTGGCCTTGTCCCTGACGGCCGACCAGATGGTCAGTGCCTTGTTGGATGCTGAGCCCCCCATACTCTATTCCGAGTATGATCCTACCAGACCCTTCAGTGAAGCTTCGATGATGGGCTTACTGACCAACCTGGCAGACAGGGAGCTGGTTCACATGATCAACTGGGCGAAGAGGGTGCCAGGCTTTGTGGATTTGACCCTCCATGATCAGGTCCACCTTCTAGAATGTGCCTGGCTAGAGATCCTGATGATTGGTCTCGTCTGGCGCTCCATGGAGCACCCAGTGAAGCTACTGTTTGCTCCTAACTTGCTCTTGGACAGGAACCAGGGAAAATGTGTAGAGGGCATGGTGGAGATCTTCGACATGCTGCTGGCTACATCATCTCGGTTCCGCATGATGAATCTGCAGGGAGAGGAGTTTGTGTGCCTCAAATCTATTATTTTGCTTAATTCTGGAGTGTACACATTTCTGTCCAGCACCCTGAAGTCTCTGGAAGAGAAGGACCATATCCACCGAGTCCTGGACAAGATCACAGACACTTTGATCCACCTGATGGCCAAGGCAGGCCTGACCCTGCAGCAGCAGCACCAGCGGCTGGCCCAGCTCCTCCTCATCCTCTCCCACATCAGGCACATGAGTAACAAAGGCATGGAGCATCTGTACAGCATGAAGTGCAAGAACGTGGTGCCCCTCTATGACCTGCTGCTGGAGATGCTGGACGCCCACCGCCTACATGCGCCCACTAGCCGTGGAGGGGCATCCGTGGAGGAGACGGACCAAAGCCACTTGGCCACTGCGGGCTCTACTTCATCGG",
    "subject_seq": "TCGCGCGTTTCGGTGATGACGGTGAAAACCTCTGACACATGCAGCTCCCGGAGACGGTCACAGCTTGTCTGTAAGCGGATGCCGGGAGCAGACAAGCCCGTCAGGGCGCGTCAGCGGGTGTTGGCGGGTGTCGGGGCTGGCTTAACTATGCGGCGTTTAAACTTAGCAGATGCGCGCACCTGCGTTGTTACCACAACTCTTATGAGGCCCGCGGACAGCATCAAACTGTAAGATTCCGCCACATTTTATACACTCTGGTCCTTTAACTGGCAAACCTTCGGGCGTAATGCCCAATTTTTCGCCTTTGTCTTTTGCCTTTTTCACTTCACGTGCTTCTGGTACATACTTGCAATTTATACAGTGATGACCGCTGAATTTGTATCTTCCATAGCATCTAGCACATACTCGATTTTTACCACTCCAATCTTTATAAAAATACTTGATTCCCTTTCTGGGACAAGCAACACAGTGTTTTAGATTCTTTTTTTGTGATATTTTAAGCTGTTCTCCCACACAGCAGCCTCGACATGATTTCACTTCTATTTTGTTGCCAAGCAAGAAATTTTTATGGCCTTCTATCGTAAGCCCATATACAGTACTCTCACCCTGGAAATCATCCGTGAAGCTGAAATATACGGGTTCCCTTTTTATAATTGGCGGAACTTCTCTTGTTTTGTGACCACTTCGACAATATGACAAAACATTCTGTGAAGTTGTTCCCCCAGCATCAGAGCAGATTGTACTGAGAGTGCACCGGCGCGCCAGATCTGTTTAGCTTGCCTCGTCCCCGCCGGGTCACCCGGCCAGCGACATGGGGGCCCAGAATACCCTCCTTGACAGTCTTGACGTGCGCAGCTCAGGGGCATGATGTGACTGTCGCCCGTACATTTAGCCCATACATCCCCATGTATAATCATTTGCATCCATACATTTTGATGGCCGCACGGCGCGAAGCAAAAATTACGGCTCCTCGCTGCAGACCTGCGAGCAGGGAAACGCTCCCCTCACAGACGCGTTGAATTGTCCCCACGCCGCGCCCCTGTAGAGAAATATAAAAGGTTAGGATTTGCCACTGAGGTTCTTCTTTCATATACTTCCTTTTAAAATCTTGCTAGGATACAGTTCTCACATCACATCCGAACATAAACAACCATGGGTAAGGAAAAGACTCACGTTTCGAGGCCGCGATTAAATTCCAACATGGATGCTGATTTATATGGGTATAAATGGGCTCGCGATAATGTCGGGCAATCAGGTGCGACAATCTATCGATTGTATGGGAAGCCCGATGCGCCAGAGTTGTTTCTGAAACATGGCAAAGGTAGCGTTGCCAATGATGTTACAGATGAGATGGTCAGACTAAACTGGCTGACGGAATTTATGCCTCTTCCGACCATCAAGCATTTTATCCGTACTCCTGATGATGCATGGTTACTCACCACTGCGATCCCCGGCAAAACAGCATTCCAGATATTAGAAGAATATCCTGATTCAGGTGAAAATATTGTTGATGCGCTGGCAGTGTTCCTGCGCCGGTTGCATTCGATTCCTGTTTGTAATTGTCCTTTTAACAGCGATCGCGTATTTCGTCTCGCTCAGGCGCAATCACGAATGAATAACGGTTTGGTTGATGCGAGTGATTTTGATGACGAGCGTAATGGCTGGCCTGTTGAACAAGTCTGGAAAGAAATGCATAAGCTTTTGCCATTCTCACCGGATTCAGTCGTCACTCATGGTGATTTCTCACTTGATAACCTTATTTTTGACGAGGGGAAATTAATAGGTTGTATTGATGTTGGACGAGTCGGAATCGCAGACCGATACCAGGATCTTGCCATCCTATGGAACTGCCTCGGTGAGTTTTCTCCTTCATTACAGAAACGGCTTTTTCAAAAATATGGTATTGATAATCCTGATATGAATAAATTGCAGTTTCATTTGATGCTCGATGAGTTTTTCTAATCAGTACTGACAATAAAAAGATTCTTGTTTTCAAGAACTTGTCATTTGTATAGTTTTTTTATATTGTAGTTGTTCTATTTTAATCAAATGTTAGCGTGATTTATATTTTTTTTCGCCTCGACATCATCTGCCCAGATGCGAAGTTAAGTGCGCAGAAAGTAATATCATGCGTCAATCGTATGTGAATGCTGGTCGCTATACTGCTGTCGATTCGATACTAACGCCGCCATCCAGTGTCCGCCAGGGTTTTCCCAGTCACGACGCCTCTACCTTGCAGACCCATATAATATAATAACTAAATAAGTAAATAAGACACACGCGAGAACATATATACACAATTACAGTAACAATAACAAGAGGACAGATACTACCAAAATGTGTGGGGAAGCGGGTAAGCTGCCACAGCAATTAATGCACAACATTTAACCTACATTCTTCCTTATCGGATCCTCAAAACCCTTAAAAACATATGCCTCACCCTAACATATTTTCCAATTAACCCTCAATATTTCTCTGTCACCCGGCCTCTATTTTCCATTTTCTTCTTTACCCGCCACGCGTTTTTTTCTTTCAAATTTTTTTCTTCCTTCTTCTTTTTCTTCCACGTCCTCTTGCATAAATAAATAAACCGTTTTGAAACCAAACTCGCCTCTCTCTCTCCTTTTTGAAATATTTTTGGGTTTGTTTGATCCTTTCCTTCCCAATCTCTCTTGTTTAATATATATTCATTTATATCACGCTCTCTTTTTATCTTCCTTTTTTTCCTCTCTCTTGTATTCTTCCTTCCCCTTTCTACTCAAACCAAGAAGAAAAAGAAAAGGTCAATCTTTGTTAAAGAATAGGATCTTCTACTACATCAGCTTTTAGATTTTTCACGCTTACTGCTTTTTTCTTCCCAAGATCGAAAATTTACTGAATTAACAGGGCCCCCCCTCGAGGTCGACGGTATCGATAAGCTTGAAGCAAGCCTCCTGAAAGATGGGTACCCGCCCATATGCTTGCCCTGTCGAGTCCTGCGATCGCCGCTTTTCTCGCCACGCCAATCTTACCCGCCATATCCGCATCCATACCGGTCAGAAGCCCTTCCAGTGTCGAATCTGCATGCGTAACTTCAGTCGTAATGCGAACCTTGTGCGCCACATCCGCACCCACACAGGATCCCAAAAGCCGTTCCAATGTCGGATCTGTATGCGGAACTTTAGTCGAAAGGCCGACCTGAGGCGTCACATTCGCACGCACACCGGCGAGAAGCCTTTTGCCTGTGACATTTGTGGGAGGAAGTTTGCCAGGAAGGGCGACCTCAAGAGGCATACCAAAATCCATACAGGTGGCGGAGGCACACCTGCAGCTGCGTCGACTCTAGAGGATCCATCTGCTGGAGACATGAGAGCTGCCAACCTTTGGCCAAGCCCGCTCATGATCAAACGCTCTAAGAAGAACAGCCTGGCCTTGTCCCTGACGGCCGACCAGATGGTCAGTGCCTTGTTGGATGCTGAGCCCCCCATACTCTATTCCGAGTATGATCCTACCAGACCCTTCAGTGAAGCTTCGATGATGGGCTTACTGACCAACCTGGCAGACAGGGAGCTGGTTCACATGATCAACTGGGCGAAGAGGGTGCCAGGCTTTGTGGATTTGACCCTCCATGATCAGGTCCACCTTCTAGAATGTGCCTGGCTAGAGATCCTGATGATTGGTCTCGTCTGGCGCTCCATGGAGCACCCAGTGAAGCTACTGTTTGCTCCTAACTTGCTCTTGGACAGGAACCAGGGAAAATGTGTAGAGGGCATGGTGGAGATCTTCGACATGCTGCTGGCTACATCATCTCGGTTCCGCATGATGAATCTGCAGGGAGAGGAGTTTGTGTGCCTCAAATCTATTATTTTGCTTAATTCTGGAGTGTACACATTTCTGTCCAGCACCCTGAAGTCTCTGGAAGAGAAGGACCATATCCACCGAGTCCTGGACAAGATCACAGACACTTTGATCCACCTGATGGCCAAGGCAGGCCTGACCCTGCAGCAGCAGCACCAGCGGCTGGCCCAGCTCCTCCTCATCCTCTCCCACATCAGGCACATGAGTAACAAAGGCATGGAGCATCTGTACAGCATGAAGTGCAAGAACGTGGTGCCCCTCTATGACCTGCTGCTGGAGATGCTGGACGCCCACCGCCTACATGCGCCCACTAGCCGTGGAGGGGCATCCGTGGAGGAGACGGACCAAAGCCACTTGGCCACTGCGGGCTCTACTTCATCGG"
  },
  ...]:
 """
```

