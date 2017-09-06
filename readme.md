This is a wrapper for initiating blast searches...


# Blast Setup
https://www.ncbi.nlm.nih.gov/books/NBK52640/

Run "install_blash.sh" to install blash

## Usage

To be used as a dropin for BioPython blast command lines, which I found
confusing to use..

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
b.parse_results()

# send your JSON to your other apps
```


## Status

20170908 - This repo is currently in beta...