import re

r = \
'''
# BLASTN 2.6.0+
# Query: 
# Database: tests/data/blast_results/db
# Fields: query acc., subject acc., score, evalue, bit score, alignment length, identical, gap opens, gaps, query length, q. start, q. end, subject length, s. start, s. end, subject strand, query seq, subject seq
# 105 hits found
some other shit
'''

print(r)

g = re.search(
    '#\s*(?P<blast_ver>.+)\n' +
    '# Query:\s*(?P<query>.*)\n' +
    '# Database:\s*(?P<database>.+)\n' +
    '# Fields:\s*(?P<fields>.+)',
    r)
print(g.groupdict())
