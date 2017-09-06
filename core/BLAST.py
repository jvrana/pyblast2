# More information: https://www.ncbi.nlm.nih.gov/guide/howto/run-blast-local/
# BLAST User Manual: https://www.ncbi.nlm.nih.gov/books/NBK1762/
import tempfile
from .seqio import *
from .utils import *
# custom output format: https://www.ncbi.nlm.nih.gov/books/NBK279682/
class BLAST(object):
    blast_config = {
        "outfmt": "\"7 qacc sacc score evalue bitscore" +
         "length nident gapopen gaps qlen qstart qend" +
                  "slen sstart send sstrand qseq sseq\""
    }

    def __init__(self,
                 name,
                 templates_directory,
                 query,
                 db_output_directory,
                 out,
                 **config):
        self.name = name
        self.templates_directory = os.path.abspath(templates_directory)
        self.query = os.path.abspath(query)
        self.config = BLAST.blast_config
        if config is not None:
            self.config.update(config)
        self.dbdir = os.path.abspath(db_output_directory)
        self.db = os.path.join(db_output_directory, name)
        self.out = out
        self.validate_files()

    def validate_files(self):
        def _is_file(f):
            return os.path.isfile(os.path.abspath(f))

        def _is_dir(d):
            return os.path.isdir(os.path.abspath(d))

        errors = []
        for f in [self.query]:
            if not _is_file(f):
                errors.append("File not found: {}".format(f))
        for d in [self.dbdir, self.templates_directory]:
            if not _is_dir(d):
                errors.append("Directory not found {}".format(d))
        if len(errors) > 0:
            raise ValueError("\n".join(errors))

    def create_config(self):
        d = {
            "db": self.db,
            "out": self.out,
            "query": self.query,
        }
        d.update(self.config)
        return d

    def run(self):
        run_cmd("blastn", **self.create_config())
        with open(self.out, 'rU') as handle:
            self.results = handle.read()

    def concat_templates(self):
        """
        Concatenates sequences in self.templates_directory into a .fsa file while saving important metadata
        :return:
        """
        out = self.db + '.fsa'
        fasta, seqs, metadata = concat_seqs(self.templates_directory, out, savemeta=True)
        self.db_input_metadata = metadata
        self.seqs = seqs
        return out, seqs, metadata

    def makedb(self):
        """
        Concatenates sequencing files from the db_in_dir directory and
        makes a database from the resulting fasta file
        :return: output_path to blast database
        """
        out, seqs, metadata = self.concat_templates()
        return self.fasta_to_db(out)

    def fasta_to_db(self, fasta):
        """
        Makes a blast database from a fasta file
        :param fasta: path to input fasta file
        :param title: title to name database (e.g. "db")
        :param database_output: path to output directory (e.g. "/Users/databases/db")
        :return: path to output directory
        """
        cmd_str = "makeblastdb -in {input} -dbtype nucl -title {title} -out {output}"
        run_cmd("makeblastdb", dbtype="nucl", title=self.name, out=self.db, **{"in": fasta})
        self.db_input_path = fasta
        return self.db

    def parse_results(self, delim=','):
        '''
        This parses a tabulated blast result
        :return:
        '''

        def cleanup_fields(fields, replacements=None):
            if replacements is None:
                replacements = {
                    ('.', ''),
                    (' ', '_'),
                    ('%', 'perc')
                }
            for i, f in enumerate(fields):
                for r in replacements:
                    fields[i] = f.replace(r[0], r[1])
            return fields

        def extract_metadata(r, delim=','):
            g = re.search(
                '#\s*(?P<blast_ver>.+)\n' +
                '# Query:\s*(?P<query>.*)\n' +
                '# Database:\s*(?P<database>.+)\n' +
                '# Fields:\s*(?P<fields>.+)',
                r)
            meta = g.groupdict()
            # clean up fields
            meta['fields'] = re.split('\s*{}\s*'.format(delim), meta['fields'])
            meta['fields'] = cleanup_fields(meta['fields'])
            return g.groupdict()

        def extract_raw_matches(r):
            return re.findall('\n([^#].*)', r)

        def validate_matches(raw_matches, fields):
            match_dicts = []
            for m in raw_matches:
                values = [str_to_f_to_i(v) for v in m.split('\t')]
                match_dicts.append(dict(list(zip(fields, values))))
            return match_dicts
        # print(self.results)
        results = self.results
        if results.strip() == '':
            return {}
        meta = extract_metadata(results)
        fields = meta['fields']
        raw_matches = extract_raw_matches(results)
        match_dicts = validate_matches(raw_matches,fields)
        return match_dicts

    def __str__(self):
        return "{}".format(self.create_config())




class Aligner(BLAST):
    """
    Wrapper for running BLAST. Creates temporary databases and parses results
    """

    def __init__(self, name, templates, design, **additional_params):
        db_out_dir = tempfile.mkdtemp()
        results_out = tempfile.mktemp(dir=db_out_dir)

        super(Aligner, self).__init__(
            name,
            templates,
            design,
            db_out_dir,
            results_out,
            **additional_params
        )

    def run(self, contig_type):
        self.makedbfromdir()
        self.runblast()
        self.contig_container = self.parse_results(contig_type=contig_type)
        return self.contig_container