from os.path import join

from pyblast import BioBlastFactory
from pyblast.utils import load_fasta_glob
from pyblast.utils import load_genbank_glob
from pyblast.utils import make_linear


def test_bioblast_factory_init(here):
    subjects = load_genbank_glob(
        join(here, "data/test_data/genbank/templates/*.gb"), force_unique_ids=True
    )
    queries = load_genbank_glob(
        join(here, "data/test_data/genbank/designs/*.gb"), force_unique_ids=True
    )
    primers = load_fasta_glob(join(here, "data/test_data/primers/*.fasta"))

    factory = BioBlastFactory()
    factory.add_records(make_linear(primers), "primers")
    factory.add_records(queries, "queries")
    factory.add_records(subjects, "subjects")

    primer_blaster = factory("primers", "queries")
    template_blaster = factory("subjects", "queries")

    primer_results = primer_blaster.blastn_short()

    template_results = template_blaster.blastn()

    print(len(primer_results))
    print(len(template_results))


def test_validate_rc(here):
    queries = load_genbank_glob(
        join(here, "data/test_data/genbank/designs/*.gb"), force_unique_ids=True
    )

    templates = make_linear(queries[:1])
    queries = make_linear([templates[0].reverse_complement()])

    factory = BioBlastFactory()
    factory.add_records(queries, "queries")
    factory.add_records(templates, "templates")

    blaster = factory("templates", "queries")

    results = blaster.blastn()

    assert results[0]["subject"]["strand"] == -1
    assert results[0]["query"]["strand"] == 1
    assert results[0]["subject"]["start"] == len(queries[0])
