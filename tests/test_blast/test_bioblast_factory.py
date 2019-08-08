from pyblast import BioBlastFactory
from pyblast.utils import load_genbank_glob, load_fasta_glob, make_circular, make_linear
from os.path import join


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

    primer_results = primer_blaster.quick_blastn_short()

    template_results = template_blaster.quick_blastn()

    print(len(primer_results))
    print(len(template_results))
