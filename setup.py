from distutils.core import setup

setup(
    name='pyblast',
    version='0.0.1',
    packages=['pyblast', 'blast_bin', 'tests'],
    url='',
    license='MIT',
    author='Justin D. Vrana',
    author_email='justin.vrana@gmail.com',
    description='Python wrapper for running BLAST locally',
    install_requires=['biopython'],
)
