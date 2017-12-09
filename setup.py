from distutils.core import setup

setup(
    name='pyblast',
    version='0.0.2a',
    packages=['pyblast', 'blast_bin', 'tests'],
    url='',
    license='MIT',
    author='Justin D. Vrana',
    author_email='justin.vrana@gmail.com',
    description='Python wrapper for running BLAST locally',
    install_requires=['biopython'],
    entry_points={
        'console_scripts': [
            'install_pyblast = blast_bin.install_blast:main',
        ]
    }
)
