import os
import re
from distutils.core import setup

tests_require = [
    'pytest',
    'pytest-cov'
]

install_requires = [
    'marshmallow==3.0.0b13'
]

def parse_version_file():
    here = os.path.abspath(os.path.dirname(__file__))
    ver_dict = {}
    with open(os.path.join(here, 'pyblast', '__version__.py'), 'r') as f:
        for line in f.readlines():
            m = re.match('__(\w+)__\s*=\s*(.+)', line)
            if m:
                ver_dict[m.group(1)] = m.group(2)
    return ver_dict


ver = parse_version_file()

setup(
    title=ver['title'],
    name='pyblast',
    version=ver['version'],
    packages=['pyblast', 'pyblast.blast_bin', 'tests', 'pyblast.utils'],
    url=ver['url'],
    license='MIT',
    author=ver['author'],
    author_email='justin.vrana@gmail.com',
    description='Python wrapper for running BLAST locally',
    install_requires=install_requires,
    tests_require=tests_require,
    entry_points={
        'console_scripts': [
            'pyblast = pyblast.cli:main',
        ]
    }
)
