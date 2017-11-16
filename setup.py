#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
    BioDownloader: a Command Line Tool for downloading protein structures,
    protein sequences and multiple sequence alignments.
    Copyright (C) 2017  Fábio Madeira

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""


from setuptools import setup
from setuptools import find_packages

__author__ = 'Fábio Madeira'
__email__ = 'fabiomadeira@me.com'
__version__ = '0.1.0'


def gather_dependencies():
    import os
    with open('requirements.txt', 'r') as f_in:
        return [l for l in f_in.read().rsplit(os.linesep) 
                if l and not l.startswith("#")]
DEPENDENCIES = gather_dependencies()


setup(
    # Basic package information.
    name='biodownloader',
    version=__version__,
    packages=find_packages(),

    # Packaging options.
    include_package_data=True,

    # Package dependencies
    # should always match the entries in requirements.txt
    install_requires=DEPENDENCIES,

    entry_points={
        "console_scripts": [
            "BioDownloader=biodownloader.cli:downloads",
        ]
    },

    # tests
    test_suite="tests.test_biodownloader",
    tests_require=['mock'],

    # Metadata for PyPI.
    author=__author__,
    author_email=__email__,
    license='GPLv3',
    url='http://github.com/biomadeira/biodownloader/tree/master',
    keywords='download pdb uniprot sifts cath pfam python cli',
    description=('BioDownloader: a Command Line Tool for downloading '
                 'protein structures, protein sequences and '
                 'multiple sequence alignments.'),
    long_description=open('README.rst').read(),
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Internet',
        'Topic :: Scientific/Engineering :: Bio-informatics',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
)

