BioDownloader
=============

|Build Status| |Coverage Status| |License|

A Command Line Tool for downloading protein structures, protein sequences
and multiple sequence alignments.

Setup
~~~~~

Easy install from github using pip.

.. code:: bash

    $ pip install --upgrade http://github.com/biomadeira/BioDownloader/zipball/master



If you want to mess up with the source code.

.. code:: bash

    $ git clone https://github.com/biomadeira/BioDownloader.git
    $ cd BioDownloader
    $ sudo python setup.py install



Quickstart
~~~~~~~~~~

Printing help information...

::

   $ BioDownloader -h
   Usage: BioDownloader [OPTIONS] COMMAND1 [ARGS]... [COMMAND2 [ARGS]...]...

     BioDownloader: a Command Line Tool for downloading protein structures,
     protein sequences and multiple sequence alignments.

         $ BioDownloader COMMAND --help for additional help

   Options:
     --version      Show the version and exit.
     -h, --help     Show this message and exit.

   Commands:
     cath     Multiple sequence alignments (fasta) from...
     pdb      Macromolecular structures from the PDBe.
     pfam     Multiple sequence alignments (fasta) from...
     sifts    SIFTS xml structure-sequence mappings from...
     uniprot  Sequences (fasta) and sequence annotations in...


Printing help information for one of the available commands...

::

   $ BioDownloader uniprot -h
   Usage: BioDownloader uniprot [OPTIONS] IDS...

     Sequences (fasta) and sequence annotations in SwissProt (txt) or GFF (gff)
     format from the UniProt.

     Pass one or more accession IDs (e.g. 'P00439' or 'P00439 P12345').

   Options:
     --fasta        UniProt sequence in fasta format (expects UniProt ID).
     --gff          UniProt record in gff format (expects UniProt ID).
     --txt          UniProt record in txt format (expects UniProt ID).
     -v, --verbosity LVL  Either CRITICAL, ERROR, WARNING, INFO or DEBUG
     --override     Overrides any existing file, if available.
     --output TEXT  Directory path to which the files will be written.
     -h, --help     Show this message and exit.


Downloading a bunch of structure files...

.. code:: bash

    # Downloads structures in PDB and mmCIF format
    $ BioDownloader pdb --pdb --mmcif 2pah 3pah 4pah


Changing where the files will be downloaded to...

.. code:: bash

    # Downloads a UniProt sequence in FASTA and sequence annotations in GFF
    $ BioDownloader uniprot --fasta --gff --output /path/to/output/dir/ P00439



Dependencies
~~~~~~~~~~~~

|Python: versions|

See the necessary `requirements`_ for this module.

Contributing and Bug tracking
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Feel free to fork, clone, share and distribute. If you find any bugs or
issues please log them in the `issue tracker`_.

Contributors:  
 - `Fábio Madeira`_
 - `Felix Simkovic`_


License
~~~~~~~

GNU General Public License v3 (GPLv3). See `license`_ for details.

.. _requirements: https://github.com/biomadeira/BioDownloader/blob/master/requirements.txt
.. _issue tracker: https://github.com/biomadeira/BioDownloader/issues
.. _license: https://github.com/biomadeira/BioDownloader/blob/master/LICENSE.md
.. _Fábio Madeira: https://github.com/biomadeira
.. _Felix Simkovic: https://github.com/fsimkovic

.. |Build Status| image:: https://secure.travis-ci.org/biomadeira/BioDownloader.png?branch=master
   :target: http://travis-ci.org/biomadeira/BioDownloader
.. |Coverage Status| image:: https://coveralls.io/repos/github/biomadeira/BioDownloader/badge.svg?branch=master
   :target: https://coveralls.io/github/biomadeira/BioDownloader?branch=master
.. |License| image:: http://img.shields.io/badge/license-GPLv3-brightgreen.svg?style=flat
   :target: https://github.com/biomadeira/BioDownloader/blob/master/LICENSE.md
.. |Python: versions| image:: https://img.shields.io/badge/python-3.4,_3.5,_3.6-blue.svg?style=flat
   :target: http://travis-ci.org/biomadeira/BioDownloader
