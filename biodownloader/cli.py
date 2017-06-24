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

import click

import biodownloader

from biodownloader.fetchers import download_structure_from_pdbe
from biodownloader.fetchers import download_sifts_from_ebi
from biodownloader.fetchers import download_data_from_uniprot
from biodownloader.fetchers import download_alignment_from_cath
from biodownloader.fetchers import download_alignment_from_pfam

from biodownloader.config import config

# TODO variants in vcf, json, etc.
# TODO add logging


# https://github.com/pallets/click/issues/108
def add_common(options):
    def _add_options(func):
        for option in reversed(options):
            func = option(func)
        return func
    return _add_options


common_options = [
    click.option('--verbose', '-v', 'verbose', flag_value=2, default=1,
                 help="Verbosity level (via logging)"),
    click.option('--override', 'override',
                 multiple=False, help='Overrides any existing file, if available.',
                 default=False, is_flag=True, required=False),
    click.option('--output', 'output_dir', multiple=False, required=False,
                 help='Directory path to which the files will be written.'),
]


common_arguments = [
    click.argument('ids', nargs=-1, required=True),
]


# main application
@click.group(chain=True,
             context_settings={'help_option_names': ['-h', '--help']})
@click.version_option(version=biodownloader.__version__)
@add_common(common_options)
def downloads(**kwargs):
    """
    BioDownloader: a Command Line Tool for downloading protein
    structures, protein sequences and multiple sequence alignments.

        $ BioDownloader COMMAND --help for additional help
    """
    pass


@downloads.command('pdb')
@click.option('--pdb', 'pdb', multiple=False,
              help='Structure in PDB format (expects PDB ID).',
              default=False, is_flag=True, required=False)
@click.option('--mmcif', 'mmcif', multiple=False,
              help='Structure in mmCIF format (expects PDB ID).',
              default=False, is_flag=True, required=False)
@click.option('--bio', 'bio', multiple=False,
              help=('Preferred BioUnit instead of the asymmetric unit. '
                    'This option only works paired with --mmcif'),
              default=False, is_flag=True, required=False)
@add_common(common_options)
@add_common(common_arguments)
def pdb(ids, pdb=False, mmcif=False, bio=False,
        override=False, output_dir=None, verbose=1):
    """
    Macromolecular structures from the PDBe.

    Pass one or more accession IDs (e.g. '2pah' or '2pah 3kic').
    """

    file_downloader(ids, pdb=pdb, mmcif=mmcif, bio=bio, sifts=False,
                    fasta=False, gff=False, txt=False, cath=False, pfam=False,
                    override=override, output_dir=output_dir, verbose=verbose)


@downloads.command('sifts')
@click.option('--sifts', 'sifts', multiple=False,
              help='SIFTS xml format (expects PDB ID).',
              default=False, is_flag=True, required=False)
@add_common(common_options)
@add_common(common_arguments)
def sifts(ids, sifts=False,
          override=False, output_dir=None, verbose=1):
    """
    SIFTS xml structure-sequence mappings from the EBI.

    Pass one or more accession IDs (e.g. '2pah' or '2pah 3kic').
    """

    file_downloader(ids, pdb=False, mmcif=False, bio=False, sifts=sifts,
                    fasta=False, gff=False, txt=False, cath=False, pfam=False,
                    override=override, output_dir=output_dir, verbose=verbose)


@downloads.command('uniprot')
@click.option('--fasta', 'fasta', multiple=False,
              help='UniProt sequence in fasta format (expects UniProt ID).',
              default=False, is_flag=True, required=False)
@click.option('--gff', 'gff', multiple=False,
              help='UniProt record in gff format (expects UniProt ID).',
              default=False, is_flag=True, required=False)
@click.option('--txt', 'txt', multiple=False,
              help='UniProt record in txt format (expects UniProt ID).',
              default=False, is_flag=True, required=False)
@add_common(common_options)
@add_common(common_arguments)
def uniprot(ids, fasta=False, gff=False, txt=False,
            override=False, output_dir=None, verbose=1):
    """
    Sequences (fasta) and sequence annotations in SwissProt (txt) or
    GFF (gff) format from the UniProt.

    Pass one or more accession IDs (e.g. 'P00439' or 'P00439 P12345').
    """

    file_downloader(ids, pdb=False, mmcif=False, bio=False, sifts=False,
                    fasta=fasta, gff=gff, txt=txt, cath=False, pfam=False,
                    override=override, output_dir=output_dir, verbose=verbose)


@downloads.command('cath')
@click.option('--cath', 'cath', multiple=False,
              help=('CATH Funfam alignment in fasta format '
                    '(expects a CATH <Superfamily>_<Funfam> ID).'),
              default=False, is_flag=True, required=False)
@add_common(common_options)
@add_common(common_arguments)
def cath(ids, cath=False,
         override=False, output_dir=None, verbose=1):
    """
    Multiple sequence alignments (fasta) from CATH.

    Pass one or more accession IDs (e.g. '1.50.10.100_1318').
    """

    file_downloader(ids, pdb=False, mmcif=False, bio=False, sifts=False,
                    fasta=False, gff=False, txt=False, cath=cath, pfam=False,
                    override=override, output_dir=output_dir, verbose=verbose)


@downloads.command('pfam')
@click.option('--pfam', 'pfam', multiple=False,
              help=('Pfam alignment in Stockholm format '
                    '(expects a Pfam ID).'),
              default=False, is_flag=True, required=False)
@add_common(common_options)
@add_common(common_arguments)
def pfam(ids, pfam=False,
         override=False, output_dir=None, verbose=1):
    """
    Multiple sequence alignments (fasta) from Pfam.

    Pass one or more accession IDs (e.g. 'PF08124').
    """

    file_downloader(ids, pdb=False, mmcif=False, bio=False, sifts=False,
                    fasta=False, gff=False, txt=False, cath=False, pfam=pfam,
                    override=override, output_dir=output_dir, verbose=verbose)


# FIXME add verbosity (via logging)
def file_downloader(ids, pdb=False, mmcif=False, bio=False, sifts=False,
                    fasta=False, gff=False, txt=False, cath=False, pfam=False,
                    override=False, output_dir=None, verbose=1):

    if output_dir is not None:
        config.db_root = output_dir
    for pid in ids:

        if pdb:
            download_structure_from_pdbe(pid, pdb=True,
                                         override=override)
        if mmcif:
            download_structure_from_pdbe(pid, pdb=False, bio=False,
                                         override=override)
        if bio:
            download_structure_from_pdbe(pid, pdb=False, bio=True,
                                         override=override)
        if sifts:
            download_sifts_from_ebi(pid, override=override)

        if fasta:
            download_data_from_uniprot(pid, file_format="fasta",
                                       override=override)
        if gff:
            download_data_from_uniprot(pid, file_format="gff",
                                       override=override)
        if txt:
            download_data_from_uniprot(pid, file_format="txt",
                                       override=override)
        if cath:
            download_alignment_from_cath(pid, max_sequences=20000,
                                         override=override)
        if pfam:
            download_alignment_from_pfam(pid, override=override)


if __name__ == '__main__':
    downloads()
