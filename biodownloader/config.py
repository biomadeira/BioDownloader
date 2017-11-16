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

config_defaults = dict()

# Absolute path working dir
config_defaults["db_root"] = "."
# PDBx/mmCIF files
config_defaults["db_pdbx"] = "."
# SIFTS xml files
config_defaults["db_sifts"] = "."
# UniProt dir
config_defaults["db_uniprot"] = "."
# CATH dir
config_defaults["db_cath"] = "."
# Pfam dir
config_defaults["db_pfam"] = "."
# Pickled objects
config_defaults["db_pickled"] = "."

# UniProt HTTP
config_defaults["http_uniprot"] = "http://www.uniprot.org/uniprot/"
# PDBe HTTP
config_defaults["http_pdbe"] = "http://www.ebi.ac.uk/pdbe/"
# CATH HTTP
config_defaults["http_cath"] = "http://www.cathdb.info/version/v4_1_0/"
# SIFTS FTP
config_defaults["ftp_sifts"] = "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/"
# Pfam HTTP
config_defaults["http_pfam"] = "http://pfam.xfam.org/"


class Config(object):
    def __init__(self, config):
        self.config = config
        self._populate_attributes()

    def _populate_attributes(self):
        for key in self.config:
            value = self.config[key]
            setattr(self, key, value)


config = Config(config=config_defaults)
