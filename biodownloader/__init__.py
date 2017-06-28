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

import sys
import logging

logging.getLogger("biodownloader").addHandler(logging.NullHandler())
logging.captureWarnings(True)
logging.basicConfig(stream=sys.stderr, level=logging.INFO,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s')


__author__ = "Fábio Madeira"
__copyright__ = "Copyright 2017, Fábio Madeira"
__credits__ = ["Fábio Madeira"]
__license__ = "GNU GPLv3"
import biodownloader.version
__version__ = biodownloader.version.__version__
__maintainer__ = "Fábio Madeira"
__email__ = "fabiomadeira@me.com"
__status__ = "beta"
