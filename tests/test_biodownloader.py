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

import os
import re
import sys
import json
import logging
import unittest
import requests
import responses
from click.testing import CliRunner

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
try:
    from mock import patch, MagicMock
except ImportError:
    from unittest.mock import patch, MagicMock

from biodownloader.fetchers import (fetch_from_url_or_retry,
                                    fetch_summary_properties_pdbe,
                                    get_preferred_assembly_id,
                                    download_structure_from_pdbe,
                                    download_sifts_from_ebi,
                                    download_data_from_uniprot,
                                    download_alignment_from_cath,
                                    download_alignment_from_pfam)

from biodownloader.cli import downloads, file_downloader

from biodownloader.config import config as c

from biodownloader.version import __version__

cwd = os.path.abspath(os.path.dirname(__file__))


def response_mocker(kwargs, base_url, endpoint_url, status=200,
                    content_type='application/json', post=False, data=None):
    """
    Generates a mocked requests response for a given set of
    kwargs, base url and endpoint url
    """

    url = re.sub('\{\{(?P<m>[a-zA-Z_]+)\}\}', lambda m: "%s" % kwargs.get(m.group(1)),
                 base_url + endpoint_url)
    with responses.RequestsMock() as rsps:
        if post:
            rsps.add(responses.POST, url,
                     body=b'{"data": "some json formatted output"}',
                     status=status, content_type='application/json')
            response = requests.post(url, data=data)

        elif content_type == 'application/json':
            rsps.add(responses.GET, url,
                     body=b'{"data": "some json formatted output"}',
                     status=status, content_type='application/json')
            response = requests.get(url)
        elif content_type == 'text/plain':
            rsps.add(responses.GET, url,
                     body="Some text-based content\n spanning multiple lines",
                     status=status, content_type='text/plain')
            response = requests.get(url)
        else:
            rsps.add(responses.GET, url,
                     body=b"Some other binary stuff...",
                     status=status, content_type='application/octet-stream')
            response = requests.get(url)
    return response


@patch("biodownloader.config.config.db_root", cwd)
class TestBioDownloader(unittest.TestCase):
    """

    This suit doesn't generally check the data fetched because it will change overtime.
    Here, we only test whether the endpoints still exist or not!

    """

    def setUp(self):
        """Initialize the framework for testing."""

        self.uniprotid = "P00439"
        self.pdbid = "2pah"
        self.cathid = "1.50.10.100_1318"
        self.pfamid = "PF08124"
        self.config = c
        self.tmp = c.db_root
        self.fetch_from_url_or_retry = fetch_from_url_or_retry
        self.fetch_summary_properties_pdbe = fetch_summary_properties_pdbe
        self.get_preferred_assembly_id = get_preferred_assembly_id
        self.download_structure_from_pdbe = download_structure_from_pdbe
        self.download_sifts_from_ebi = download_sifts_from_ebi
        self.download_data_from_uniprot = download_data_from_uniprot
        self.download_alignment_from_cath = download_alignment_from_cath
        self.download_alignment_from_pfam = download_alignment_from_pfam
        self.downloads = downloads
        self.file_downloader = file_downloader

    def tearDown(self):
        """Remove testing framework."""

        self.uniprotid = None
        self.pdbid = None
        self.cathid = None
        self.pfamid = None
        self.config = None
        self.tmp = None
        self.fetch_from_url_or_retry = None
        self.fetch_summary_properties_pdbe = None
        self.get_preferred_assembly_id = None
        self.download_structure_from_pdbe = None
        self.download_sifts_from_ebi = None
        self.download_data_from_uniprot = None
        self.download_alignment_from_cath = None
        self.download_alignment_from_pfam = None
        self.downloads = None
        self.file_downloader = None

    def test_loading_config_defaults(self):
        config = self.config
        self.assertTrue(hasattr(config, 'db_pdbx'))
        self.assertEqual(config.db_pdbx, ".")
        self.assertTrue(hasattr(config, 'db_sifts'))
        self.assertEqual(config.db_sifts, ".")
        self.assertFalse(hasattr(config, 'test'))

    def test_updating_config_defaults(self):
        config = self.config
        config.db_root = 'new_value/'
        self.assertNotEqual(config.db_root, "random_dir/")
        self.assertEqual(config.db_root, "new_value/")

    def test_fetch_from_url_or_retry_get_text(self):
        # mocked requests
        identifier = "P00439"
        base_url = c.http_uniprot
        endpoint_url = "{}.fasta".format(identifier)
        response = response_mocker(kwargs={"P00439.fasta"}, base_url=base_url,
                                   endpoint_url="",
                                   content_type='text/plain')
        self.fetch_from_url_or_retry = MagicMock(return_value=response)
        url = base_url + endpoint_url
        r = self.fetch_from_url_or_retry(url, json=True,
                                         header={'text/plain'}).content
        self.assertEqual(str(r, 'utf-8'),
                         "Some text-based content\n spanning multiple lines")

    def test_fetch_from_url_or_retry_get_json(self):
        # mocked requests
        identifier = "2pah"
        base_url = c.api_pdbe
        endpoint_url = "pdb/entry/summary/"
        response = response_mocker(kwargs={identifier}, base_url=base_url,
                                   endpoint_url=endpoint_url,
                                   content_type='application/json')
        self.fetch_from_url_or_retry = MagicMock(return_value=response)
        url = base_url + endpoint_url + identifier
        r = self.fetch_from_url_or_retry(url, json=True,
                                         header={'application/json'}).json()
        self.assertEqual(r, json.loads('{"data": "some json formatted output"}'))

    def test_fetch_from_url_or_retry_get_binary(self):
        # mocked requests
        identifier = "P00439"
        base_url = c.http_uniprot
        endpoint_url = "{}.fasta".format(identifier)
        response = response_mocker(kwargs={"P00439.fasta"}, base_url=base_url,
                                   endpoint_url="",
                                   content_type='application/octet-stream')
        self.fetch_from_url_or_retry = MagicMock(return_value=response)
        url = base_url + endpoint_url
        r = self.fetch_from_url_or_retry(url, json=True,
                                         header={'application/octet-stream'},
                                         retry_in=None, wait=0,
                                         n_retries=10, stream=False).content
        self.assertEqual(r, b"Some other binary stuff...")

    def test_fetch_from_url_or_retry_post_json(self):
        # mocked requests
        identifier = "1csb, 2pah"
        base_url = c.api_pdbe
        endpoint_url = "pdb/entry/summary/"
        response = response_mocker(kwargs={}, base_url=base_url,
                                   endpoint_url=endpoint_url,
                                   content_type='application/octet-stream',
                                   post=True, data=identifier)
        self.fetch_from_url_or_retry = MagicMock(return_value=response)
        url = base_url + endpoint_url + identifier
        r = self.fetch_from_url_or_retry(url, json=True, post=True, data=identifier,
                                         header={'application/octet-stream'},
                                         retry_in=None, wait=0,
                                         n_retries=10, stream=False).json()
        self.assertEqual(r, json.loads('{"data": "some json formatted output"}'))

    def test_fetch_from_url_or_retry_get_404(self):
        # mocked requests
        identifier = "P00439"
        base_url = c.http_uniprot
        endpoint_url = "{}.fasta".format(identifier)
        response = response_mocker(kwargs={"P00439.fasta"}, base_url=base_url,
                                   endpoint_url="",
                                   content_type='text/plain', status=404)
        self.fetch_from_url_or_retry = MagicMock(return_value=response)
        url = base_url + endpoint_url
        r = self.fetch_from_url_or_retry(url, json=True, header={'text/plain'},
                                         retry_in=None, wait=0,
                                         n_retries=10, stream=False)
        self.assertEqual(r.status_code, 404)
        self.assertFalse(r.ok)

    def test_fetch_from_url_or_retry_get_500(self):
        # mocked requests
        identifier = "P00439"
        base_url = c.http_uniprot
        endpoint_url = "{}.fasta".format(identifier)
        response = response_mocker(kwargs={"P00439.fasta"}, base_url=base_url,
                                   endpoint_url="",
                                   content_type='text/plain', status=500)
        self.fetch_from_url_or_retry = MagicMock(return_value=response)
        url = base_url + endpoint_url
        r = self.fetch_from_url_or_retry(url, json=True, header={'text/plain'},
                                         retry_in=(500,), wait=1,
                                         n_retries=10, stream=False)
        self.assertEqual(r.status_code, 500)
        self.assertFalse(r.ok)

    def test_summary_properties_pdbe(self):
        r = self.fetch_summary_properties_pdbe(self.pdbid)
        self.assertTrue(r.ok)

    def test_summary_properties_cached(self):
        pickled = os.path.join(c.db_root, c.db_pickled, self.pdbid+"_sp.pkl")
        self.assertFalse(os.path.isfile(pickled))
        r = self.fetch_summary_properties_pdbe(self.pdbid, cached=True)
        self.assertTrue(r.ok)
        self.assertTrue(os.path.isfile(pickled))
        os.remove(pickled)

    def test_preferred_assembly_pdbe_1(self):
        r = self.get_preferred_assembly_id(self.pdbid)
        self.assertEqual("1", r)

    def test_download_structure_from_pdbe_pdb_1(self):
        self.download_structure_from_pdbe(self.pdbid, pdb=True)
        os.remove(os.path.join(c.db_root, c.db_pdbx, self.pdbid+".pdb"))

    def test_download_structure_from_pdbe_pdb_2(self):
        self.file_downloader([self.pdbid], pdb=True)
        os.remove(os.path.join(c.db_root, c.db_pdbx, self.pdbid+".pdb"))

    def test_download_structure_from_pdbe_mmcif_1(self):
        self.download_structure_from_pdbe(self.pdbid, pdb=False)
        os.remove(os.path.join(c.db_root, c.db_pdbx, self.pdbid+".cif"))

    def test_download_structure_from_pdbe_mmcif_2(self):
        self.file_downloader([self.pdbid], mmcif=True)
        os.remove(os.path.join(c.db_root, c.db_pdbx, self.pdbid+".cif"))

    def test_download_structure_from_pdbe_mmcif_bio_1(self):
        self.download_structure_from_pdbe(self.pdbid, pdb=False, bio=True)
        os.remove(os.path.join(c.db_root, c.db_pdbx, self.pdbid+"_bio.cif"))

    def test_download_structure_from_pdbe_mmcif_bio_2(self):
        # REM: We are downloading both the mmCIF and assembly!!!
        self.file_downloader([self.pdbid], mmcif=True, bio=True)
        os.remove(os.path.join(c.db_root, c.db_pdbx, self.pdbid+".cif"))
        os.remove(os.path.join(c.db_root, c.db_pdbx, self.pdbid+"_bio.cif"))

    def test_download_sifts_from_ebi_1(self):
        self.download_sifts_from_ebi(self.pdbid)
        os.remove(os.path.join(c.db_root, c.db_sifts, self.pdbid+".xml"))
    
    def test_download_sifts_from_ebi_2(self):
        self.file_downloader([self.pdbid], sifts=True)
        os.remove(os.path.join(c.db_root, c.db_sifts, self.pdbid+".xml"))

    def test_download_data_from_uniprot_fasta_1(self):
        self.download_data_from_uniprot(self.uniprotid, file_format="fasta")
        os.remove(os.path.join(c.db_root, c.db_uniprot, self.uniprotid+".fasta"))
    
    def test_download_data_from_uniprot_fasta_2(self):
        self.file_downloader([self.uniprotid], fasta=True)
        os.remove(os.path.join(c.db_root, c.db_uniprot, self.uniprotid+".fasta"))

    def test_download_data_from_uniprot_gff_1(self):
        self.download_data_from_uniprot(self.uniprotid, file_format="gff")
        os.remove(os.path.join(c.db_root, c.db_uniprot, self.uniprotid+".gff"))
    
    def test_download_data_from_uniprot_gff_2(self):
        self.file_downloader([self.uniprotid], gff=True)
        os.remove(os.path.join(c.db_root, c.db_uniprot, self.uniprotid+".gff"))

    def test_download_data_from_uniprot_txt_1(self):
        self.download_data_from_uniprot(self.uniprotid, file_format="txt")
        os.remove(os.path.join(c.db_root, c.db_uniprot, self.uniprotid+".txt"))
    
    def test_download_data_from_uniprot_txt_2(self):
        self.file_downloader([self.uniprotid], txt=True)
        os.remove(os.path.join(c.db_root, c.db_uniprot, self.uniprotid+".txt"))

    def test_download_alignment_from_cath_1(self):
        self.download_alignment_from_cath(self.cathid)
        os.remove(os.path.join(c.db_root, c.db_cath, self.cathid+".fasta"))

    def test_download_alignment_from_cath_2(self):
        self.file_downloader([self.cathid], cath=True)
        os.remove(os.path.join(c.db_root, c.db_cath, self.cathid+".fasta"))

    def test_download_alignment_from_pfam_1(self):
        self.download_alignment_from_pfam(self.pfamid)
        os.remove(os.path.join(c.db_root, c.db_pfam, self.pfamid+".sth"))
    
    def test_download_alignment_from_pfam_2(self):
        self.file_downloader([self.pfamid], pfam=True)
        os.remove(os.path.join(c.db_root, c.db_pfam, self.pfamid+".sth"))

    def test_cli_version(self):
        runner = CliRunner()
        result = runner.invoke(self.downloads, ['--version'])
        self.assertEqual(result.exit_code, 0)
        self.assertEqual("downloads, version {}\n"
                         "".format(__version__), result.output)

    def test_cli_help(self):
        runner = CliRunner()
        result = runner.invoke(self.downloads, ['-h'])
        self.assertEqual(result.exit_code, 0)
        result = runner.invoke(self.downloads, ['--help'])
        self.assertEqual(result.exit_code, 0)

    def test_cli_pdb_pdb(self):
        runner = CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(self.downloads, ['pdb', '--pdb', '--output',
                                                    self.tmp, self.pdbid])
            self.assertEqual(result.exit_code, 0)

    def test_cli_pdb_mmcif(self):
        runner = CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(self.downloads, ['pdb', '--mmcif', '--output',
                                                    self.tmp, self.pdbid])
            self.assertEqual(result.exit_code, 0)

    def test_cli_pdb_mmcif_bio(self):
        runner = CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(self.downloads, ['pdb', '--mmcif', '--bio',
                                                    '--output', self.tmp, self.pdbid])
            self.assertEqual(result.exit_code, 0)

    def test_cli_sifts_sifts(self):
        runner = CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(self.downloads, ['sifts', '--sifts', '--output',
                                                    self.tmp, self.pdbid])
            self.assertEqual(result.exit_code, 0)

    def test_cli_uniprot_fasta(self):
        runner = CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(self.downloads, ['uniprot', '--fasta', '--output',
                                                    self.tmp, self.uniprotid])
            self.assertEqual(result.exit_code, 0)

    def test_cli_uniprot_txt(self):
        runner = CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(self.downloads, ['uniprot', '--txt', '--output',
                                                    self.tmp, self.uniprotid])
            self.assertEqual(result.exit_code, 0)

    def test_cli_uniprot_gff(self):
        runner = CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(self.downloads, ['uniprot', '--gff', '--output',
                                                    self.tmp, self.uniprotid])
            self.assertEqual(result.exit_code, 0)

    def test_cli_cath_cath(self):
        runner = CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(self.downloads, ['cath', '--cath', '--output',
                                                    self.tmp, self.cathid])
            self.assertEqual(result.exit_code, 0)

    def test_cli_pfam_pfam(self):
        runner = CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(self.downloads, ['pfam', '--pfam', '--output',
                                                    self.tmp, self.pfamid])
            self.assertEqual(result.exit_code, 0)

    def test_cli_pdb_pdb_override(self):
        runner = CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(self.downloads, ['pdb', '--pdb', '--override',
                                                    '--output',
                                                    self.tmp, self.pdbid])
            self.assertEqual(result.exit_code, 0)

if __name__ == '__main__':
    logging.basicConfig(stream=sys.stderr)
    logging.getLogger("biodownloader").setLevel(logging.DEBUG)
    suite = unittest.TestLoader().loadTestsFromTestCase(TestBioDownloader)
    unittest.TextTestRunner(verbosity=2).run(suite)
