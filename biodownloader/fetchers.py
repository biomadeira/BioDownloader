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
import sys
import gzip
import time
import shutil
import urllib
import pickle
import requests

from biodownloader.config import config


def flash(message):
    """
    Flashes a message out.
    :param message: input message str()
    """
    print(str(message))
    sys.stdout.flush()


def fetch_from_url_or_retry(url, json=True, header=None, post=False, data=None,
                            retry_in=None, wait=1, n_retries=10, stream=False, **params):
    """
    Fetch an url using Requests or retry fetching it if the server is
    complaining with retry_in error. There is a limit to the number of retries.

    Retry code examples: 429, 500 and 503

    :param url: url to be fetched as a string
    :param json: json output
    :param header: dictionary
    :param post: boolean
    :param data: dictionary: only if post is True
    :param retry_in: http codes for retrying
    :param wait: sleeping between tries in seconds
    :param n_retries: number of retry attempts
    :param stream: boolean
    :param params: request.get kwargs.
    :return: url content
    """

    if retry_in is None:
        retry_in = ()
    else:
        assert type(retry_in) is tuple or type(retry_in) is list

    if header is None:
        header = {}
    else:
        assert type(header) is dict

    if json:
        header.update({"Content-Type": "application/json"})
    else:
        if "Content-Type" not in header:
            header.update({"Content-Type": "text/plain"})

    if post:
        if data is not None:
            assert type(data) is dict or type(data) is str
            response = requests.post(url, headers=header, data=data)
        else:
            return None
    else:
        response = requests.get(url, headers=header, params=params, stream=stream)

    if response.ok:
        return response
    elif response.status_code in retry_in and n_retries >= 0:
        time.sleep(wait)
        return fetch_from_url_or_retry(url, json, header, post, data, retry_in, wait,
                                       (n_retries - 1), stream, **params)
    else:
        try:
            response.raise_for_status()
        except requests.exceptions.HTTPError as e:
            message = '{}: Unable to retrieve {} for {}'.format(response.status_code,
                                                                url, str(e))
            flash(message)


def fetch_summary_properties_pdbe(identifier, cached=False, retry_in=(429,)):
    """
    Queries the PDBe API to get summary validation report.

    :param identifier: PDB ID
    :param cached: (boolean) if True, stores a json file locally
    :param retry_in: http code for retrying connections
    :return: response object
    """

    pickled = "{}{}{}_sp.pkl".format(config.db_root, config.db_pickled, identifier)
    if cached and os.path.isfile(pickled):
        response = pickle.load(open(pickled, 'rb'))
    else:
        pdbe_endpoint = "pdb/entry/summary/"
        url = config.api_pdbe + pdbe_endpoint + identifier
        response = fetch_from_url_or_retry(url, retry_in=retry_in, json=True)
        if cached:
            with open(pickled, 'wb') as output:
                pickle.dump(response, output, -1)
    return response


def fetch_preferred_assembly_id(identifier):
    """
    Gets the preferred assembly id for the given PDB ID, from the PDBe API.

    :param identifier: PDB ID
    :return: (str)
    """

    # getting the preferred biological assembly from the PDBe API
    pref_assembly = "1"
    try:
        data = fetch_summary_properties_pdbe(identifier)
    except Exception as e:
        message = "Something went wrong for {}... {}".format(identifier, e)
        flash(message)
    try:
        if data is not None:
            data = data.json()
            nassemblies = data[identifier][0]["assemblies"]
            if len(nassemblies) > 1:
                for entry in nassemblies:
                    if entry["preferred"]:
                        pref_assembly = entry["assembly_id"]
                        break
            else:
                pref_assembly = data[identifier][0]["assemblies"][0]["assembly_id"]
    except Exception as e:
        pref_assembly = "1"
        message = "Something went wrong for {}... {}".format(identifier, e)
        flash(message)

    bio_best = str(pref_assembly)
    return bio_best


def download_structure_from_pdbe(identifier, pdb=False, bio=False, override=False):
    """
    Downloads a structure from the PDBe to the filesystem.

    :param identifier: (str) PDB ID
    :param pdb: (boolean) PDB formatted if True, otherwise mmCIF format
    :param bio: (boolean) if true downloads the preferred Biological Assembly
    :param override: (boolean)
    :return: (side effects)
    """

    if pdb:
        filename = "{}.pdb".format(identifier)
    else:
        if bio:
            filename = "{}_bio.cif.gz".format(identifier)
        else:
            filename = "{}.cif".format(identifier)

    outputfile = "{}{}{}".format(config.db_root, config.db_pdbx, filename)
    os.makedirs("{}{}".format(config.db_root, config.db_pdbx), exist_ok=True)
    if not os.path.exists(outputfile.rstrip('.gz')) or override:
        if pdb:
            url = config.http_pdbe + ("entry-files/download/"
                                      "pdb{}.ent".format(identifier))
        else:
            if bio:
                pref = fetch_preferred_assembly_id(identifier=identifier)
                endpoint = ("static/entry/download/"
                            "{}-assembly-{}.cif.gz".format(identifier, pref))
                # atom lines only?
                # endpoint = "{}-assembly-{}_atom_site.cif.gz".format(pdbid, pref)
                url = config.http_pdbe + endpoint
            else:
                # original mmCIF?
                # url = config_defaults.pdbe_download + "{}.cif".format(pdbid)
                url = config.http_pdbe + ("entry-files/download/"
                                          "{}_updated.cif".format(identifier))
        r = fetch_from_url_or_retry(url, stream=True)
        if r.ok:
            with open(outputfile, 'wb') as f:
                r.raw.decode_content = True
                shutil.copyfileobj(r.raw, f)

            if filename.endswith('.gz'):
                with gzip.open(outputfile, 'rb') as infile, \
                        open(outputfile.replace('.gz', ''), 'wb') as outfile:
                    shutil.copyfileobj(infile, outfile)
                    os.remove(outputfile)
    return


def download_sifts_from_ebi(identifier, override=False):
    """
    Downloads a SIFTS xml from the EBI FTP to the filesystem.

    :param identifier: (str) PDB ID
    :param override: (boolean)
    :return: (side effects)
    """

    filename = "{}.xml.gz".format(identifier)

    outputfile = "{}{}{}".format(config.db_root, config.db_sifts, filename)
    os.makedirs("{}{}".format(config.db_root, config.db_sifts), exist_ok=True)
    if not os.path.exists(outputfile.rstrip('.gz')) or override:

        url = config.ftp_sifts + "{}.xml.gz".format(identifier)
        try:
            urllib.urlretrieve(url, outputfile)
        except AttributeError:
            # python 3.5+
            urllib.request.urlretrieve(url, outputfile)
        except IOError as e:
            message = 'Unable to retrieve {} for {}'.format(url, str(e))
            flash(message)

        if filename.endswith('.gz'):
            with gzip.open(outputfile, 'rb') as infile, \
                    open(outputfile.replace('.gz', ''), 'wb') as outfile:
                shutil.copyfileobj(infile, outfile)
                os.remove(outputfile)
    return


def download_data_from_uniprot(identifier, file_format="fasta", override=False):
    """
    Downloads a UniProt fasta, gff or txt to the filesystem.

    :param identifier: (str) UniProt ID
    :param file_format: (str) endpoint
    :param override: (boolean)
    :return: (side effects)
    """

    file_format = file_format.lstrip('.')
    if file_format in ['txt', 'fasta', 'gff']:
        filename = "{}.{}".format(identifier, file_format)

        outputfile = "{}{}{}".format(config.db_root, config.db_uniprot, filename)
        os.makedirs("{}{}".format(config.db_root, config.db_uniprot), exist_ok=True)
        if not os.path.exists(outputfile) or override:

            url = config.http_uniprot + "{}.{}".format(identifier, file_format)
            r = fetch_from_url_or_retry(url, stream=True)
            if r.ok:
                with open(outputfile, 'wb') as f:
                    r.raw.decode_content = True
                    shutil.copyfileobj(r.raw, f)
    else:
        raise ValueError("File format {} is not currently implemented..."
                         "".format(file_format))
    return


def download_alignment_from_cath(identifier, max_sequences=200, override=False):
    """
    Downloads a MSA in fasta format from CATH to the filesystem.

    :param identifier: (str) CATH ID (<Superfamily>_<Funfam>)
    :param max_sequences: (str) Maximum number of sequences (default = 200)
    :param override: (boolean)
    :return: (side effects)
    """

    if '_' in identifier:
        filename = "{}.fasta".format(identifier)
        superfamily, funfam = identifier.split('_')[0], identifier.split('_')[1]
        outputfile = "{}{}{}".format(config.db_root, config.db_cath, filename)
        os.makedirs("{}{}".format(config.db_root, config.db_cath), exist_ok=True)
        if not os.path.exists(outputfile) or override:

            url = config.http_cath + "superfamily/{}/funfam/{}/" \
                                     "files/seed_alignment.fasta" \
                                     "?max_sequences={}".format(superfamily, funfam,
                                                                max_sequences)
            r = fetch_from_url_or_retry(url, stream=True)
            if r.ok:
                with open(outputfile, 'wb') as f:
                    r.raw.decode_content = True
                    shutil.copyfileobj(r.raw, f)
    else:
        raise ValueError("Expected CATH  ID but got {}..."
                         "".format(identifier))
    return


def download_alignment_from_pfam(identifier, alignment_size="seed",
                                 override=False):
    """
    Downloads a MSA in Stockholm format from Pfam to the filesystem.

    :param identifier: (str) PFam ID
    :param alignment_size: (str) either "seed" or "full"
    :param override: (boolean)
    :return: (side effects)
    """

    filename = "{}.sth".format(identifier)

    outputfile = "{}{}{}".format(config.db_root, config.db_pfam, filename)
    os.makedirs("{}{}".format(config.db_root, config.db_pfam), exist_ok=True)
    if not os.path.exists(outputfile) or override:

        url = config.http_pfam + "family/{}/alignment/{}".format(identifier,
                                                                 alignment_size)
        r = fetch_from_url_or_retry(url, stream=True)
        if r.ok:
            with open(outputfile, 'wb') as f:
                r.raw.decode_content = True
                shutil.copyfileobj(r.raw, f)
    return


if __name__ == '__main__':
    pass
