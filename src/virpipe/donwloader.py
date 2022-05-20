#!/usr/bin/python3
import gzip
from pathlib import Path
import requests
import re
import sys
import os
from argparse import Namespace
import shutil
import tarfile

from virpipe.arguments_loader import *

def download_url(url, target_dir : Path, decompress=False, untar=False):
    target_dir.mkdir(parents=True, exist_ok=True)
    timeout_sec = 5
    try:
        r = requests.get(url, allow_redirects=True, timeout=timeout_sec, stream=True)
    except requests.exceptions.Timeout:
        print(f"There was no response for {url} made within {timeout_sec} seconds. Try again after a while.", file=sys.stderr)
        return

    if 'Content-Disposition' in r.headers.keys():
        filename = get_filename_from_cd(r.headers.get('Content-Disposition'))
    else:
        filename = url.split("/")[-1]
    
    filename_path = Path(str(target_dir), filename)
    if filename_path.exists():
        print(f"The file {str(filename_path)} already exists. Please remove it if you want to redownload it.", file=sys.stderr)
    else:
        print(f"Downloading {url} into {target_dir}  ...")
        wget.download(url, str(filename_path))

        if decompress:
            with gzip.open(filename_path, 'rb') as f_in:
                with open(os.path.splitext(filename_path)[0], 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(filename_path)

        if untar:
            tar = tarfile.open(os.path.splitext(filename_path)[0], "r:")
            tar.extractall(str(target_dir))
            os.remove(os.path.splitext(filename_path)[0])

def get_filename_from_cd(cd):
    """
    Get filename from content-disposition
    """
    if not cd:
        return None
    fname = re.findall('filename=(.+)', cd)

    if len(fname) == 0:
        return None

    return fname[0]

class DBDownloader:
    def __init__(self, args):
        self.argsloader : ArgsLoader = self.load_args_to_loader(args)

    def load_args_to_loader(self, args : Namespace) -> ArgsLoader:
        argsloader = ArgsLoader()

        for attr in filter(lambda a: not a.startswith('_'), dir(args) ):
            if type(getattr(args, attr)) == list:
                argsloader.add(ListArg(attr, getattr(args, attr)))
            elif type(getattr(args, attr)) == str:
                    argsloader.add(ValueArg(attr, getattr(args, attr)))
            else:
                print("This case cannot happen.", file=sys.stderr)
                exit(1)

        return argsloader

    def download(self):
        target_dir = self.argsloader["dir"]
        
        if self.argsloader.has("db"):
            self.db_download(self.argsloader['db'], target_dir)
        
    def db_download(self, db_list, target_dir):
        
        for db in db_list:
            if db == "accession2taxid":
                a2t_list = ["https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz",
                            "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz",
                            "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz"
                            ]

                for a2t in a2t_list:
                    download_url(a2t, Path(target_dir, "accession2taxid"))

            elif db == "refseq-viral":
                fna_list = ["https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz",
                            "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.1.genomic.fna.gz",
                            "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.3.1.genomic.fna.gz",
                            "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.4.1.genomic.fna.gz",
                            ]

                for fna in fna_list:
                    download_url(fna, Path(target_dir, db), decompress=True)
            else:

                url_d = {
                        'kraken2-standard' : 'https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20210517.tar.gz',
                        'kraken2-viral' : 'https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20210517.tar.gz',
                        'rvdb' : 'https://rvdb.dbi.udel.edu/download/C-RVDBv23.0.fasta.gz',
                        'rvdb-prot' : 'https://rvdb-prot.pasteur.fr/files/U-RVDBv23.0-prot.fasta.xz',
                        'taxdump' : 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
                }

                if url_d[db].endswith('.tar.gz') or url_d[db].endswith('.tgz'):
                    download_url(url_d[db], Path(target_dir, db), decompress=True, untar=True)
                else:
                    download_url(url_d[db], Path(target_dir, db), decompress=True)
