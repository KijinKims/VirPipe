import sys
import os
import argparse
from pathlib import PurePath, Path
import warnings
import gzip

from Bio import SeqIO

from virpipe.config import *

def is_existing_nonempty_file(filename):
    if not os.path.exists(filename):
        return False
    elif os.stat(filename).st_size == 0:
        return False
    else:
        return True

def openfile(filename):
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')

def is_fasta(filename):
    with openfile(filename) as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)

def is_fastq(filename):
    with openfile(filename) as handle:
        fastq = SeqIO.parse(handle, "fastq")
        try: return any(fastq)
        except Exception:
            return False

def valid_file(filename):
    if not is_existing_nonempty_file(filename):
        return False, "The file %s does not exists or is empty!" % filename
    else:
        return True, ''

def valid_fastq(filename):
    if not is_existing_nonempty_file(filename):
        return False, "The file %s does not exists or is empty!" % filename
    else:
        if not is_fastq(filename):
            return False, "The file %s is invalid fastq!" % filename
        else:
            return True, ''

def valid_fasta(filename):
    if not is_existing_nonempty_file(filename):
        return False, "The file %s does not exists or is empty!" % filename
    else:
        if not is_fasta(filename):
            return False, "The file %s is invalid fasta!" % filename
        else:
            return True, ''

def valid_fasta_list(filename):
    if not is_existing_nonempty_file(filename):
        return False, "The file %s does not exists or is empty!" % filename

    with open(filename, 'r') as f:
        for line in f.readlines():
            stripped_line = line.strip()
            tf, _ = valid_fasta(stripped_line)
            if not tf:
                return False, "The file %s contains the line %s that is invalid fasta!" % (filename, stripped_line)
    return True, ''

def valid_dir(dirname):
    if not os.path.exists(dirname):
        return False, "The directory %s does not exists or is empty!" % dirname
    else:
        return True, ''

def parser_path_check(parser, func, path):
    
    tf, message = func(path)

    if tf:
        return path
    else:
        parser.error(message)

class Parser:
    def __init__(self, argv):

        nxf_script_dir = str(PurePath(os.path.dirname(os.path.realpath(__file__)), "modules"))

        if os.environ.get('VP_DB'):
            pkgs_dir = os.environ.get('VP_DB')
        else:
            pkgs_dir = str(PurePath(Path.home(),"VP_DB"))
            os.environ['VP_DB'] = pkgs_dir
        

        parser = argparse.ArgumentParser(prog='virpipe', description='%(prog)s is a command line program for detection and analysis of virus sequence from sequencing read input.')
        subparsers = parser.add_subparsers(dest='task', required=True, title='tasks', description='valid tasks')

        shared_parser = argparse.ArgumentParser(add_help=False, argument_default=argparse.SUPPRESS)
        shared_parser.add_argument('--prefix', '-p', nargs='*')
        shared_parser.add_argument('--file-input', '-f', nargs='?', type=lambda x: parser_path_check(shared_parser, valid_file, x))
        shared_parser.add_argument('--outdir', '-o', nargs='*')
        shared_parser.add_argument('--resume', '-r', action='store_true', default=False)
        shared_parser.add_argument('--modules-dir', nargs='?', default=nxf_script_dir)
        shared_parser.add_argument('--config', '-c', nargs='?')
        shared_parser.add_argument('--with-report', nargs='?')
        shared_parser.add_argument('--with-trace', nargs='?')
        shared_parser.add_argument('--with-timeline', nargs='?')

        platform_parser = argparse.ArgumentParser(add_help=False, argument_default=argparse.SUPPRESS)
        platform_parser.add_argument('--platform', nargs='?', choices=['illumina', 'nanopore'], required=True)

        input_parser = argparse.ArgumentParser(add_help=False, argument_default=argparse.SUPPRESS)
        input_parser.add_argument('--fastq', nargs='*', type=lambda x: parser_path_check(input_parser, valid_file, x))
        input_parser.add_argument('--fastq2', nargs='*', type=lambda x: parser_path_check(input_parser, valid_fastq, x))

        postassembly_input_parser = argparse.ArgumentParser(add_help=False, argument_default=argparse.SUPPRESS)
        postassembly_input_parser.add_argument('--fasta', nargs='*', type=lambda x: parser_path_check(postassembly_input_parser, valid_fasta, x))

        qc_parser = subparsers.add_parser('qc', parents=[shared_parser, platform_parser, input_parser], argument_default=argparse.SUPPRESS)
        
        preprocess_parser = subparsers.add_parser('preprocess', parents=[shared_parser, platform_parser, input_parser], argument_default=argparse.SUPPRESS)
        preprocess_parser.add_argument('--illumina-min-read-quality', '-iq', nargs='?', type=int)
        preprocess_parser.add_argument('--nanopore-min-read-quality', '-nq', nargs='?', type=int)
        preprocess_parser.add_argument('--nanopore-min-read-length', '-nl', nargs='?', type=int)
        
        remove_host_parser = subparsers.add_parser('remove-host', parents=[shared_parser, platform_parser, input_parser], argument_default=argparse.SUPPRESS)
        remove_host_parser.add_argument('--host-genome', nargs='*', type=lambda x: parser_path_check(remove_host_parser, valid_fasta, x), required=True)

        map_parser = subparsers.add_parser('map', parents=[shared_parser, platform_parser, input_parser], argument_default=argparse.SUPPRESS)
        map_parser.add_argument('--ref', nargs='*', type=lambda x: parser_path_check(map_parser, valid_fasta, x))
        map_parser.add_argument('--dir-ref', '-dr', nargs='?', type=lambda x: parser_path_check(map_parser, valid_dir, x))
        map_parser.add_argument('--min-map-out-avg-dep', '-dep', nargs='?', type=float)

        taxclassify_parser = subparsers.add_parser('taxclassify', parents=[shared_parser, platform_parser, input_parser], argument_default=argparse.SUPPRESS)
        taxclassify_parser.add_argument('--kraken2-db', nargs='?', type=lambda x: parser_path_check(taxclassify_parser, valid_dir, x))

        assembly_parser = subparsers.add_parser('assembly', parents=[shared_parser, platform_parser, input_parser], argument_default=argparse.SUPPRESS)
        assembly_parser.add_argument('--tool', '-t', nargs='?', choices=['spades', 'megahit', 'canu', 'flye'])
        assembly_parser.add_argument('--min-contig-length', '-l', nargs='?', type=float)
        
        polish_parser = subparsers.add_parser('polish', parents=[shared_parser, postassembly_input_parser], argument_default=argparse.SUPPRESS)
        polish_parser.add_argument('--reads', nargs='*', type=lambda x: parser_path_check(polish_parser, valid_fastq, x))

        blast_parser = subparsers.add_parser('blast', parents=[shared_parser, postassembly_input_parser], argument_default=argparse.SUPPRESS)
        blast_parser.add_argument('--tool', '-t', nargs='*', choices=['blastn' ,'megablast'], default=['blastn' ,'megablast'])
        blast_parser.add_argument('--blast-db-dir', '-d', nargs='?', type=lambda x: parser_path_check(blast_parser, valid_dir, x))
        blast_parser.add_argument('--blast-db-name', '-n', nargs='?')
        blast_parser.add_argument('--min-blast-aln-len', nargs='?', type=float)
        blast_parser.add_argument('--taxonomizr-db', nargs='?', type=lambda x: parser_path_check(blast_parser, valid_file, x))

        zoonosis_parser = subparsers.add_parser('zoonosis', parents=[shared_parser, postassembly_input_parser], argument_default=argparse.SUPPRESS)

        all_parser = subparsers.add_parser('all', parents=[shared_parser, platform_parser, input_parser], argument_default=argparse.SUPPRESS)
        all_parser.add_argument('--host-genome', nargs='*', type=lambda x: parser_path_check(all_parser, valid_fasta, x))
        all_parser.add_argument('--assembly-tool', nargs='?', choices=['spades', 'megahit', 'canu', 'flye'])
        all_parser.add_argument('--skip-qc', action='store_true', default=False)
        all_parser.add_argument('--skip-preprocess', action='store_true', default=False)
        all_parser.add_argument('--skip-remove-host', action='store_true', default=False)
        all_parser.add_argument('--skip-map', action='store_true', default=False)
        all_parser.add_argument('--skip-taxclassify', action='store_true', default=False)
        all_parser.add_argument('--skip-assembly', action='store_true', default=False)
        all_parser.add_argument('--skip-blast', action='store_true', default=False)
        all_parser.add_argument('--skip-zoonosis', action='store_true', default=False)
    
        consensus_parser = subparsers.add_parser('consensus', parents=[shared_parser, platform_parser, input_parser], argument_default=argparse.SUPPRESS)
        consensus_parser.add_argument('--ref', nargs='*', type=lambda x: parser_path_check(consensus_parser, valid_fasta, x), required=True)
        consensus_parser.add_argument('--low-cov-threshold', nargs='?', type=int)
        consensus_parser.add_argument('--variant-quality-threshold', nargs='?', type=int)
        consensus_parser.add_argument('--variant-depth-threshold', nargs='?', type=int)

        self.args = parser.parse_args(argv)

    def get_args(self):
        return vars(self.args)