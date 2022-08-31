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

def valid_fastx(filename):
    if not is_existing_nonempty_file(filename):
        return False, "The file %s does not exists or is empty!" % filename
    else:
        if not is_fasta(filename) and not is_fastq(filename):
            return False, "The file %s is invalid fastx!" % filename
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

        nxf_script_dir = str(PurePath(os.path.dirname(os.path.realpath(__file__)), "nextflow_modules"))

        if os.environ.get('VP_DB'):
            pkgs_dir = os.environ.get('VP_DB')
        else:
            pkgs_dir = str(PurePath(Path.home(),"VP_DB"))
            os.environ['VP_DB'] = pkgs_dir
        

        parser = argparse.ArgumentParser(prog='virpipe', description='%(prog)s is a command line program for detection and analysis of Bunyavirus sequence.')
        subparsers = parser.add_subparsers(dest='task', required=True, title='tasks', description='valid tasks')

        shared_parser = argparse.ArgumentParser(add_help=False, argument_default=argparse.SUPPRESS)
        shared_parser.add_argument('--prefix', '-p', nargs='*')
        shared_parser.add_argument('--file_input', '-f', nargs='?', type=lambda x: parser_path_check(shared_parser, valid_file, x))
        shared_parser.add_argument('--outdir', '-o', nargs='*')
        shared_parser.add_argument('--resume', '-r', action='store_true', default=False)
        shared_parser.add_argument('--background', '-bg', action='store_true', default=False)
        shared_parser.add_argument('--nextflow_binary', nargs='?', default='nextflow')
        shared_parser.add_argument('--nextflow_config', '-c', nargs='?')
        shared_parser.add_argument('--nextflow_modules_dir', nargs='?', default=nxf_script_dir)

        platform_parser = argparse.ArgumentParser(add_help=False, argument_default=argparse.SUPPRESS)
        platform_parser.add_argument('--platform', nargs='?', choices=['illumina', 'nanopore'], required=True)

        input_parser = argparse.ArgumentParser(add_help=False, argument_default=argparse.SUPPRESS)
        input_parser.add_argument('-x', nargs='*', type=lambda x: parser_path_check(input_parser, valid_file, x))
        input_parser.add_argument('-x2', nargs='*', type=lambda x: parser_path_check(input_parser, valid_fastq, x))

        post_assembly_input_parser = argparse.ArgumentParser(add_help=False, argument_default=argparse.SUPPRESS)
        post_assembly_input_parser.add_argument('-x', nargs='*', type=lambda x: parser_path_check(post_assembly_input_parser, valid_fasta, x))

        general_input_parser = argparse.ArgumentParser(add_help=False, argument_default=argparse.SUPPRESS)
        general_input_parser.add_argument('-x', nargs='*', type=lambda x: parser_path_check(general_input_parser, valid_file, x))
        
        qc_parser = subparsers.add_parser('qc', parents=[shared_parser, platform_parser, input_parser], argument_default=argparse.SUPPRESS)

        filter_parser = subparsers.add_parser('filter')
        filter_subparsers = filter_parser.add_subparsers(dest='subtask', title='subtasks', description='valid subtasks', required=True)
        
        filter_reads_parser = filter_subparsers.add_parser('reads', parents=[shared_parser, platform_parser, input_parser], argument_default=argparse.SUPPRESS)
        filter_reads_parser.add_argument('--illumina_min_read_quality', '-iq', nargs='?', type=int)
        filter_reads_parser.add_argument('--nanopore_min_read_quality', '-nq', nargs='?', type=int)
        filter_reads_parser.add_argument('--nanopore_min_read_length', '-nl', nargs='?', type=int)
        
        filter_host_parser = filter_subparsers.add_parser('host', parents=[shared_parser, platform_parser, input_parser], argument_default=argparse.SUPPRESS)
        filter_host_parser.add_argument('--host_genome', nargs='*', type=lambda x: parser_path_check(filter_host_parser, valid_fasta, x))

        filter_map_parser = filter_subparsers.add_parser('map', parents=[shared_parser, general_input_parser], argument_default=argparse.SUPPRESS)
        filter_map_parser.add_argument('--min_map_out_avg_cov', '-cov', nargs='?', type=float)

        filter_contigs_parser = filter_subparsers.add_parser('contigs', parents=[shared_parser, post_assembly_input_parser], argument_default=argparse.SUPPRESS)
        filter_contigs_parser.add_argument('--min_contig_length', '-l', nargs='?', type=float)
        
        filter_blast_parser = filter_subparsers.add_parser('blast', parents=[shared_parser, general_input_parser], argument_default=argparse.SUPPRESS)
        filter_blast_parser.add_argument('--min_blast_aln_len', nargs='?', type=float)

        map_parser = subparsers.add_parser('map', parents=[shared_parser, platform_parser, input_parser], argument_default=argparse.SUPPRESS)
        map_parser.add_argument('--ref', nargs='*', type=lambda x: parser_path_check(map_parser, valid_fasta, x))
        map_parser.add_argument('--file_ref', '-fr', nargs='?', type=lambda x: parser_path_check(map_parser, valid_fasta_list, x))
        map_parser.add_argument('--dir_ref', '-dr', nargs='?', type=lambda x: parser_path_check(map_parser, valid_dir, x))

        tax_classify_parser = subparsers.add_parser('tax_classify', parents=[shared_parser, platform_parser, input_parser], argument_default=argparse.SUPPRESS)
        tax_classify_parser.add_argument('--tool', '-t', nargs='*', choices=['kraken2'], default=['kraken2'])
        tax_classify_parser.add_argument('--kraken2_db', nargs='?', type=lambda x: parser_path_check(tax_classify_parser, valid_dir, x))

        assembly_parser = subparsers.add_parser('assembly', parents=[shared_parser, platform_parser, input_parser], argument_default=argparse.SUPPRESS)
        assembly_parser.add_argument('--tool', '-t', nargs='*', choices=['spades', 'megahit', 'canu', 'flye'], required=True)

        polish_parser = subparsers.add_parser('polish', parents=[shared_parser, post_assembly_input_parser], argument_default=argparse.SUPPRESS)
        polish_parser.add_argument('--tool', '-t', nargs='*', choices=['racon', 'medaka'], default=['racon', 'medaka'])
        polish_parser.add_argument('--reads', nargs='*', type=lambda x: parser_path_check(polish_parser, valid_fastx, x))

        post_assembly_parser = subparsers.add_parser('post_assembly')
        post_assembly_subparsers = post_assembly_parser.add_subparsers(dest='subtask', title='subtasks', description='valid subtasks', required=True)

        blast_parser = post_assembly_subparsers.add_parser('blast', parents=[shared_parser, post_assembly_input_parser], argument_default=argparse.SUPPRESS)
        blast_parser.add_argument('--tool', '-t', nargs='*', choices=['blastn' ,'megablast'], default=['blastn' ,'megablast'])
        blast_parser.add_argument('--blast_db_dir', '-d', nargs='?', type=lambda x: parser_path_check(blast_parser, valid_dir, x))
        blast_parser.add_argument('--blast_db_name', '-n', nargs='?')

        zoonosis_parser = post_assembly_subparsers.add_parser('zoonosis', parents=[shared_parser, post_assembly_input_parser], argument_default=argparse.SUPPRESS)
        zoonosis_parser.add_argument('--tool', '-t', nargs='*', choices=['zoonotic_rank'], default='zoonotic_rank')

        post_assembly_all_parser = post_assembly_subparsers.add_parser('all', parents=[shared_parser, post_assembly_input_parser], argument_default=argparse.SUPPRESS)
        post_assembly_all_parser.add_argument('--blast_db_dir', '-d', nargs='?', type=lambda x: parser_path_check(post_assembly_all_parser, valid_dir, x))
        post_assembly_all_parser.add_argument('--blast_db_name', '-n', nargs='?')
        post_assembly_all_parser.add_argument('--taxonomizr_db', nargs='?', type=lambda x: parser_path_check(post_assembly_all_parser, valid_file, x))
        
        report_parser = subparsers.add_parser('report')
        report_subparsers = report_parser.add_subparsers(dest='subtask', title='subtasks', description='valid subtasks', required=True)

        report_blast_parser = report_subparsers.add_parser('blast', parents=[shared_parser, general_input_parser], argument_default=argparse.SUPPRESS)
        report_blast_parser.add_argument('--taxonomizr_db', nargs='?', type=lambda x: parser_path_check(report_blast_parser, valid_file, x))

        end_to_end_parser = subparsers.add_parser('end_to_end', parents=[shared_parser, platform_parser, input_parser], argument_default=argparse.SUPPRESS)
        end_to_end_parser.add_argument('--host_genome', nargs='*', type=lambda x: parser_path_check(end_to_end_parser, valid_fasta, x))

        consensus_parser = subparsers.add_parser('consensus', parents=[shared_parser, platform_parser, input_parser], argument_default=argparse.SUPPRESS)
        consensus_parser.add_argument('--ref', nargs='*', type=lambda x: parser_path_check(consensus_parser, valid_fasta, x))
        consensus_parser.add_argument('--low_cov_threshold', nargs='?', type=int)
        consensus_parser.add_argument('--variant_quality_threshold', nargs='?', type=int)
        consensus_parser.add_argument('--variant_depth_threshold', nargs='?', type=int)

        database_parser = subparsers.add_parser('database')
        database_subparsers = database_parser.add_subparsers(dest='subtask', title='subtasks', description='valid subtasks', required=True)
        
        build_shared_parser = argparse.ArgumentParser(add_help=False, argument_default=argparse.SUPPRESS)
        build_shared_parser.add_argument('--nextflow_binary', nargs='?', default='nextflow')
        build_shared_parser.add_argument('--nextflow_config', '-c', nargs='?')
        build_shared_parser.add_argument('--nextflow_modules_dir', nargs='?', default=nxf_script_dir)

        build_parser = database_subparsers.add_parser('build', argument_default=argparse.SUPPRESS)
        build_parsers = build_parser.add_subparsers(dest='program', title='programs', description='name of program that DB will be used for', required=True)

        build_blast_parser = build_parsers.add_parser('blast', parents = [build_shared_parser], argument_default=argparse.SUPPRESS)
        build_blast_parser.add_argument('--in', '-i', nargs='?', type=lambda x: parser_path_check(build_blast_parser, valid_fasta, x), required=True)
        build_blast_parser.add_argument('--outdir', nargs='?', default='%s/blast/rvdb' % pkgs_dir)
        build_blast_parser.add_argument('--dbname', nargs='?', default='rvdb')
        build_blast_parser.add_argument('--title', nargs='?', required=True)
        build_blast_parser.add_argument('--dbtype', nargs='?', choices=['prot', 'nucl'], required=True)
        build_blast_parser.add_argument('--accession2taxid', type=lambda x: parser_path_check(build_blast_parser, valid_file, x))
        build_blast_parser.add_argument('--parse_seqids', action='store_true')

        build_taxonomizr_parser = build_parsers.add_parser('taxonomizr', parents = [build_shared_parser], argument_default=argparse.SUPPRESS)
        build_taxonomizr_parser.add_argument('--outdir', '-o', nargs='?', default='%s/taxonomizr'% pkgs_dir)
        build_taxonomizr_parser.add_argument('--nucl_gb_a2t', nargs='?', type=lambda x: parser_path_check(build_taxonomizr_parser, valid_file, x))
        build_taxonomizr_parser.add_argument('--nucl_wgs_a2t', nargs='?', type=lambda x: parser_path_check(build_taxonomizr_parser, valid_file, x))
        build_taxonomizr_parser.add_argument('--prot_a2t', nargs='?', type=lambda x: parser_path_check(build_taxonomizr_parser, valid_file, x))
        build_taxonomizr_parser.add_argument('--taxonnodes', nargs='?', type=lambda x: parser_path_check(build_taxonomizr_parser, valid_file, x))
        build_taxonomizr_parser.add_argument('--taxonnames', nargs='?', type=lambda x: parser_path_check(build_taxonomizr_parser, valid_file, x))

        build_kraken2_parser = build_parsers.add_parser('kraken2', parents = [build_shared_parser], argument_default=argparse.SUPPRESS)
        build_kraken2_parser.add_argument('--db', nargs='?', type=str, required=True)
        build_kraken2_parser.add_argument('--library', nargs='?', choices=['archaea', 'bacteria', 'plasmid', 'viral', 'human', 'fungi', 'plant', 'protozoa', 'nr', 'nt', 'UniVec', 'UniVec_Core'])
        build_kraken2_parser.add_argument('--sequence', nargs='*', type=lambda x: parser_path_check(build_kraken2_parser, valid_fasta, x))

        self.args = parser.parse_args(argv)

    def get_args(self):
        return self.args