import argparse
import os

def check_positive_int(value):
    try:
        value = int(value)
        if value <= 0:
            raise argparse.ArgumentTypeError("{} is not a positive integer".format(value))
    except ValueError:
        raise Exception("{} is not an integer".format(value))
    return value

def check_positive_float(value):
    try:
        value = float(value)
        if value <= 0.0:
            raise argparse.ArgumentTypeError("{} is not a positive float".format(value))
    except ValueError:
        raise Exception("{} is not an float".format(value))
    return value

class Parser:
    def __init__(self, argv):
        
        parser = argparse.ArgumentParser(prog='virpipe', description='%(prog)s is a command line program for identification of virus sequence from sequencing reads.')
        subparsers = parser.add_subparsers(dest='task', required=True, title='tasks', description='valid tasks')

        common_parser = argparse.ArgumentParser(add_help=False, argument_default=argparse.SUPPRESS)
        common_parser.add_argument('--image', '-i', required=True)
        common_parser.add_argument('--prefix', '-p', required=True)
        common_parser.add_argument('--outdir', '-o')
        common_parser.add_argument('--resume', '-r', action='store_true', default=False)

        platform_parser = argparse.ArgumentParser(add_help=False, argument_default=argparse.SUPPRESS)
        platform_parser.add_argument('--platform', nargs='?', choices=['illumina', 'nanopore'], required=True)

        fastq_input_parser = argparse.ArgumentParser(add_help=False, argument_default=argparse.SUPPRESS)
        fastq_input_parser.add_argument('--fastq', required=True, type=os.path.abspath)
        fastq_input_parser.add_argument('--fastq2', type=os.path.abspath)

        fasta_input_parser = argparse.ArgumentParser(add_help=False, argument_default=argparse.SUPPRESS)
        fasta_input_parser.add_argument('--fasta', required=True, type=os.path.abspath)

        qc_parser = subparsers.add_parser('qc', parents=[common_parser, platform_parser, fastq_input_parser], argument_default=argparse.SUPPRESS)
        qc_parser.add_argument('--nanopore-min-read-quality', nargs='?', type=check_positive_int, default=8)
        qc_parser.add_argument('--nanopore-min-read-length', nargs='?', type=check_positive_int, default=200)

        preprocess_parser = subparsers.add_parser('preprocess', parents=[common_parser, platform_parser, fastq_input_parser], argument_default=argparse.SUPPRESS)
        preprocess_parser.add_argument('--illumina-min-read-quality', nargs='?', type=check_positive_int, default=8)
        preprocess_parser.add_argument('--nanopore-min-read-quality', nargs='?', type=check_positive_int, default=8)
        preprocess_parser.add_argument('--nanopore-min-read-length', nargs='?', type=check_positive_int, default=200)
        
        remove_host_parser = subparsers.add_parser('remove-host', parents=[common_parser, platform_parser, fastq_input_parser], argument_default=argparse.SUPPRESS)
        remove_host_parser.add_argument('--host-genome', nargs='?', required=True, type=os.path.abspath)

        map_parser = subparsers.add_parser('map', parents=[common_parser, platform_parser, fastq_input_parser], argument_default=argparse.SUPPRESS)
        map_parser.add_argument('--save-map-bam', action='store_true', default=False)
        map_parser.add_argument('--ref', nargs='*', type=os.path.abspath)
        map_parser.add_argument('--dir-ref', nargs='?', type=os.path.abspath)
        map_parser.add_argument('--min-avg-cov', nargs='?', type=check_positive_float, default=1.0)

        classify_taxonomy_parser = subparsers.add_parser('classify-taxonomy', parents=[common_parser, platform_parser, fastq_input_parser], argument_default=argparse.SUPPRESS)
        classify_taxonomy_parser.add_argument('--kraken2-confidence-threshold', nargs='?', type=check_positive_float, default=0.1)
        classify_taxonomy_parser.add_argument('--kraken2-db', nargs='?', type=os.path.abspath)
        classify_taxonomy_parser.add_argument('--centrifuge-db', nargs='?', type=os.path.abspath)

        assemble_parser = subparsers.add_parser('assemble', parents=[common_parser, platform_parser, fastq_input_parser], argument_default=argparse.SUPPRESS)
        assemble_parser.add_argument('--assembly-tool', nargs='?', choices=['spades', 'megahit', 'canu', 'flye'])
        assemble_parser.add_argument('--min-contig-length', nargs='?', type=check_positive_int, default=600)

        polish_parser = subparsers.add_parser('polish', parents=[common_parser, fasta_input_parser], argument_default=argparse.SUPPRESS)
        polish_parser.add_argument('--reads', required=True, type=os.path.abspath)

        blast_parser = subparsers.add_parser('blast', parents=[common_parser, fasta_input_parser], argument_default=argparse.SUPPRESS)
        blast_parser.add_argument('--blast-db', nargs='?', required=True, type=os.path.abspath)
        blast_parser.add_argument('--min-evalue', nargs='?', default='1.0e-5')
        blast_parser.add_argument('--min-blast-aln-len', nargs='?', type=check_positive_int, default=100)
        blast_parser.add_argument('--taxonomizr-db', nargs='?', required=True, type=os.path.abspath)

        zoonotic_rank_parser = subparsers.add_parser('zoonotic-rank', parents=[common_parser, fasta_input_parser], argument_default=argparse.SUPPRESS)

        all_parser = subparsers.add_parser('all', parents=[common_parser, platform_parser, fastq_input_parser], argument_default=argparse.SUPPRESS)
        all_parser.add_argument('--save-preprocessed', action='store_true', default=False)
        all_parser.add_argument('--illumina-min-read-quality', nargs='?', type=check_positive_int, default=8)
        all_parser.add_argument('--nanopore-min-read-quality', nargs='?', type=check_positive_int, default=8)
        all_parser.add_argument('--nanopore-min-read-length', nargs='?', type=check_positive_int, default=200)
        all_parser.add_argument('--save-host-removed', action='store_true', default=False)
        all_parser.add_argument('--host-genome', nargs='?', type=os.path.abspath)
        all_parser.add_argument('--save-map-bam', action='store_true', default=False)
        all_parser.add_argument('--ref', nargs='*', type=os.path.abspath)
        all_parser.add_argument('--dir-ref', nargs='?', type=os.path.abspath)
        all_parser.add_argument('--min-avg-cov', nargs='?', type=check_positive_float, default=1.0)
        all_parser.add_argument('--kraken2-confidence-threshold', nargs='?', type=check_positive_float, default=0.1)
        all_parser.add_argument('--kraken2-db', nargs='?', type=os.path.abspath)
        all_parser.add_argument('--centrifuge-db', nargs='?', type=os.path.abspath)
        all_parser.add_argument('--assembly-tool', nargs='?', choices=['spades', 'megahit', 'canu', 'flye'])
        all_parser.add_argument('--min-contig-length', nargs='?', type=check_positive_int, default=600)
        all_parser.add_argument('--reads', type=os.path.abspath)
        all_parser.add_argument('--blast-db', nargs='?', type=os.path.abspath)
        all_parser.add_argument('--min-evalue', nargs='?', default='1.0e-5')
        all_parser.add_argument('--min-blast-aln-len', nargs='?', type=check_positive_int, default=100)
        all_parser.add_argument('--taxonomizr-db', nargs='?', type=os.path.abspath)

        all_parser.add_argument('--skip-qc', action='store_true', default=False)
        all_parser.add_argument('--skip-preprocess', action='store_true', default=False)
        all_parser.add_argument('--skip-remove-host', action='store_true', default=False)
        all_parser.add_argument('--skip-map', action='store_true', default=False)
        all_parser.add_argument('--skip-classify-taxonomy', action='store_true', default=False)
        all_parser.add_argument('--skip-assemble', action='store_true', default=False)
        all_parser.add_argument('--skip-polish', action='store_true', default=False)
        all_parser.add_argument('--skip-blast', action='store_true', default=False)
        all_parser.add_argument('--do-zoonotic-rank', action='store_true', default=False)

        all_parser.add_argument('--args-from-file', nargs='?', type=os.path.abspath)

        self.args = parser.parse_args(argv)

    def get_args(self):
        return vars(self.args) # return as dictionary