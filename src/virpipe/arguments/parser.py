import argparse
import os

def check_positive_int(value):
    try:
        value = int(value)
        if value < 0:
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
        parser.add_argument('--version', action='version', version='%(prog)s 1.0')
        subparsers = parser.add_subparsers(dest='task', required=True, title='tasks', description='valid tasks')

        common_parser = argparse.ArgumentParser(add_help=False, argument_default=argparse.SUPPRESS)
        common_parser.add_argument('--image', '-i', required=True, help="The name of the Docker image from which a container for running will be instantiated.")
        common_parser.add_argument('--prefix', '-p', required=True, help="Run name. All outputs are named with this prefix.")
        common_parser.add_argument('--outdir', '-o', help="Name of output directory. All outputs will be located under this directory. If the given path doesn't exist, it will be newly created.")
        common_parser.add_argument('--resume', '-r', action='store_true', default=False, help="If a run has been imperfectly completed, you can resume it from where your run left off by giving this option with the same command line before.")

        platform_parser = argparse.ArgumentParser(add_help=False, argument_default=argparse.SUPPRESS)
        platform_parser.add_argument('--platform', nargs='?', choices=['illumina', 'nanopore'], required=True, help="The platform through which the input file(s) was sequenced.")

        fastq_input_parser = argparse.ArgumentParser(add_help=False, argument_default=argparse.SUPPRESS)
        fastq_input_parser.add_argument('--fastq', required=True, type=os.path.abspath, help="Input fastq file.")
        fastq_input_parser.add_argument('--fastq2', type=os.path.abspath, help="Second input file. This option is only used for Illumina paired-end reads.")

        fasta_input_parser = argparse.ArgumentParser(add_help=False, argument_default=argparse.SUPPRESS)
        fasta_input_parser.add_argument('--fasta', required=True, type=os.path.abspath, help="Input fasta file.")

        qc_parser = subparsers.add_parser('qc', parents=[common_parser, platform_parser, fastq_input_parser], argument_default=argparse.SUPPRESS, help="The task qc is for producing a quality control report for input reads.")
        qc_parser.add_argument('--nanopore-min-read-quality', nargs='?', type=check_positive_int, default=8, help="For nanopore reads, QC statistics are calculated only with the reads with the quality score over this value.")
        qc_parser.add_argument('--nanopore-min-read-length', nargs='?', type=check_positive_int, default=200, help="For nanopore reads, QC statistics are calculated only with the reads with the read length over this value.")

        preprocess_parser = subparsers.add_parser('preprocess', parents=[common_parser, platform_parser, fastq_input_parser], argument_default=argparse.SUPPRESS, help="The task preprocess is for filtering input reads depending on their quality or read length and trimming adapter sequences.")
        preprocess_parser.add_argument('--illumina-min-read-quality', nargs='?', type=check_positive_int, default=8, help="Threshold for read quality filtering for illumina reads.")
        preprocess_parser.add_argument('--nanopore-min-read-quality', nargs='?', type=check_positive_int, default=8, help="Threshold for read quality filtering for nanopore reads.")
        preprocess_parser.add_argument('--nanopore-min-read-length', nargs='?', type=check_positive_int, default=200, help="Threshold for read length filtering for nanopore reads.")
        
        remove_host_parser = subparsers.add_parser('remove-host', parents=[common_parser, platform_parser, fastq_input_parser], argument_default=argparse.SUPPRESS, help="The task remove-host is for getting rid of host-derived reads from the input.")
        remove_host_parser.add_argument('--host-genome', nargs='?', required=True, type=os.path.abspath, help="Host genome sequence in FASTA format.")

        map_parser = subparsers.add_parser('map', parents=[common_parser, platform_parser, fastq_input_parser], argument_default=argparse.SUPPRESS, help="The task map is for mapping input reads onto given reference sequence(s).")
        map_parser.add_argument('--save-map-bam', action='store_true', default=False, help="The flag to indicate whether intermediate BAM files generated during map task will be saved in the output directory.")
        map_parser.add_argument('--ref', nargs='*', type=os.path.abspath, help="Reference sequence(s) in FASTA format. Multiple files can be given by delimiting them with a blank. For example, --ref virus1.fasta virus2.fasta")
        map_parser.add_argument('--dir-ref', nargs='?', type=os.path.abspath, help="Reference sequence(s) given as a directory. All .fasta and .fa files under this directory will be given as references.")
        map_parser.add_argument('--min-avg-cov', nargs='?', type=check_positive_float, default=1.0, help="Threshold for average coverage filtering")

        classify_taxonomy_parser = subparsers.add_parser('classify-taxonomy', parents=[common_parser, platform_parser, fastq_input_parser], argument_default=argparse.SUPPRESS, help="The task classify-taxonomy is for classifying each read into a taxonomy using the k-mer matching strategy.")
        classify_taxonomy_parser.add_argument('--kraken2-confidence-threshold', nargs='?', type=check_positive_float, default=0.1, help="For Kraken2, the classifications with scores less than the threshold will be regarded as unclassified.")
        classify_taxonomy_parser.add_argument('--kraken2-db', nargs='?', type=os.path.abspath, help="Database for Kraken2. It should be the path of a directory containing files such as .k2d")
        classify_taxonomy_parser.add_argument('--centrifuge-db', nargs='?', type=os.path.abspath, help="Database for Centrifuge. it should be the directory + basename of indexes, such as /path/to/dir/centrifuge_db/abv")

        assemble_parser = subparsers.add_parser('assemble', parents=[common_parser, platform_parser, fastq_input_parser], argument_default=argparse.SUPPRESS, help="The task assemble is for constructing genomes from input reads by de novo assembly.")
        assemble_parser.add_argument('--assembly-tool', nargs='?', choices=['spades', 'megahit', 'canu', 'flye'], help="The assembler to be used for assembly.")
        assemble_parser.add_argument('--min-contig-length', nargs='?', type=check_positive_int, default=600, help="The contigs shorter than this value will be discarded.")

        polish_parser = subparsers.add_parser('polish', parents=[common_parser, fasta_input_parser], argument_default=argparse.SUPPRESS, help="The task polish is for correcting errors in contigs assembled from nanopore reads using the reads again.")
        polish_parser.add_argument('--reads', required=True, type=os.path.abspath, help="FASTQ file containing nanopore reads.")

        blast_parser = subparsers.add_parser('blast', parents=[common_parser, fasta_input_parser], argument_default=argparse.SUPPRESS, help="The task blast is for running BLAST with input .fasta file.")
        blast_parser.add_argument('--blast-db', nargs='?', required=True, type=os.path.abspath, help="Database for BLAST. It should be the directory + basename of indexes such as /path/to/dir/blast_db/refseq_viral.")
        blast_parser.add_argument('--min-evalue', nargs='?', default='1.0e-5', help="Threshold for E-value ranged from 0 to 1.")
        blast_parser.add_argument('--min-blast-aln-len', nargs='?', type=check_positive_int, default=100, help="Threshold for alignment length.")
        blast_parser.add_argument('--taxonomizr-db', nargs='?', required=True, type=os.path.abspath, help="Database for taxonomizr.")

        zoonotic_rank_parser = subparsers.add_parser('zoonotic-rank', parents=[common_parser, fasta_input_parser], argument_default=argparse.SUPPRESS, help="The task zoonotic-rank is for assessing zoonotic risk")

        all_parser = subparsers.add_parser('all', parents=[common_parser, platform_parser, fastq_input_parser], argument_default=argparse.SUPPRESS, help="The task all executes all tasks available in VirPipe in a reasonable order.")
        all_parser.add_argument('--save-preprocessed', action='store_true', default=False, help="For preprocess tasks, the flag to indicate whether preprocessed FASTQ will be saved in the output directory.")
        all_parser.add_argument('--illumina-min-read-quality', nargs='?', type=check_positive_int, default=8, help="For preprocess tasks, threshold for read quality filtering for illumina reads.")
        all_parser.add_argument('--nanopore-min-read-quality', nargs='?', type=check_positive_int, default=8, help="For qc and preprocess tasks, threshold for read quality filtering for nanopore reads.")
        all_parser.add_argument('--nanopore-min-read-length', nargs='?', type=check_positive_int, default=200, help="For qc and preprocess tasks, threshold for read length filtering for illumina reads.")
        all_parser.add_argument('--save-host-removed', action='store_true', default=False, help="For remove-host task, the flag to indicate whether host-removed FASTQ will be saved in the output directory.")
        all_parser.add_argument('--host-genome', nargs='?', type=os.path.abspath, help="For remove-host task, Host genome sequence in FASTA format.")
        all_parser.add_argument('--save-map-bam', action='store_true', default=False, help="For map task, the flag to indicate whether intermediate BAM files generated during map task will be saved in the output directory.")
        all_parser.add_argument('--ref', nargs='*', type=os.path.abspath, help="For map task, reference sequence(s) in FASTA format. Multiple files can be given by delimiting them with a blank. For example, --ref virus1.fasta virus2.fasta")
        all_parser.add_argument('--dir-ref', nargs='?', type=os.path.abspath, help="For map task, reference sequence(s) given as a directory. All .fasta and .fa files under this directory will be given as references.")
        all_parser.add_argument('--min-avg-cov', nargs='?', type=check_positive_float, default=1.0, help="For map task, threshold for average coverage filtering")
        all_parser.add_argument('--kraken2-confidence-threshold', nargs='?', type=check_positive_float, default=0.1, help="For Kraken2 in classify-taxonomy task, the classifications with scores less than the threshold will be regarded as unclassified.")
        all_parser.add_argument('--kraken2-db', nargs='?', type=os.path.abspath, help="For classify-taxonomy task, database for Kraken2. It should be the path of a directory containing files such as .k2d")
        all_parser.add_argument('--centrifuge-db', nargs='?', type=os.path.abspath, help="For classify-taxonomy task, Database for Centrifuge. it should be the directory + basename of indexes, such as /path/to/dir/centrifuge_db/abv")
        all_parser.add_argument('--assembly-tool', nargs='?', choices=['spades', 'megahit', 'canu', 'flye'], help="For assemble task, The assembler to be used for assembly.")
        all_parser.add_argument('--min-contig-length', nargs='?', type=check_positive_int, default=600, help="For assemble task, the contigs shorter than this value will be discarded.")
        all_parser.add_argument('--blast-db', nargs='?', type=os.path.abspath, help="For blast task, database for BLAST. It should be the directory + basename of indexes such as /path/to/dir/blast_db/refseq_viral.")
        all_parser.add_argument('--min-evalue', nargs='?', default='1.0e-5', help="For blast task, threshold for E-value ranged from 0 to 1.")
        all_parser.add_argument('--min-blast-aln-len', nargs='?', type=check_positive_int, default=100, help="For blast task, threshold for alignment length.")
        all_parser.add_argument('--taxonomizr-db', nargs='?', type=os.path.abspath, help="For blast task, database for taxonomizr.")

        all_parser.add_argument('--skip-qc', action='store_true', default=False, help="The flag to indicate whether qc task will be skipped.")
        all_parser.add_argument('--skip-preprocess', action='store_true', default=False, help="The flag to indicate whether preprocess task will be skipped.")
        all_parser.add_argument('--skip-remove-host', action='store_true', default=False, help="The flag to indicate whether remove-host task will be skipped.")
        all_parser.add_argument('--skip-map', action='store_true', default=False, help="The flag to indicate whether map task will be skipped.")
        all_parser.add_argument('--skip-classify-taxonomy', action='store_true', default=False, help="The flag to indicate whether classify-taxonomy task will be skipped.")
        all_parser.add_argument('--skip-assemble', action='store_true', default=False, help="The flag to indicate whether assemble task will be skipped.")
        all_parser.add_argument('--skip-polish', action='store_true', default=False, help="The flag to indicate whether polish task will be skipped.")
        all_parser.add_argument('--skip-blast', action='store_true', default=False, help="The flag to indicate whether blast task will be skipped.")
        all_parser.add_argument('--do-zoonotic-rank', action='store_true', default=False, help="The flag to indicate whether zoonotic-rank task will be done.")

        all_parser.add_argument('--args-from-file', nargs='?', type=os.path.abspath)

        self.args = parser.parse_args(argv)

    def get_args(self):
        return vars(self.args) # return as dictionary