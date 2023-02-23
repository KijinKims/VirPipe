import os
import pyfastx
import warnings

class InputError(Exception):
    def __init__(self, message):
        self.message = message
        
def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'

class Validator():
    def __init__(self, args):

        self.args = args
        self.task = args.get('task')

    def validate(self):
        self.validate_common_params()
        self.validate_input_file()

        if self.task == 'all':
            self.validate_all()
        elif self.task == 'qc':
            self.validate_qc()
        elif self.task == 'preprocess':
            self.validate_preprocess()
        elif self.task == 'remove-host':
            self.validate_remove_host()
        elif self.task == 'classify-taxonomy':
            self.validate_classify_taxonomy()
        elif self.task == 'assemble':
            self.validate_assemble()
        elif self.task == 'polish':
            self.validate_polish()
        elif self.task == 'map':
            self.validate_map()
        elif self.task == 'blast':
            self.validate_blast()
        elif self.task == 'zoonotic-rank':
            self.validate_zoonotic_rank()

        return self.args

    def validate_common_params(self):

        # --outdir copied from --prefix, if not given
        if not self.args.get('outdir'):
            self.args['outdir'] = self.args.get('prefix')

        # --config file check
        if self.args.get('config'):
            self.file_exists(self.args.get('config'))

    def validate_input_file(self):

        platform = self.args.get('platform')

        # if --platform given
        if platform:
            # nanopore
            if platform == 'nanopore':
                self.is_fastq(self.args.get('fastq'))
                if self.args.get('fastq2'):
                    raise InputError("--fastq2 is not required for nanopore platform.")

            # illumina
            elif platform == 'illumina':
                self.is_fastq(self.args.get('fastq'))
                if not self.args.get('fastq2'):
                    raise InputError("Only paired input is allowed for illumina platform. --fastq2 is required.")
                self.is_fastq(self.args.get('fastq2'))
        else:
            self.is_fasta(self.args.get('fasta'))

    def file_exists(self, filename):
        if os.path.exists(filename):
            return True
        else:
            raise FileNotFoundError("The file %s does not exists." % filename)

    def dir_exists(self, dirname):
        if os.path.exists(dirname):
            return True
        else:
            raise FileNotFoundError("The directory %s does not exists." % dirname)

    def is_fasta(self, filename):
        pyfastx.Fasta(filename, build_index=False)
        return True

    def is_fastq(self, filename):
        pyfastx.Fastq(filename, build_index=False)
        return True

    def validate_qc(self):
        pass

    def validate_preprocess(self):
        pass

    def validate_remove_host(self):
        if not self.args.get('host_genome'):
            raise InputError("--host-genome is required.  If you don't want to remove host genome, use --skip-remove-host.")
        self.is_fasta(self.args.get('host_genome'))

    def validate_classify_taxonomy(self):
        platform = self.args.get('platform')

        if platform == 'nanopore':
            if self.args.get('centrifuge_db'):
                self.dir_exists(os.path.dirname(self.args.get('centrifuge_db')))
            else:
                raise InputError("--centrifuge-db is required for nanopore platform. If you don't want to classify taxonomy, use --skip-classify-taxonomy.")

        elif platform == 'illumina':
            if self.args.get('kraken2_db'):
                self.dir_exists(self.args.get('kraken2_db'))
            else:
                raise InputError("--kraken2-db is required for illumina platform. If you don't want to classify taxonomy, use --skip-classify-taxonomy.")
                

    def validate_assemble(self):
        platform = platform = self.args.get('platform')

        if self.args.get('assembly_tool'):
            if platform == 'nanopore' and self.args.get('assembly_tool') == 'spades':
                raise InputError("--assembly-tool can be chosen from ['megahit', 'canu', 'flye'] for nanopore platform")
            elif platform == 'illumina' and self.args.get('assembly_tool') != 'spades':
                raise InputError("--assembly-tool can be chosen from ['spades'] for illumina platform")
            
        else:
            platform = self.args.get('platform')

            if platform == 'nanopore':
                warnings.warn("--assembly-tool is not given. Flye will be used as a default.")
                self.args['assembly_tool'] = 'flye'
            else:
                warnings.warn("--assembly-tool is not given. Spades will be used as a default.")
                self.args['assembly_tool'] = 'spades'

    def validate_polish(self):
        self.is_fastq(self.args.get('reads'))

    def validate_map(self):
        if self.args.get('ref'):
            for ref in self.args.get('ref'):
                self.is_fasta(ref)
        if self.args.get('dir_ref'):
            self.dir_exists(self.args.get('dir_ref'))

    def validate_blast(self):
        if not self.args.get('blast_db'):
            raise InputError("--blast-db is required. If you don't want to run blast, use --skip-blast.")
        self.dir_exists(os.path.dirname(self.args.get('blast_db')))

        if not self.args.get('taxonomizr_db'):
            raise InputError("--taxonomizr-db is required. If you don't want to run blast, use --skip-blast.")
        self.file_exists(self.args.get('taxonomizr_db'))
        if is_gz_file(self.args.get('taxonomizr_db')):
            raise InputError("taxonomizr db is in gzipped format. It should be decompressed.")

    def validate_zoonotic_rank(self):
        pass

    def parse_args_from_file(self, filename):
        int_k = ['illumina_min_read_quality',
                 'illumina_min_read_length',
                 'nanopore_min_read_quality',
                 'min_contig_length',
                 'min_blast_aln_len']
        
        float_k = ['min_avg_cov',
                   'kraken2_confidence_threshold']

        with open(filename) as f_in:
            lines = filter(None, (line.rstrip() for line in f_in)) # remove empty line
            lines = filter(line.startswith('#') for line in lines) # remove comment line
            for line in lines:
                arg_key, arg_value = line.split('=')
                k = arg_key.strip()
                v = arg_value.strip()

                if k.startswith('save') or k.startswith('skip') or k.startswith('do') or k =='resume':
                    self.args[k] = bool(v)
                elif k in int_k:
                    self.args[k] = int(v)
                elif k in float_k:
                    self.args[k] = float(v)
                else:
                    self.args[k] = k

    def validate_all(self):
        if self.args.get('args_from_file'):
            self.file_exists(self.args.get('args_from_file'))
            self.parse_args_from_file(self.args.get('args_from_file'))

        # no polish step for illumina input
        if self.args.get('platform') == 'illumina':
            self.args['skip_polish'] = True

        if not self.args.get('skip_qc'):
            self.validate_qc()
        if not self.args.get('skip_preprocess'):
            self.validate_preprocess()
        if not self.args.get('skip_remove_host'):
            self.validate_remove_host()
        if not self.args.get('skip_classify_taxonomy'):
            self.validate_classify_taxonomy()
        if not self.args.get('skip_assemble'):
            self.validate_assemble()
        if not self.args.get('skip_polish'):
            pass
        if not self.args.get('skip_map'):
            self.validate_map()
        if not self.args.get('skip_blast'):
            self.validate_blast()
        if self.args.get('do_zoonotic_rank'):
            self.validate_zoonotic_rank()