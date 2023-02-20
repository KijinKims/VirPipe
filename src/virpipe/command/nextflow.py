import os

def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'

class CommandConverter():
        def __init__(self, args):

            self.comp_list = []
            
            self.args = args
            self.prefix = self.args['prefix']
            base_script_dir = "/home/user/nextflow_scripts"
            self.comp_list.append("nextflow run")
            script = f"{base_script_dir}/{self.args['task']}.nf"
            config = f"-c {base_script_dir}/{self.args['task']}.config"
            self.comp_list.append(script)
            self.comp_list.append(config)

            for k, v in args.items():
                if k == 'image':
                    continue
                elif k.startswith('save') or k.startswith('skip') or k.startswith('do'):
                     if v:
                        self.comp_list.append(f"--{k}")
                elif k == 'resume':
                    if v:
                        self.comp_list.append(f"-{k}")
                elif 'db' in k:
                    if k in ['centrifuge_db', 'blast_db']:
                        _, tail = os.path.split(v)
                        self.comp_list.append(self.bind_database(k, v)+f"/{tail}")
                    else:
                        self.comp_list.append(self.bind_database(k, v))
                elif k in ['fastq', 'fastq2', 'fasta', 'host_genome', 'reads', 'ref', 'dir-ref']:
                    self.comp_list.append(self.bind_input(k, v))
                else:
                    self.comp_list.append(f"--{k} {v}")

        def bind_database(self, name, path):
            base_database_dir = "/home/user/databases"
            if name == 'kraken2_db':
                return f"--{name} {base_database_dir}/kraken2-db"
            elif name == 'centrifuge_db':
                return f"--{name} {base_database_dir}/centrifuge-db"
            elif name == 'blast_db':
                return f"--{name} {base_database_dir}/blast-db"
            elif name == 'taxonomizr_db':
                return f"--{name} {base_database_dir}/taxonomizr-db.sql"

        def bind_input(self, name, path):
            base_run_dir = "/home/user/input"
            if name == 'fastq':
                if is_gz_file(path):
                    return f"--{name} {base_run_dir}/{self.prefix}.fastq.gz"
                else:
                    return f"--{name} {base_run_dir}/{self.prefix}.fastq"
            elif name == 'fastq2':
                if is_gz_file(path):
                    return f"--{name} {base_run_dir}/{self.prefix}_2.fastq.gz"
                else:
                    return f"--{name} {base_run_dir}/{self.prefix}_2.fastq"
            elif name == 'fasta':
                return f"--{name} {base_run_dir}/{self.prefix}.fasta"
            elif name == 'host_genome':
                return f"--{name} {base_run_dir}/{self.prefix}_host_genome.fasta"
            elif name == 'reads':
                if is_gz_file(path):
                    return f"--{name} {base_run_dir}/{self.prefix}_reads.fastq.gz"
                else:
                    return f"--{name} {base_run_dir}/{self.prefix}_reads.fastq.gz"
            elif name == 'ref':
                return f"--{name} " + ' '.join([f"{base_run_dir}/{self.prefix}_ref_{i+1}.fasta" for i in range(len(self.args['ref']))])
            elif name == 'dir_ref':
                return f"--{name} {base_run_dir}/dir-ref"

        def generate_nextflow_script_command(self):
            return ' '.join(self.comp_list)