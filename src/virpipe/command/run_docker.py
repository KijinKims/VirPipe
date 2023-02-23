import os
import sys

def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'

class DockerRunner():
    def __init__(self, args):
        self.args = args
        self.prefix = args.get('prefix')
        self.image = args.get('image')
        self.volumes = {os.getcwd(): {'bind': '/home/user/run', 'mode': 'rw'}}
        self.bind_input()
        self.bind_database()

    def bind_input(self):
        base_input_dir = "/home/user/input"
        for name, path in self.args.items():
            mount_path = ''
            if name == 'ref':
                for i, path in enumerate(self.args['ref']):
                    mount_path = f"{base_input_dir}/{self.prefix}_ref_{i+1}.fasta"
                    self.volumes[path] = {'bind': mount_path, 'mode': 'ro'}
                continue

            if name == 'fastq':
                if is_gz_file(path):
                    mount_path = f"{base_input_dir}/{self.prefix}.fastq.gz"
                else:
                    mount_path = f"{base_input_dir}/{self.prefix}.fastq"
            elif name == 'fastq2':
                if is_gz_file(path):
                    mount_path = f"{base_input_dir}/{self.prefix}_2.fastq.gz"
                else:
                    mount_path = f"{base_input_dir}/{self.prefix}_2.fastq"
            elif name == 'fasta':
                mount_path = f"{base_input_dir}/{self.prefix}.fasta"
            elif name == 'host_genome':
                mount_path = f"{base_input_dir}/{self.prefix}_host_genome.fasta"
            elif name == 'reads':
                if is_gz_file(path):
                    mount_path = f"{base_input_dir}/{self.prefix}_reads.fastq.gz"
                else:
                    mount_path = f"{base_input_dir}/{self.prefix}_reads.fastq"
            elif name == 'dir_ref':
                mount_path = f"{base_input_dir}/dir-ref"

            if mount_path != '':
                self.volumes[path] = {'bind': mount_path, 'mode': 'ro'}

    def bind_database(self):
        base_database_dir = "/home/user/databases"
        for name, path in self.args.items():
            mount_path = ''
            if name == 'kraken2_db':
                mount_path = f"{base_database_dir}/kraken2-db"
            elif name == 'centrifuge_db':
                mount_path = f"{base_database_dir}/centrifuge-db"
            elif name == 'blast_db':
                mount_path = f"{base_database_dir}/blast-db"
            elif name == 'taxonomizr_db':
                mount_path = f"{base_database_dir}/taxonomizr-db.sql"

            if mount_path != '':
                if name in ['centrifuge_db', 'blast_db']:
                    head, _ = os.path.split(path)
                    self.volumes[head] = {'bind': mount_path, 'mode': 'ro'}
                else:
                    self.volumes[path] = {'bind': mount_path, 'mode': 'ro'}

    def run(self, nextflow_command):
        client = docker.from_env()

        container = client.containers.run(image=self.image,
                              command=nextflow_command,
                              detach=True,
                              volumes=self.volumes,
                              working_dir='/home/user/run')
        
        out = container.logs(stdout=True, stderr=False, stream=True)
        err = container.logs(stdout=False, stderr=True, stream=True)

        for line in out:
            if line.decode().strip():
                print(line.decode().strip())
            
        for line in err:
            if line.decode().strip():
                print(line.decode().strip(), file=sys.stderr)
            
        container.remove()