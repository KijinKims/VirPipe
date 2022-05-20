sample_specific_args = [
                        'prefix',
                        'outdir',
                        'x',
                        'x2',
                        'y',
                        'host_genome',
                        'reads',
                        ]

file_args = [
            'x',
            'x2',
            'y',
            'host_genome',
            'reads',
            ]

not_forwarded_to_nxf_args = [
                            'nextflow_binary', 
                            'nextflow_modules_dir',
                            'include', 
                            'exclude', 
                            'task', 
                            'subtask', 
                            'program', 
                            'analysis', 
                            ]

download_db_list = [
                    'kraken2-standard', 
                    'kraken2-viral',
                    'refseq-viral',
                    'rvdb', 
                    'rvdb-prot', 
                    'taxdump',
                    'accession2taxid',
                    ]