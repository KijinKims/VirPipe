mkdir databases

# download BLAST database
curl --output blast_ref_viruses_rep_genomes.tar.gz https://ftp.ncbi.nlm.nih.gov/blast/db/v5/v5/ref_viruses_rep_genomes.tar.gz
mkdir -p databases/blast_ref_viruses_rep_genomes && tar -xzvf blast_ref_viruses_rep_genomes.tar.gz -C databases/blast_ref_viruses_rep_genomes

# download Kraken2 database
curl --output k2_viral_20221209.tar.gz https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20221209.tar.gz
mkdir -p databases/k2_viral && tar -zxvf k2_viral_20221209.tar.gz -C databases/k2_viral

# download centrifuge database
curl --output centrifuge_viral_230213.tar.gz https://zenodo.org/api/files/814e2a69-a24c-4810-80a1-06fe3ce3bec9/centrifuge_viral_230213.tar.gz
mkdir -p databases/centrifuge_db && tar -xzvf blast_refseq_viral_230213.tar.gz -C databases/centrifuge_db

# download taxonomizr database
curl --output databases/accessionTaxa.sql.gz https://zenodo.org/api/files/814e2a69-a24c-4810-80a1-06fe3ce3bec9/accessionTaxa.sql.gz
gunzip databases/accessionTaxa.sql.gz