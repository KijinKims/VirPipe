# download BLAST database
curl --output blast_refseq_viral_230213.tar.gz https://zenodo.org/api/files/402dcb71-c95b-40fd-9ef4-6d3906594c61/blast_refseq_viral_230213.tar.gz
mkdir -p blast_refseq_viral && tar -xzvf blast_refseq_viral_230213.tar.gz -C blast_refseq_viral

# download Kraken2 database
curl --output k2_viral_20221209.tar.gz https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20221209.tar.gz
mkdir -p k2_viral && tar -zxvf k2_viral_20221209.tar.gz -C k2_viral

# download centrifuge database
curl --output centrifuge_viral_230213.tar.gz https://zenodo.org/api/files/402dcb71-c95b-40fd-9ef4-6d3906594c61
mkdir -p centrifuge_viral && tar -xzvf blast_refseq_viral_230213.tar.gz -C centrifuge_refseq_viral

# download taxonomizr database
curl --output accessionTaxa.sql.gz https://zenodo.org/api/files/402dcb71-c95b-40fd-9ef4-6d3906594c61/accessionTaxa.sql.gz 
gunzip accessionTaxa.sql.gz