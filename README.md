# VirPipe

## QuickStart

Requirement: Docker, conda(or miniconda)

```bash
git clone https://github.com/KijinKims/virpipe.git
cd virpipe
conda env create -f environment.yml
conda activate virpipe
pip install src/dist/virpipe-1.0.tar.gz
```

```bash
sh utils/docker_pull_list_of_images.sh docker_images.list
```

```bash
tar -xzvf vp-db.tar.gz && mv VP_DB $HOME/
```

```bash
#test
virpipe end_to_end --platform nanopore -x input.fastq --prefix test
```

