# The build-stage image:
FROM continuumio/miniconda3 AS build

# Install prerequisites:
RUN apt-get update && apt-get install -y curl bc libgsl-dev

# Change ownership of conda environments to all:
RUN chmod -R 777 /opt/conda

# Add non-root user:
ARG USER_ID
ARG GROUP_ID

RUN addgroup --gid $GROUP_ID user
RUN adduser --disabled-password --gecos '' --uid $USER_ID --gid $GROUP_ID user
USER user

# Copy configurations for conda:
COPY environment_virpipe.yml .
COPY environment_zoonotic_rank.yml .
COPY environment_centrifuge.yml .

# Install the package for VirPipe:
RUN conda install mamba -n base -c conda-forge
RUN mamba env create -f environment_virpipe.yml
ENV PATH /opt/conda/envs/virpipe/bin:$PATH

# Install r package that cannot be installed with conda:
RUN /bin/bash -c "source activate virpipe"
RUN R -e "install.packages(c(\"MetaComp\",\"taxonomizr\"), repos = \"http://cran.us.r-project.org\")"

# Install independent environment for zoonotic_rank:
RUN mamba env create -f environment_zoonotic_rank.yml
ENV PATH /opt/conda/envs/virpipe_zoonotic_rank/bin:$PATH

# Install r package that cannot be installed with conda:
RUN /bin/bash -c "source activate virpipe_zoonotic_rank"
RUN R -e "install.packages(\"betacal\", repos = \"http://cran.us.r-project.org\")"

# Create independent environment for centrifuge:
RUN mamba env create -f environment_centrifuge.yml

# Configure directory structure:
RUN mkdir /home/user/databases

# Create a folder for input files
RUN mkdir /home/user/input

# Create a folder for running
RUN mkdir /home/user/run

# Install additional tools:
COPY bamcov /home/user/bamcov

# Download krona taxonomy database
#WORKDIR /opt/conda/envs/virpipe/opt/krona/taxonomy
#RUN wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz.md5
#RUN wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz && md5sum -c taxdump.tar.gz.md5
# or Copy it from local
#COPY taxdump.tar.gz /opt/conda/envs/virpipe/opt/krona/taxonomy

RUN /bin/bash -c "source activate virpipe"
WORKDIR /opt/conda/envs/virpipe/opt/krona
RUN /bin/bash ktUpdateTaxonomy.sh
# If taxonomy database is downloaded manually previously
#RUN /bin/bash ktUpdateTaxonomy.sh --only-build

# Copy nextflow scripts:
COPY --chown=$USER_ID:$GROUP_ID nextflow_scripts/ /home/user/nextflow_scripts/
COPY --chown=$USER_ID:$GROUP_ID custom_scripts/ /home/user/custom_scripts/
