# Next-Gen Dorado Basecalling and ITS Barcoding Pipeline
# Version 0.1
# 2024-04-15
# Z. Geurin

# Update this to latest runtime when needed; Current 2024-04-15
# https://hub.docker.com/r/nvidia/cuda
FROM nvidia/cuda:12.4.1-runtime-ubuntu22.04

# set docker shell to bash
SHELL ["/bin/bash", "-c"]

# Install any additional system dependencies
RUN apt-get update && apt-get install -y \
    build-essential wget git pigz nano

# create non-root user. -m creates home directory. -s sets shell to bash
# RUN useradd -ms /bin/bash doradocker 

# Install conda (miniforge)
RUN wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" && \
    bash Miniforge3-$(uname)-$(uname -m).sh -b -p /miniconda3/ && \
    rm Miniforge3-$(uname)-$(uname -m).sh
ENV PATH="${PATH}:/miniconda3/bin"
RUN conda config --add channels bioconda && \
    conda update -n base conda && \
    conda config --set auto_activate_base false && \
    conda config --set solver libmamba && \
    conda init bash 

# symlink python3 to python
RUN ln -s /usr/bin/python3 /usr/bin/python

# Set up various conda env
RUN conda create -n samtools -c bioconda "samtools>=1.19"

RUN conda create -n porechop -c bioconda porechop

RUN conda create -n pycoqc -c bioconda pycoqc plotly h5py tqdm cython pyabpoa

# Get latest NGSpeciesID from GitHub (conda releases are outdated)
# Note: this pip install is local to the conda NGSpeciesID env
RUN conda create -n NGSpeciesID python=3.6 pip medaka=0.11.5 openblas=0.3.3 spoa racon minimap2 && \
    conda run -n NGSpeciesID pip install NGSpeciesID


# Install minibar demultiplexing software
RUN wget https://raw.githubusercontent.com/calacademy-research/minibar/master/minibar.py && \
    chmod 775 minibar.py

# Install Dorado from GitHub
RUN conda create -n dorado && \
    conda run -n dorado pip install pod5 && \
    wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.6.0-linux-x64.tar.gz && \
    tar -xvf dorado-0.6.0-linux-x64.tar.gz && \
    rm dorado-0.6.0-linux-x64.tar.gz && \
    mv dorado-0.6.0-linux-x64 /miniconda3/envs/dorado/bin/

# for some reason, mv and cp only move the subcontents of dorado-0.6.0-linux-x64 (the bin and lib subfolders) to /miniconda3/envs/dorado/bin/
ENV PATH="${PATH}:/miniconda3/envs/dorado/bin/bin"


# Set up working directories within the container, copy necessary scripts and supporting files
