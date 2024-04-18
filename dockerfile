# Next-Gen Dorado Basecalling and ITS Barcoding Pipeline
# Version 0.1
# 2024-04-15
# Z. Geurin

# Update this to latest runtime when needed; Current 2024-04-15
# https://hub.docker.com/r/nvidia/cuda
FROM nvidia/cuda:12.4.1-runtime-ubuntu22.04

# Set docker shell to bash for convenience
SHELL ["/bin/bash", "-c"]

# Install any additional system dependencies
RUN apt-get update && apt-get install -y \
    build-essential wget git pigz nano

# Set up working directories within the container, copy necessary scripts and supporting files
RUN mkdir /data /NGSpeciesID
COPY /scripts/ /NGSpeciesID/
RUN chmod +x /NGSpeciesID/*.py

# Copy in the Conda environment files
COPY /conda_envs/ /conda_envs/


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

# Install pip packages at the system level
RUN pip install biopython  edlib

# Set up various conda env from files
RUN conda env create -f /conda_envs/samtools_environment.yaml

RUN conda env create -f /conda_envs/porechop_environment.yaml

RUN conda env create -f /conda_envs/pycoqc_environment.yaml

RUN conda env create -f /conda_envs/NGSpeciesID_environment.yaml

# Install minibar demultiplexing software
RUN wget https://raw.githubusercontent.com/calacademy-research/minibar/master/minibar.py && \
    chmod +x minibar.py && \
    mv minibar.py /NGSpeciesID/

# Install Dorado from GitHub
RUN conda create -n dorado && \
    conda run -n dorado pip install pod5 && \
    wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.6.0-linux-x64.tar.gz && \
    tar -xvf dorado-0.6.0-linux-x64.tar.gz && \
    rm dorado-0.6.0-linux-x64.tar.gz && \
    mv dorado-0.6.0-linux-x64 /miniconda3/envs/dorado/bin/

# Add dorado location to PATH
# for some reason, mv and cp only move the subcontents (the bin and lib subfolders) to /miniconda3/envs/dorado/bin/
ENV PATH="${PATH}:/miniconda3/envs/dorado/bin/bin"


WORKDIR /data
