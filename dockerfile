# Next-Gen Dorado Basecalling and ITS Barcoding Pipeline
# Version 0.1
# 2024-04-15
# Z. Geurin

# Update this to latest runtime when needed; Current 2024-04
# https://hub.docker.com/r/nvidia/cuda
FROM nvidia/cuda:12.4.1-runtime-ubuntu22.04

# Install any additional system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    nvidia-container-toolkit 

# Install conda (miniforge)
RUN wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" && \
    bash Miniforge3-$(uname)-$(uname -m).sh && \
    conda config --add channels bioconda && \
    conda update -n base conda && \
    conda config --set auto_activate_base false && \
    conda config --set solver libmamba

# symlink python3 to python
RUN sudo ln -s /usr/bin/python3 /usr/bin/python

# Set up various conda env
RUN mamba create -n samtools -c bioconda samtools>=1.19

RUN mamba create -n porechop -c bioconda porechop

RUN mamba create -n pycoqc -c bioconda pycoqc plotly h5py tqdm cython pyabpoa

# Get latest NGSpeciesID from GitHub (conda releases are outdated)
# Note: this pip install is local to the conda NGSID env
RUN mamba create -n NGSID && \
    mamba activate NGSID && \
    pip install git+https://github.com/ksahlin/NGSpeciesID.git  \
    mamba deactivate

# Set up development environment within the container
WORKDIR /ONT_data
