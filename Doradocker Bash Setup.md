```
Next-Gen Dorado Basecalling and ITS Barcoding Pipeline Setup (BASH)
Version 0.2
2024-07-22
Z. Geurin
```
## Verify Nvidia CUDA driver is installed
* If using Linux, get the latest GPU driver: https://www.nvidia.com/Download/index.aspx?lang=en-us
* If using WSL2, follow these instructions to install the CUDA toolkit (NOT GPU Drivers): https://developer.nvidia.com/cuda-downloads?target_os=Linux&target_arch=x86_64&Distribution=WSL-Ubuntu&target_version=2.0&target_type=deb_network
* Run this command to verify CUDA driver is detected and working:
`nvidia-smi`

## Install any additional system dependencies
```bash
apt-get update && apt-get install -y \
    build-essential wget git pigz nano unzip
```
_________________________________________________________________________________

## Get git files
```bash
cd ~
wget https://github.com/zg404/Doradocker/archive/main.zip
unzip main.zip
rm main.zip
```

## Set up working directories within the container, copy necessary scripts and supporting files
```bash
mkdir ~/NGSpeciesID
cp ~/Doradocker-main/scripts/ ~/NGSpeciesID/
chmod +x /NGSpeciesID/*.py
```
## Copy in the Conda environment files
```bash
cp /doradocker/conda_envs/ ~/conda_envs/
```

## Install conda (miniforge)
```bash
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" && \
    bash Miniforge3-$(uname)-$(uname -m).sh -b -p ~/miniconda3/ && \
    rm Miniforge3-$(uname)-$(uname -m).sh
echo export PATH="${PATH}:~/miniconda3/bin"
conda config --add channels bioconda && \
    conda update -n base conda && \
    conda config --set auto_activate_base false && \
    conda config --set solver libmamba && \
    conda init bash 
```
## symlink python3 to python
```bash
ln -s /usr/bin/python3 /usr/bin/python
```

## Install pip packages at the system level
```bash
pip install biopython edlib pandas numpy matplotlib
```

## Set up various conda env from files
```bash
conda env create -f ~/Doradocker-main/conda_envs/samtools_environment.yaml

conda env create -f ~/Doradocker-main/conda_envs/porechop2_environment.yaml

conda env create -f ~/Doradocker-main/conda_envs/sequali_environment.yaml

conda env create -f ~/Doradocker-main/conda_envs/NGSpeciesID_environment.yaml
```

## Install minibar demultiplexing software
```bash
wget https://raw.githubusercontent.com/calacademy-research/minibar/master/minibar.py && \
    chmod +x minibar.py && \
    mv minibar.py ~/NGSpeciesID/
```
## Install Dorado from GitHub; current as of 2024-06-20
```bash
conda create -n dorado && \
    conda run -n dorado pip install pod5 && \
    wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.7.2-linux-x64.tar.gz && \
    tar -xvf dorado-0.7.2-linux-x64.tar.gz && \
    rm dorado-0.7.2-linux-x64.tar.gz && \
    mv dorado-0.7.2-linux-x64 /miniconda3/envs/dorado/bin/
```

## Add dorado location to PATH
```bash
echo export ENV PATH="${PATH}:/miniconda3/envs/dorado/bin/"
```
