# Doradocker
Docker image for Fungal ITS Barcoding via Oxford Nanopore sequencing, [Dorado basecaller](https://github.com/nanoporetech/dorado), demultiplex by [Minibar](https://github.com/calacademy-research/minibar), and read consensus generation by [NGSpeciesID](https://github.com/ksahlin/NGSpeciesID/).

## Installation
Follow these basic steps. More in depth instructions coming soon.
1. Download Docker for your system: https://docs.docker.com/get-docker/
2. Clone the github with dockerfile plus extras  
`git clone https://github.com/zg404/Doradocker`  
`cd Doradocker`  
or download and extract a zip of the repo:  
`wget https://github.com/zg404/Doradocker/archive/refs/heads/main.zip`  
3. Open a terminal in the git folder and build the docker image using the dockerfile. This may take several minutes. Progress can be monitored in the terminal.    
`docker build -t doradocker .`  
4. The docker image should now be ready to use in creating on-demand containers for processing data. Follow the [Doradocker_Run_Commands](https://github.com/zg404/Doradocker/blob/main/Doradocker_Run_Commands.md) doc for the complete pipeline instructions and ready-to-use commands. 


## To Do:
1. Complete the prep files needed for workflow (eg, Index.txt generation)
2. Maybe break up dockerfile for better caching
3. Create alternative image for CUDA-less systems? Could base on miniconda image.
4. Reduce docker image file size.
5. Incorporate cleanup commands for intermediate files (bam and fastq)
6. Implement workflow scripting (snakemake)
