# Doradocker
Docker image for Fungal ITS Barcoding via latest Oxford Nanopore Tech and Dorado basecaller

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
4. Create a folder on the host (main) OS that contains your raw ONT data, in fast5/pod5 format, along with the Index.txt and Primer.txt file for demux
5. Run an instance (container) of the image. On this first run, must specify to use all GPUs for CUDA access. Change the bind mount (-v) to the data folder from previous step
`docker run -it --gpus all -v ~/Desktop/data:/data -rm --name doradocker-live doradocker`  
When finished, use `exit` command to quit the Docker container  
6. This will start a bash terminal in a transient container. The container, and all other files outside of the data folder, will be removed when finished. Rerun the command to create a new instance on-demand.
7. Follow the Doradocker_Run_Commands doc for the complete pipeline instructions


## To Do:
1. Maybe break up dockerfile for better caching
2. Add support for Singularity?  
  `singularity build [IMAGE NAME].sif docker-daemon://[IMAGE NAME]:latest`
3. Create alternative image for CUDA-less systems? Could base on miniconda image.
4. Set up as non-root user?
5. Reduce docker image file size. Add anything dockerignore?
6. Incorporate cleanup commands for intermediate files (bam and fastq)
