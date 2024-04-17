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
4. Run an instance (container) of the image. On this first run, must specify to use all GPUs for CUDA access.  
`docker run -it --gpus all --name doradocker-live doradocker`  
When finished, use `exit` command to quit the Docker container  
5. To re-open the Docker container, first start the container, then open an interactive shell  
`docker start doradocker-live`  
`docker exec -it doradocker-live /bin/bash`  
When finished, use `exit` command to exit the container's shell. Then stop the container:  
`docker stop doradocker-live`


## To Do:
1. Add volume or bind mounts for access to host data
2. Configure working directories
3. Maybe break up dockerfile for better caching
4. Use transient containers: add `-rm` to the run command
5. Add support for Singularity?  
`singularity build [IMAGE NAME].sif docker-daemon://[IMAGE NAME]:latest`
