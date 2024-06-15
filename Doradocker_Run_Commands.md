# Doradocker Run Commands
Follow these commands for the complete pipeline; including Dorado Basecalling, read clean up, Minibar demultiplex, and NGSpeciesID consensus
In most cases, the commands can be copied directly into the Docker container terminal
##### Author: Z. Geurin | Version 0.1, 2024-04-15

## File and Folder Preparation
* Prepare working directory with raw data and sample sheet index file:
  * **Raw ONT reads** in fast5/pod5 format, in a **folder labeled `fast5/` or `pod5s/`**
  * **Index.txt** sample sheet serving as the key to sample:index combinations
* This folder will be mounted to Docker container without having to migrate data in or out of the container. The default mount point is `/data/` in the root directory.
* **NOTE**: any files/data outside of the `/data/` folder within the container will be wiped after closing the container. By default, the container terminal starts in the `/data/` folder, and all following pipeline operations take place within this folder. Files can be monitored in real time on the host system (outside the container) as the workflow progresses.



## Docker Initialization
* Create new transient container by using the `docker run` command each time you wish to process run data.
* Change the folder pathway (`~/Desktop/data`) in the `-v` parameter to match your prepared folder on the host filesystem. Do not change the `:data` part of the parameter.
* **NOTE**: the `--rm` argument removes this container instance when finished. This is easier than restarting old instances and prevents clutter from building up.
```
docker run -it --gpus all -v ~/Desktop/data:/data --rm --name doradocker-live doradocker
```
This will start an interactive bash terminal in a new container instance of the Doradocker image. When finished, use `exit` command to quit the Docker container.  

**Reminder again**: the container, and all other files outside of the `/data/` folder, will be removed when finished. If for some reason you perform additional work outside of the data folder, copy the files into the `/data/` folder before exiting the container. You can quickly check the folder on the host system to ensure you have all data needed before exiting.

## Dorado Basecalling
* Verify GPU and CUDA driver are detected. Command should output a table with the CUDA driver stated in the top-right, and your GPU model listed.
```
nvidia-smi
```
Dorado takes Pod5 input instead of the legacy Fast5 format. While it's easy to convert from Fast5 to Pod5, it's best to output Pod5 directly from MinKnow.
* Convert Fast5 to Pod5 (if necessary)
```
conda activate dorado
pod5 convert fast5 ./fast5/*.fast5 --output ./pod5s/ --one-to-one ./fast5
conda deactivate
```

* Basecall in duplex mode, super accurate model
* Duplex basecalling forces adapter trimming as part of the algorithm.
* The resulting bam file contains reads that are tagged with their origin (ie, simplex vs duplex)
```
conda activate dorado
dorado duplex sup ./pod5s > bamcalls.bam
```

### ONT Run Summary QC 
* Generate run summary file; compatible with QC tools that use traditional Guppy summary file
```
dorado summary bamcalls.bam > dorado_raw_QC.txt
conda deactivate
```
* Perform qc plots with pycoQC. PycoQC outputs a single html file with interactive figures for the run summary info.
```
conda activate pycoqc
pycoQC -f dorado_raw_QC.txt -o QCreport.html
conda deactivate
```

## Samtools Read Processing
* By default, the duplex reads also output a duplicate simplex read in the dataset, tagged with "dx:-1".
* Split out duplex "dx:1" and simplex "dx:0" reads, leaving out the duplicate "dx:-1" reads
```
conda activate samtools
samtools view -d dx:1 -o ./bamcalls.duplex.bam ./bamcalls.bam
samtools view -d dx:0 -o ./bamcalls.simplex.bam ./bamcalls.bam
```

* Merge duplex and simplex back into a single bam file
```
samtools merge -o combinedcalls.bam ./bamcalls.simplex.bam ./bamcalls.duplex.bam
```

* Remove extra long reads (not ITS length)
```
samtools view -e 'length(seq)<1500'  -O BAM -o combinedcalls.short.bam combinedcalls.bam
```

* Convert to fastq for downstream compatibility. “-T dx” ensures the duplex tag is preserved in the fastq header
```
samtools fastq -T dx combinedcalls.short.bam > combinedcalls.fastq
conda deactivate
```

## Porechop Read Trimming
* Trim ONT adapters using Porechop
```
conda activate porechop
porechop -i combinedcalls.fastq -o combinedcalls.chopped.fastq
conda deactivate
```

## Minibar Demultiplexing
* Move the prepared NGSpeciesID folder into the working directory
* Move/copy the combinedcalls.chopped.fastq and Index.txt files to the NGSpeciesID folder
```
mv /NGSpeciesID/ .
cp ./Index.txt ./NGSpeciesID/
cp combinedcalls.chopped.fastq ./NGSpeciesID/
cd ./NGSpeciesID
python3 ./minibar.py -F Index.txt combinedcalls.chopped.fastq
rm combinedcalls.chopped.fastq
rm sample_unk.fastq
mv sample_Multiple_Matches.fastq sample_Multiple_Matches.fq
```

## NGSpeciesID Consensus Sequence Generation
* Run NGSpeciesID in a loop on all fastq files (output by minibar) in the working directory
* Be sure to change the thread count based on your CPU (--t 12)
```
conda activate NGSpeciesID
for file in *.fastq; do
bn=`basename $file .fastq`
NGSpeciesID --ont --consensus --t 12 --sample_size 500 --m 730 --s 400 --medaka --primer_file primers.txt --fastq $file --outfolder ${bn}
done
conda deactivate
```

### Summarize Script
* Run a script to prepare the summary folder for upload to MycoMap
```
python summarize.py .
```
