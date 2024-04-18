# Doradocker Run Commands
Follow these commands for the complete pipeline; including Dorado Basecalling, read clean up, Minibar demultiplex, and NGSpeciesID consensus
In most cases, the commands can be copied directly into the Docker container terminal
##### Author: Z. Geurin | Version 0.1, 2024-04-15

## File and Folder Preparation
* Prepare working directory with raw data and accessory files

## Docker Initialization
* Build docker image (if not already performed)
```
docker build -t doradocker .
```
* Create new transient instance (MAIN RUN COMMAND TO START THE CONTAINER; START HERE WHEN PROCESSING NEW DATA)
* NOTE: the "--rm" param removes this instance when finished. This is easier than restarting old instances and prevents clutter
```
docker run -it --gpus all -v ~/Desktop/data:/data --rm --name doradocker-live doradocker
```

# Dorado Basecalling
* Verify GPU and CUDA driver is detected
```
nvidia-smi
```

* Convert Fast5 to Pod5 (if necessary)
```
conda activate dorado
pod5 convert fast5 ./fast5/*.fast5 --output ./pod5s/ --one-to-one ./fast5
conda deactivate
```

* Basecall in duplex mode, super accurate model
```
conda activate dorado
dorado duplex sup ./pod5s > bamcalls.bam
```

### Run Summary QC 
* Generate summary file; should be compatible with QC tools that use Guppy summary file
```
dorado summary bamcalls.bam > dorado_raw_QC.txt
conda deactivate
```

### Samtools Read Handling
* Split out for duplex and simplex reads, respectively
```
conda activate samtools
samtools view -d dx:1 -o ./bamcalls.duplex.bam ./bamcalls.bam
samtools view -d dx:0 -o ./bamcalls.simplex.bam ./bamcalls.bam
```

* Merge back into a single bam file, leaving behind the redundant duplex reads
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

### Porechop Read Trimming
* Trim ONT adapters using Porechop
```
conda activate porechop
porechop --discard_middle -i combinedcalls.fastq -o combinedcalls.clean.fastq
conda deactivate
```

### QC Reports
* Perform qc plots with pycoQC
```
conda activate pycoqc
dorado summary combinedcalls.short.bam > dorado_clean_QC.txt
pycoQC -f dorado_clean_QC.txt -o QCreport.html
conda deactivate
```

### Minibar Demultiplexing
* Move the prepared NGSpeciesID folder into the working directory
* Move/copy the combinedcalls.clean.fastq and Index.txt files to the NGSpeciesID folder
```
mv /NGSpeciesID/ .
cp ./Index.txt ./NGSpeciesID/
cp combinedcalls.clean.fastq ./NGSpeciesID/
cd ./NGSpeciesID
python3 ./minibar.py -F Index.txt combinedcalls.clean.fastq
rm combinedcalls.clean.fastq
rm sample_unk.fastq
mv sample_Multiple_Matches.fastq sample_Multiple_Matches.fq
```

### NGSpeciesID Consensus Sequence Generation
```
conda activate NGSpeciesID
for file in *.fastq; do
bn=`basename $file .fastq`
NGSpeciesID --ont --consensus --sample_size 500 --m 730 --s 400 --medaka --primer_file primers.txt --fastq $file --outfolder ${bn}
done
conda deactivate
```

### Summarize Script
```
python summarize.py .
```