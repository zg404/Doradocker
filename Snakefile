# Dorado pipeline for ONT ITS Barcoding Data Analysis
# Author: Zach Geurin, 2024-07
# Version: 0.1

configfile: "config.yaml"

# TODO: Add schema validation
# validate(config, "config.schema.yaml")

FAST5_to_POD5 = config["FAST5"]
INPUT_DIR = config["INPUT_DIR"]
SAMPLE_INDEX = config["SAMPLE_INDEX"]
PRIMERS = config["PRIMERS"]
ADAPTER_TRIM = config["ADAPTERS"]
THREADS = config["THREADS"]


if FAST5_to_POD5:
  rule POD5_CONVERT:
      input:
          "/data/fast5/{sample}.fast5",
      output:
          "/data/pod5/{sample}.pod5",
      shell:
          "pod5 convert fast5 {input} --output {output} --one-to-one {INPUT_DIR}",

rule DORADO_BASECALL:
    input:
        "/data/pod5/{sample}.pod5",
    output:
        "/data/bamcalls.bam",
    shell:
        "dorado duplex sup {input} > {output}"

rule SAMTOOLS_DEDUPE:
    input:
        "/data/bamcalls.bam",
    output:
        "/data/combinedcalls.bam",
    shell:
        "conda activate samtools",
        "samtools view -d dx:1 -o ./bamcalls.duplex.bam {input}",
        "samtools view -d dx:0 -o ./bamcalls.simplex.bam {input}",
        "samtools merge -o {output} ./bamcalls.simplex.bam ./bamcalls.duplex.bam",
        "conda deactivate"


rule SAMTOOLS_2FASTQ:
    input:
        "/data/combinedcalls.bam",
    output:
        "/data/combinedcalls.fastq",
    shell:
        "conda activate samtools",
        "samtools view -e 'length(seq)<1500'  -O BAM -o combinedcalls.short.bam {input}",
        "samtools fastq -T dx combinedcalls.short.bam > {output}",
        "conda deactivate"

rule PORECHOP_TRIM:
    input:
        "/data/combinedcalls.fastq",
    output:
        "/data/combinedcalls.chopped.fastq",
    shell:
        "conda activate porechop2",
        "porechop_abi --discard_middle --discard_database --custom_adapters /NGSpeciesID/ONT_v14_Adapters.txt -i {input} -o {output}",
        "conda deactivate"

rule SEQUALI_QC:
    input:
        "/data/bamcalls.bam",
        "combinedcalls.chopped.fastq",
    shell:
        "conda activate sequali",
        "sequali {input}",
        "conda deactivate"

rule HOUSEKEEPING:
    shell:
        "mv /NGSpeciesID/ .",
        "cp ./Index.txt ./NGSpeciesID/",
        "cp combinedcalls.chopped.fastq ./NGSpeciesID/",
        "cd ./NGSpeciesID"

rule MINIBAR_DEMUX:
    input:
        "/data/NGSpeciesID/combinedcalls.chopped.fastq",
    output:
        "/data/NGSpeciesID/{sample}.fastq",
    shell:
        "python3 ./minibar.py -F Index.txt combinedcalls.chopped.fastq",
        "rm combinedcalls.chopped.fastq",
        "rm sample_unk.fastq",
        "mv sample_Multiple_Matches.fastq sample_Multiple_Matches.fq"

rule QAQC_FILTER:
    input: 
        "/data/NGSpeciesID/{sample}.fastq",
    output:
        "/data/NGSpeciesID/Length_Filtered_Fastq/{sample}_length-filtered.fastq",
    shell:
        "python FASTQ_IQR_length_filter.py"

rule NGSPECIESID:
    input:
        "/data/NGSpeciesID/Length_Filtered_Fastq/{sample}_length-filtered.fastq",
    output:
        "/data/{sample}/",
    shell:
        "conda activate NGSpeciesID",
        "NGSpeciesID --ont --consensus --t 12 --sample_size 500 --medaka --primer_file primers.txt --fastq {input} --outfolder {output}",
        "conda deactivate"

rule SUMMARIZE:
    shell:
        "summarize.py /data/NGSpeciesID/"

