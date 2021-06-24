# cfDNA pipeline
This pipeline is intended for processing of cfDNA samples for Quynh.

This pipeline takes gzipped paired-end fastq files and a metadata file specifying references, then outputs bams and relevant information.

## Installation

1. Install [nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation).
   - Make sure you move nextflow to a directory in your PATH variable.
2. Install [docker](https://docs.docker.com/get-docker/).

## Usage
- Example command for fastqs in current directory: ```nextflow run greninger-lab/cfdna_pipeline --METADATA metadata.csv --OUTDIR output/ -resume -with-trace -with-docker ubuntu:18.04 -c ~/nextflow.config -profile cloud```


| Command  | Description |
| ---      | ---         | 
| --METADATA  | Three-column csv with header Sample_R1,Sample_R2,Ref. 
| --OUTDIR | Output folder where .bams and consensus fastas will be piped into.
| -resume  | nextflow will pick up where it left off if the previous command was interrupted for some reason.
| -with-docker ubuntu:18.04 | Runs command with Ubuntu docker.
| -with-trace | Outputs a trace.txt that shows which processes end up in which work/ folders. 
