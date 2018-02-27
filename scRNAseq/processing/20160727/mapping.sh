## script for calling mapping


#!/bin/bash

ispet=0
fastq=($(ls fastq/*.fq.gz))
genome=/lustre/jmlab/resources/genomes/subread/mm10_ERCC

source ${HOME}/Code/mapping/multi_align.sh