## script for calling mapping


#!/bin/bash

PATH=/platform/lsf/gui/2.0/bin:/platform/lsf/perf/1.2/bin:/platform/lsf/7.0/linux2.6-glibc2.3-x86_64/etc:/platform/lsf/7.0/linux2.6-glibc2.3-x86_64/bin:/usr/kerberos/bin:/usr/local/bin:/bin:/usr/bin:/usr/X11R6/bin:/usr/local/matlab/toolbox/local:/usr/local/matlab2008bclient/bin:/home/richar02/bin:/home/richar02/bin/FastQC:/lustre/jmlab/software/bin:/home/richar02/bin/subread-1.5.1-Linux-x86_64/bin:/home/richar02/bin/samtools-1.3.1


ispet=0
fastq=($(ls fastq/*.fq.gz))
genome=/lustre/jmlab/resources/genomes/subread/mm10_ERCC

source ${HOME}/Code/mapping/multi_align.sh


