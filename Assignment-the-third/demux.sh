#!/bin/bash

#SBATCH --account=bgmp                          #REQUIRED: which account to use
#SBATCH --partition=bgmp                        #REQUIRED: which partition to use
#SBATCH --cpus-per-task=8                 #optional: number of cpus, default is 1
#SBATCH --mem=16GB
#SBATCH --mail-user=amdo@uoregon.edu     #optional: if you'd like email
#SBATCH --mail-type=ALL                   #optional: must set email first, what type of email you want
#SBATCH --job-name=demux               #optional: job name
#SBATCH --output=part3_%j.out       #optional: file to store stdout from job, %j adds the assigned jobID
#SBATCH --error=part3_%j.err

f1=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz

f2=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz

f3=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz

f4=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz

I=/projects/bgmp/shared/2017_sequencing/indexes.txt


/usr/bin/time -v python ./part3.1.py -R1 $f1 -R2 $f2 -R3 $f3 -R4 $f4 -I $I



