#!/bin/bash

#SBATCH --account=bgmp                          #REQUIRED: which account to use
#SBATCH --partition=bgmp                        #REQUIRED: which partition to use
#SBATCH --cpus-per-task=8                 #optional: number of cpus, default is 1
#SBATCH --mem=16GB
#SBATCH --mail-user=amdo@uoregon.edu     #optional: if you'd like email
#SBATCH --mail-type=ALL                   #optional: must set email first, what type of email you want
#SBATCH --job-name=demux_p1               #optional: job name
#SBATCH --output=part1_%j.out       #optional: file to store stdout from job, %j adds the assigned jobID
#SBATCH --error=part1_%j.err

f1=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz

f2=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz

f3=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz

f4=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz



/usr/bin/time -v python ./part1.2.py -f $f1 -l Read1 -o Read1



/usr/bin/time -v python ./part1.2.py -f $f2 -l Index1 -o Index1



/usr/bin/time -v python ./part1.2.py -f $f3 -l Index2 -o Index2



/usr/bin/time -v python ./part1.2.py -f $f4 -l Read2 -o Read2


