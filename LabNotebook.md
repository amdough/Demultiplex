## Demultiplexing and Index Swapping 


**52 total files (FASTQ)**

48 FASTQ files that contain acceptable index pairs (read1 and read2 for 24 different index pairs), two FASTQ files with index-hopped reads-pairs, and two FASTQ files undetermined (non-matching or low quality) index-pairs.

**4 FQ files in talapas:**

/projects/bgmp/shared/2017_sequencing/

1294_S1_L008_R1_001.fastq.gz
1294_S1_L008_R2_001.fastq.gz
1294_S1_L008_R3_001.fastq.gz
1294_S1_L008_R4_001.fastq.gz

**note: gzip module in python**

First steps: 

PART 1 (all initial data exploration): 

Data exploration (zcat head, zcat wc -l) → just bash commands for this part

Phred options (+33, +64)

383, 246, 735 reads (*4 = total lines = 1,452,986,940

Zcat <file> | head -2 | tail -1 | wc 

Wc includes new line character!!! So wc total -1 for each line

Barcodes = 8 

Get per base distribution of quality scores for read 1,2, index 1,2

    X axis = base pos
    Y axis = avg QS

Python script! Not Jupyter

Don’t hardcode → use argparse duh doy (file, read length, etc)

Dont use 2D array strategy
Turn in the 4 histograms.

    What is a good quality score cutoff for index reads and biological read pairs to utilize for sample identification and downstream analysis, respectively? Justify your answer.
    How many indexes have undetermined (N) base calls? (Utilize your command line tool knowledge. Submit the command(s) you used. CHALLENGE: use a one-line command)


## PART 2 : Directions

    Pseudo code: write up a plan/strategy for this algorithm (not script)
    Someone who doesn’t know how to code should understand what we’re doing → not python syntax aka
    Write paragraph about the problem
    Describe what output would be most informative
    Write examples (unit tests)
    Include four properly formatted input FASTQ files with read pairs that cover all three categories (dual matched, index-hopped, unknown index)
    Include the appropriate number of properly formatted output FASTQ files given your input files
    If in fastq file, have: 
    Index 1 matched with index 1
    6 or more output files (properly formatted, no comments)
    Develop your algorithm using pseudocode
    Determine high level functions
    Function headers (name and parameters)
    Description/doc string – What does this function do?
    Test examples for individual functions
    Return statement
Example: If you were planning to write the function convert_phred(), you would include something like
def convert_phred(letter: str) -> int:
    '''Takes a single ASCII character (string) encoded in Phred+33 and
    returns the quality score value as an integer.'''
    return qscore
Input: I
Expected output: 40
Markdown file for part 2




Determine which files contain the indexes, and which contain the paired end reads containing the biological data of interest. Create a table and label each file with either read1, read2, index1, or index2.

R1 and R4 - biological data
R2 and R3 - indexes R1 and R2



Determine the length of the reads in each file. *remember wc -1 because of /n*

    $ zcat 1294_S1_L008_R1_001.fastq.gz | head -2 | tail -1 | wc
    102

    $ zcat 1294_S1_L008_R2_001.fastq.gz | head -2 | tail -1 | wc
    9

    $ zcat 1294_S1_L008_R3_001.fastq.gz | head -2 | tail -1 | wc
    9

    $ zcat 1294_S1_L008_R4_001.fastq.gz | head -2 | tail -1 | wc
    102


Determine the phred encoding for these data.

(base) [amdo@n0349 2017_sequencing]$  zcat 1294_S1_L008_R1_001.fastq.gz | sed -n '4~4p' | grep "[abcde]+"| wc -l

*this took too long so got a diff strategy*

```(base) [amdo@login4 2017_sequencing]$ zcat 1294_S1_L008_R1_001.fastq.gz | awk 'NR%4==0 && NR%10000==0' | grep '[a-z]' -m 1 && echo "Phred+64 likely" || echo "Phred+33 likely"```

Phred+33 likely

```(base) [amdo@login4 2017_sequencing]$ zcat 1294_S1_L008_R2_001.fastq.gz | awk 'NR%4==0 && NR%10000==0' | grep '[a-z]' -m 1 && echo "Phred+64 likely" || echo "Phred+33 likely"```


Phred+33 likely

```(base) [amdo@login4 2017_sequencing]$ zcat 1294_S1_L008_R3_001.fastq.gz | awk 'NR%4==0 && NR%10000==0' | grep '[a-z]' -m 1 && echo "Phred+64 likely" || echo "Phred+33 likely"```


Phred+33 likely


```(base) [amdo@login4 2017_sequencing]$ zcat 1294_S1_L008_R4_001.fastq.gz | awk 'NR%4==0 && NR%10000==0' | grep '[a-z]' -m 1 && echo "Phred+64 likely" || echo "Phred+33 likely"```

Phred+33 likely

**psuedocode written in [Answers.md file](Assignment-the-first/Answers.md)

## Part 1 continued

* wrote part1.1.py --> failed (initial runtimes ~ 4hours +)
    [Part1.1 Scripts and outputs](part1.1_outs)

Ex. error ouptut: 

Traceback (most recent call last):
  File "/gpfs/projects/bgmp/amdo/bioinfo/Bi622/Demultiplex/./part1.py", line 37, in <module>
    phred_score = bioinfo.convert_phred(str(char))
  File "/gpfs/projects/bgmp/amdo/bioinfo/Bi622/Demultiplex/bioinfo.py", line 25, in convert_phred
    return ord(letter)-33
           ~~~^^^^^^^^
TypeError: ord() expected a character, but string of length 2 found
Command exited with non-zero status 1
	Command being timed: "python ./part1.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -l Read4 -o Read4"
	User time (seconds): 0.96
	System time (seconds): 0.04
	Percent of CPU this job got: 242%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.41
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 60008
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 13379
	Voluntary context switches: 629
	Involuntary context switches: 6
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 1

- this reveals issues with plot generator, as well as a bug in the python script that has an issue with input type in the phred input

- found another bug in the 1.1.py script where there was essentially an infinite loop happening because of an if/else statement. removed the else so that the dictionary was always updating with the += phred score instead of resetting to zero. le oops... first version .py visible here: [first python script](part1.1_outs/part1.py)

-- made second version of .py and .sh scripts and submitted on talapas 
-- successful .png plots generated, no errors reported: 

	Command being timed: "python ./part1.2.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -l Read1 -o Read1"
	User time (seconds): 7124.24
	System time (seconds): 1.83
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:59:03
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 66880
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 48518
	Voluntary context switches: 754
	Involuntary context switches: 2849
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
	Command being timed: "python ./part1.2.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -l Index1 -o Index1"
	User time (seconds): 832.23
	System time (seconds): 0.37
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 13:59.52
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 65092
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 105890
	Voluntary context switches: 1969
	Involuntary context switches: 998
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
	Command being timed: "python ./part1.2.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -l Index2 -o Index2"
	User time (seconds): 835.62
	System time (seconds): 0.31
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 13:57.53
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 65412
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 48656
	Voluntary context switches: 655
	Involuntary context switches: 246
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
	Command being timed: "python ./part1.2.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -l Read2 -o Read2"
	User time (seconds): 7251.04
	System time (seconds): 1.87
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 2:01:07
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 68356
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 153610
	Voluntary context switches: 660
	Involuntary context switches: 3779
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0


run time was around 2 hours, though from the output R2 and R3 (indexes) only took about 15 minutes. could have probably split up the index runs and read runs to shorten it? but idk

Correct scripts: 
    [python](part1.2.py),
    [bash](part1.2.sh)
