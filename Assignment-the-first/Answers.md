# Assignment the First

## Part 1
1. Be sure to upload your Python script. Provide a link to it here:

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz |  |  |  |
| 1294_S1_L008_R2_001.fastq.gz |  |  |  |
| 1294_S1_L008_R3_001.fastq.gz |  |  |  |
| 1294_S1_L008_R4_001.fastq.gz |  |  |  |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
    2. **YOUR ANSWER HERE**
    3. **YOUR ANSWER HERE**
    
## Part 2
1. Define the problem

 We have a large amount of prepped genomic library data from multiple different samples, and need to demultiplex the data in order to perform downstream analyses like quality filtering index reads and determining rates of index swapping or undetermined index pairs. We have 24 indexed (dual matched) libraries that were submitted for sequencing, and four large input files to parse through: 

    1294_S1_L008_R1_001.fastq.gz -- representing biological read 1
    1294_S1_L008_R2_001.fastq.gz -- representing index read 1
    1294_S1_L008_R3_001.fastq.gz -- representing index read 2 (reverse-complemented)
    1294_S1_L008_R4_001.fastq.gz -- representing biological read 2

To demultiplex, we need to parse through each file simultaneously and match each index1/index2 from each read pair (R1, R4) to the list of 24 indexes. They will fall in one of three categories:   
    * correctly matched pairs --> index 1 matches index 2 and they are one of the 24 known index libraries   
    * index hopped pairs --> both are valid indexes but they do not match (and where did they mismatch from?)   
    * unknown pairs --> one or both are not on the list, or fall below quality threshold, have Ns, etc. 

We also need to:   
    * append the sequence of the index-pair to the header of BOTH reads in all of the output FASTQ files for ALL categories    
    * report:  
            the number of read-pairs with properly matched indexes (per index-pair),  
            the number of read pairs with index-hopping observed. 
            the number of read-pairs with unknown index(es).



2. Describe output

We are expecting to output:  
        48 FASTQ files containing acceptable index pairs (ex: Index A: IA_Read 1 and IA_Read 2 → for all 24 indexes)  
        2 FASTQ files with index-hopped reads and their bio data pairs  
        2 FASTQ files with undetermined (non-matching or low quality) index-pairs  
        == 52 fq output files total


3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
[R1 test file](../TEST-input_FASTQ/R1_tinytest.fq)
[R2 test file](../TEST-input_FASTQ/R2_tinytest.fq)
[R3 test file](../TEST-input_FASTQ/R3_tinytest.fq)
[R4 test file](../TEST-input_FASTQ/R4_tinytest.fq)

[output text test files](../TEST-output_FASTQ)

Mock report: 

    2 matched index pairs   
    1 index hopped pair   
    1 unknown index pair   
    x2 per read --> 8 output files   


4. Pseudocode


STORE indexes as a list.  
READ in and open all four .fastq files (R1, R2, R3, R4) in parallel  
SET header, seq/index, plus, quality


    LOOP through all four files per one record (four lines) at a time
        SET count for matched, hopped, and unknown index pairs = 0
        GET index from R2 
        GET index from R3 and reverse complement it
        CHECK both for Ns, low quality 
        IF one or both have low quality or contain N:
            THEN send to unknown
            INCREMENT unknown
        ELSE IF index 1 == index 2 and is in known index list:
            THEN send R1 to matched_R1.fq file, send R4 to matched_R2.fq file
            INCREMENT matched
        ELSE IF both are matched to known index but do not match eachother:
            THEN send R1 to hopped_R1.fq file, R4 to hopped_R2.fq file
            INCREMENT hopped 





5. High level functions. For each function, be sure to include:



        def convert_phred(letter: str) -> int:
            '''Converts a single character into a phred score'''
            return ord(letter)-33 

        test example:     
        assert convert_phred("I") == 40, "wrong phred score for 'I'"


        def reverse_complement(seq:str) -> str:
            '''Returns the reverse compliment of a DNA sequence'''
            SET dictionary to matched complement pairs: A:T, G:C, C:G, T:A
            return new_seq

        test example: 
        "GTCA" --> "TGAC"
        assert reverse_complement("GTCA")=="TGAC", "Incorrect reverse complement"


        def quality_filter(index_seq:str,qual_str:str, threshold:int = 30)-> bool:
            '''Returns True if all bases exceed threshold value, no N present'''
            IF 'Nn' in index_seq:
                return False
            ELSE
                return all(convert_phred) >= threshold 
        
        test example:
        "TATGCCC", "IIIIIII" --> True

        def qual_score(phred_score: str) -> float:
            '''Fuction that iterates through the phred_score string, calculates quality score using convert_phred function, and returns the calculated average quality score of the entire string'''
            average=0
            for i, letter in enumerate(phred_score):
                average += convert_phred(phred_score[i])
            average /= len(phred_score)
            return average

    
