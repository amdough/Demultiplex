#!/usr/bin/env python

import argparse
import bioinfo 
import gzip



def get_args():
    parser = argparse.ArgumentParser(description="per base quality score function")
    parser.add_argument("-R1", help="Specify read1 file name", type=str, required =True)
    parser.add_argument("-R2", help="Specify read2 file name", type=str, required =True)
    parser.add_argument("-R3", help="Specify read3 file name", type=str, required =True)
    parser.add_argument("-R4", help="Specify read4 file name", type=str, required =True)
    # index file
    parser.add_argument("-I", help="Specify index file name", type=str, required =True)
    return parser.parse_args()

# index = "/projects/bgmp/shared/2017_sequencing/indexes.txt"

def read_for(f):
    header = f.readline().strip()
    seq = f.readline().strip()
    plus = f.readline().strip()
    qual = f.readline().strip()
    return header, seq, plus, qual

args=get_args()
read1=args.R1
read2=args.R2
read3=args.R3
read4=args.R4
index= args.I


index_set = set()
with open(index, "r") as f:
    next(f)
    for line in f:
        last_col = line.strip().split()[-1]  # split by whitespace and take last column
        index_set.add(last_col)
# print(index_set)
    
file_dict= {}
for ind in index_set: 
    R1= open(ind + "_R1.fq", "w")
    R2= open(ind+ "_R2.fq", "w")
    file_dict[ind]= (R1,R2)

with gzip.open (read1, "rt") as r1, \
     gzip.open(read2, "rt") as r2, \
     gzip.open(read3, "rt") as r3, \
     gzip.open(read4, "rt") as r4, \
     open("unk_R1.fq", "w")as u1, \
     open ("unk_R2.fq", "w") as u2, \
     open ("hop_R1.fq", "w") as hop1, \
     open ("hop_R2.fq", "w") as hop2:
    
    hopped ={}
    matched={}
    unk_pairs = 0  

    while True: 
        head1, seq1, plus1, qual1 = read_for(r1)
        head2, seq2, plus2, qual2 = read_for(r2)
        head3, seq3, plus3, qual3 = read_for(r3)
        head4, seq4, plus4, qual4 = read_for(r4)

        if not head1:
            break

        seq3_rc= bioinfo.reverse_complement(seq3)
        index_heads = f"{seq2}-{seq3_rc}"

        # for matched pairs

        if seq2 == seq3_rc and seq2 in index_set:
            newheader1 = f"{head1} {index_heads}"
            newheader2 = f"{head4} {index_heads}"
            file_dict[seq2][0].write(f"{newheader1}\n{seq1}\n{plus1}\n{qual1}\n")
            file_dict[seq2][1].write(f"{newheader2}\n{seq4}\n{plus4}\n{qual4}\n")

            matched[index_heads] = matched.get(index_heads, 0) + 1 

            continue

        # checking for index hoppin'

        elif seq2 in index_set and seq3_rc in index_set and seq2 != seq3_rc:
            newheader1 = f"{head1} {index_heads}"
            newheader2 = f"{head4} {index_heads}"
            hop1.write(f"{newheader1}\n{seq1}\n{plus1}\n{qual1}\n")
            hop2.write(f"{newheader2}\n{seq4}\n{plus4}\n{qual4}\n")

            hopped[index_heads] = hopped.get(index_heads, 0) + 1
            
            continue 

        # unknowns

        else: 
            newheader1 = f"{head1} {index_heads}"
            newheader2 = f"{head4} {index_heads}"
            u1.write(f"{newheader1}\n{seq1}\n{plus1}\n{qual1}\n")
            u2.write(f"{newheader2}\n{seq4}\n{plus4}\n{qual4}\n")

            unk_pairs += 1 


for R1_file, R2_file in file_dict.values():
    R1_file.close()
    R2_file.close()

with open("matched.txt", "w") as matched_out:
    for index, count in matched.items():
        matched_out.write(f"{index}\t{count}\n")

with open("hopped.txt", "w") as hopped_out:
    for index, count in hopped.items():
        hopped_out.write(f"{index}\t{count}\n")

matched_lines = len(matched)
hopped_lines = len(hopped)

print(f"Matched index pairs: {matched_lines} correct matches")
print(f"Hopped index pairs: {hopped_lines} hopped indexes")
print(f"Unknown read pairs: {unk_pairs}")