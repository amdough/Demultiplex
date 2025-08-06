#!/usr/bin/env python

import argparse
import bioinfo 
import gzip
import matplotlib.pyplot as plt



def get_args():
    parser = argparse.ArgumentParser(description="per base quality score function")
    parser.add_argument("-f", help="Specify the file name", type=str, required =True)
    parser.add_argument("-l", help="Specify the read label", type=str, required=True)
    parser.add_argument("-o", help="Specify the output filename", type=str, required=True)

    return parser.parse_args()

args=get_args()
f=args.f 
l=args.l 
o=args.o

def avg_qualscores(f):
    dict_sums = {} # total qual score by base
    dict_counts = {} # total number of bases at that pos

    with gzip.open (f, "rt") as fh:    
        while True: 
            header = fh.readline()
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()

            if not (header and seq and plus and qual):
                break

            qual = qual.strip()
            for base_pos, char in enumerate(qual):
                phred_score = bioinfo.convert_phred(char)
                if base_pos not in dict_sums:
                    dict_sums[base_pos] = 0
                    dict_counts[base_pos] = 0
                
                dict_sums[base_pos]+= phred_score
                dict_counts[base_pos] += 1
    mean_scores = [dict_sums[base_pos]/dict_counts[base_pos] for base_pos in sorted(dict_sums)]
    return mean_scores


def plot_quals(mean_scores,l,o):
    if not mean_scores:
        raise ValueError("No quality scores found—check input.")
    x_ax = list(range(len(mean_scores)))
    plt.figure()
    plt.plot(x_ax, mean_scores)
    plt.xlabel("Base Position")
    plt.ylabel("Average Quality Score")
    plt.title(f"Per-base Avg Quality: {l}")
    plt.tight_layout()
    plt.savefig(o)
    plt.close()

if __name__ == "__main__":
    args = get_args()
    scores = avg_qualscores(args.f)
    plot_quals(scores,args.l,args.o)
   
