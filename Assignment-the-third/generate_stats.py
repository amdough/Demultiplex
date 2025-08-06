# source files / data
matched = "matched.txt"
hopped = "hopped.txt"
unknown_total = 30783962

# getting the match counts
matched_counts = {}
total_matches = 0
with open(matched) as fh:
    for line in fh:
        things = line.strip().split("\t")
        index_pair = things[0]
        count = things[1]
        count = int(count)
        matched_counts[index_pair] = count
        total_matches += count

# getting hopped counts
hopped_counts = {}
total_hopped = 0
with open(hopped) as fh:
    for line in fh: 
        index_pair, count = line.strip().split("\t")
        count = int(count)
        hopped_counts[index_pair] = count
        total_hopped += count

total_reads = total_matches + total_hopped + unknown_total

# starting the markdown report

def sort_by_count(item):
    return item[1] 

with open ("results.md", "w") as out:
    out.write("# Assignment the Third: Demultiplexing Summary\n\n\n")
    out.write(f"**Total reads:** {total_reads}\n\n")
    out.write(f"Total matched reads: {total_matches}\n\n")
    perc_match = (total_matches/total_reads)*100
    out.write(f"Percentage of matched reads: {perc_match:.2f}%\n\n")
    out.write("## Percentage of reads per sample: \n\n")
    out.write("| Index Pair | Count | Percentage (%) |\n")
    out.write("|------------|-------|----------------|\n")
    for index_pair, count in sorted(matched_counts.items(), key=sort_by_count, reverse=True):
        perc = (count/total_matches)*100
        out.write(f"| {index_pair} | {count} | {perc:.2f}% |\n")

    out.write("\n\n## Total Index Swaps\n\n")
    perc_swap = (total_hopped/total_reads)*100
    out.write(f"- Total hopped reads: {total_hopped}\n")
    out.write(f"- Percentage of hopped reads: {perc_swap:.2f}%\n\n")
    
    out.write("## Percentage of hopped reads per sample: \n\n")
    out.write("| Index Pair | Count | Percentage (%) |\n")
    out.write("|------------|-------|----------------|\n")
    for index_pair, count in sorted(hopped_counts.items(), key=sort_by_count, reverse=True):
        perc = (count/total_hopped)*100
        out.write(f"| {index_pair} | {count} | {perc:.2f}% |\n")

    out.write("## Unknown Reads\n")
    perc_unk = (unknown_total/total_reads)*100
    out.write(f"- Unknown read pairs: {unknown_total}\n")
    out.write(f"- Percentage of hopped reads: {perc_unk:.2f}%\n\n")

print("markdown summary written to results.md")