# Usage: python converter.py #_cols "pasted_text"
# Created on April 21, 2021
# @author Myung Chang Lee (Noah Lee)

import sys

text_array = sys.argv[2].strip().split()

num_cols = int(sys.argv[1])

macrolist_array = ['\t'.join(text_array[num_cols * i: num_cols * i + num_cols])
                   for i in range(0, int(len(text_array) / num_cols))]

with open("output.tsv", "wt") as outfile:
    outfile.write('\n'.join(macrolist_array))
