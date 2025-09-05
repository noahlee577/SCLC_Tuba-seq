"""
hamming.py
Created on July 2nd, 2019 by Myung Chang Lee (Noah Lee)
Calculates hamming distance for sgIDs
Will report ones with smaller distance than the minimum (default 2 or less)

"""

# Import sys to get arguments
import sys
# Import re for regex
import re

import pandas as pd

min_hamm = 2

def hamming(seq1, seq2):

    try:
        seq1 = seq1.strip()
        seq2 = seq2.strip()

        if len(seq1) != len(seq2):
            print("Unclear how to calculate hamming distance of", seq1, "and", seq2)
            return(-1)
    except:
        print("Unclear how to calculate hamming distance of", seq1, "and", seq2)
        return(-1)

    hamming_dist = 0

    for letter1, letter2 in zip(seq1, seq2):
        if letter1 != letter2:
            hamming_dist += 1

    return(hamming_dist)

try:
    in_filename = sys.argv[1]
    if len(sys.argv[1:]) < 1:
        print("Please do note that you can specify filename on commandline")
        print("E.g. python3 hamming.py sgID.xlsx [minimum_hamming_distance]")

    if len(sys.argv[1:]) >= 2:
        min_hamm = int(sys.argv[2])

except Exception as e:
    # If no filename was supplied, then use the default name
    print("Please specify input filename")
    print("E.g. python3 hamming.py sgID.xlsx [minimum_hamming_distance]")

with pd.ExcelFile(in_filename) as sourcefile:
    # Import all the sheets into a dictionary
    sheets = {}

    for curr_sheet in sourcefile.sheet_names:
        sheets[curr_sheet] = sourcefile.parse(curr_sheet)

    # Get first sheet with data
    sgID_data = sheets[sourcefile.sheet_names[0]]

    sgID_dict = dict()

    # Try the column label indexing first,
    # if not labeled properly then try to use numeric index
    try:
        seq_to_name = [{key:value} for key, value in zip(sgID_data.loc[:, 'sgID'], sgID_data.loc[:, 'gRNA'])]

    except:
        seq_to_name = [{key:value} for key, value in zip(sgID_data.iloc[0:,1], sgID_data.iloc[0:,0])]

    for each_dict in seq_to_name:
        for key, value in each_dict.items():
            if key in sgID_dict:
                print("DUPLICATE:", key, value, "with: ", sgID_dict[key])
                sgID_dict[key] = sgID_dict[key] + "/" + value
            else:
                sgID_dict.update(each_dict)

    # Store a list of hamming distance tuples ("gRNA1&gRNA2", HammingDistance)
    distances = [(sgID_dict[x] + "&" + sgID_dict[y], hamming(x, y)) for i,x in enumerate(sgID_dict.keys()) for y in list(sgID_dict.keys())[i+1:]]
    distance_dict = dict()

    for distance in distances:
        distance_dict.update({distance[0]:distance[1]})

    found = False

    for key, value in sorted(distance_dict.items(), key = lambda kv:(kv[1], kv[0])):
        if value <= min_hamm:
            found = True
            print("sgID Pair:", key, "   Hamming Distance:", value)

    if not found:
        print("No sgIDs have hamming distance less than or equal to", min_hamm)
