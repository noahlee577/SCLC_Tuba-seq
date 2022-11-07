"""
hamming.py
Created on July 2nd, 2019 by Myung Chang Lee (Noah Lee)
Modified on April 2nd, 2021 by Myung Chang Lee (Noah Lee) for barcode comparisons
Calculates hamming distance for two barcodes
Will output a matrix of hamming distances between barcodes (a table)

"""

import sys
import re
import pandas as pd


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
        print("E.g. python3 hamming.py sgID.xlsx")

except Exception as e:
    # If no filename was supplied, then use the default name
    print("Please specify input filename")
    print("E.g. python3 hamming.py sgID.xlsx")

with pd.ExcelFile(in_filename, engine="openpyxl") as sourcefile:
    # Import all the sheets into a dictionary
    sheets = {}

    for curr_sheet in sourcefile.sheet_names:
        # read all columns as str
        columns = sourcefile.parse(curr_sheet).columns
        converters = {column: str for column in columns}
        sheets[curr_sheet] = sourcefile.parse(
            curr_sheet, converters=converters)

    # Get first sheet with data
    sgID_data = sheets[sourcefile.sheet_names[0]]

    print("Parsing sheet", sourcefile.sheet_names[0])

    sgID_dict = dict()

    # Try the column label indexing first,
    # if not labeled properly then try to use numeric index
    try:
        for sample, BC, freq in zip(sgID_data.loc[:, 'Sample'], sgID_data.loc[:, 'BC'], sgID_data.loc[:, '#']):
            key = " ".join([sample, BC, freq])
            value = BC
            sgID_dict.update({key: value})

    except:
        for sample, BC, freq in zip(sgID_data.iloc[0:, 0], sgID_data.iloc[0:, 2], sgID_data.iloc[0:, 3]):
            key = " ".join([sample, BC, freq])
            value = BC
            sgID_dict.update({key: value})

    results_list = {}

    for key, BC1 in sgID_dict.items():
        hamming_list = []
        for BC2 in sgID_dict.values():
            hamming_list.append(hamming(BC1, BC2))

        results_list.update({key: hamming_list})

    df = pd.DataFrame(results_list, index=[key for key in sgID_dict.keys()])

    print(df)

    df.to_excel("BC Hamming Matrix.xlsx")
