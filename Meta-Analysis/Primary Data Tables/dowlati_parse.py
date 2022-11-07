"""
dowlati_parse.py
Created on June 10th, 2019 by Myung Chang Lee (Noah Lee)
Parses Supp Table 2 from Dowlati et al 2016

"""

# Import sys to get arguments
import sys
# Import re for regex
import re

import pandas as pd
import math

# Get the first argument from command line, filename
try:
    in_filename = sys.argv[1]
except Exception as e:
    # If no filename was supplied, then use the default name
    in_filename = "Dowlati et al 2016 Supp Table 2.xlsx"

    out_filename = "Dowlati et al 2016 Supp Table 2_Parsed.xlsx"

#if len(sys.argv[1:]) != 1:
#    print("Please do note that you can specify input and filename on commandline")
#    print("E.g. python3 dowlati_parse.py digest_result.txt")

mutation_key = {
    "M":"Missense",
    "N":"Nonsense",
    "J":"Frameshift",
    "B":"Splice-site",
    "I":"Insertion",
    "D":"Deletion",
    "U":"Unknown",
    "F":"Fusion",
    "S":"Splice",
    "A":"Amplification",
    "L":"Loss",
    "R":"Rearrangement"
}

data = pd.read_excel(in_filename, index_col = 0)

# Keep only the mutation data
data = data.iloc[7:,]

print(data)

gene_totals = dict()

# Get all the alterations for each gene
for row in data.itertuples():
    gene = row[0]
    gene_alterations = dict()
    for element in row[1:]:
        try:
            # The following will result in ValueError if element is not NaN
            float(element)
        except ValueError as e:
            # only count unique alterations (e.g. MMJ -> MJ)
            # Also sort alphabetically to avoid MJ <-> JM being unique
            uniques = sorted(set([letter for letter in element]))
            element = "".join(uniques)

            if element in gene_alterations:
                gene_alterations[element] += 1
            else:
                gene_alterations[element] = 1

    gene_totals.update({gene:gene_alterations})

# store parsed data to convert to df later
parsed_data = dict()

# Parse into usable format
for gene in gene_totals:
    total_count = 0
    alteration_str = ""
    first_alteration = True
    for alteration, count in gene_totals[gene].items():
        if not first_alteration:
            alteration_str += ", "

        first_alteration = False
        alteration_str += (str(count) + " ")
        first_element = True

        for element in alteration:
            if element not in mutation_key:
                print("In gene", gene + ", Unknown alteration:", element)
            else:
                if not first_element:
                    alteration_str += "+"

                alteration_str += (mutation_key[element])
                first_element = False
        total_count += count

    alteration_str = str(total_count) + "/50 patients with alterations, with " + alteration_str + " (Supp Table 2)"
    parsed_data.update({gene:alteration_str})

parsed_df = pd.DataFrame(parsed_data, index = ["Alterations"]).T

print(parsed_df)

parsed_df.to_excel(out_filename)
