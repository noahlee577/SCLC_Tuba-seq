"""
SCLC_Mutation_Counter.py
Created on May 1st, 2019 by Myung Chang Lee (Noah Lee)
Counts the total # of patients and total # patients with alterations

"""

# Import pandas for opening xlsx
import pandas as pd
# Import sys to get arguments
import sys
# Import re for regex
import re


def get_counts(row):
    alteration_count = 0
    patient_count = 0

    # Find fractional statements in each element (40/67 patients)
    for element in row:
        match_obj = re.search(r"(\d+)/(\d+)", str(element))

        # Add all fractions in the column
        while match_obj:
            # Only keep what's left after the first match in the string
            element = element[match_obj.span()[1]:]

            alteration_count += int(match_obj.group(1))
            patient_count += int(match_obj.group(2))

            match_obj = re.search(r"(\d+)/(\d+)", str(element))

    proportion = alteration_count / patient_count if patient_count > 0 else 0

    return pd.Series({'# Patients with Alterations': alteration_count, '# Total Patients Profiled': patient_count, 'Proportion of Patients with Alterations': proportion})


# Get the first argument from command line, filename
try:
    in_filename = sys.argv[1]
except Exception as e:
    # If no filename was supplied, then use the default name
    in_filename = "SCLC Mutation Summary All.xlsx"
try:
    out_filename = sys.argv[2]
except Exception as e:
    # If no output filename was specified,
    # append _cleaned to the input filename and use that as output filename
    out_filename = in_filename.split(".xlsx")[0] + "_counted.xlsx"

if len(sys.argv[1:]) != 2:
    print("Please do note that you can specify input and output filenames on commandline")
    print("E.g. python3 SCLC_Mutation_Counter.py DB.xlsx output.xlsx")
    print("Will currently take", in_filename, "and export onto", out_filename)

with pd.ExcelFile(in_filename, engine="openpyxl") as sourcefile:
    # Import all the sheets into a dictionary
    sheets = {}

    for curr_sheet in sourcefile.sheet_names:
        sheets[curr_sheet] = sourcefile.parse(curr_sheet)

    # Get second sheet with the mutation data
    mutation_data = sheets[sourcefile.sheet_names[0]]

    # get the columns with mutation counts
    study_data = mutation_data.iloc[0:, 0:]

    mutation_data = pd.concat(
        [mutation_data, study_data.apply(get_counts, axis=1)], axis=1)

    print(mutation_data)

    mutation_data.to_excel(out_filename)
