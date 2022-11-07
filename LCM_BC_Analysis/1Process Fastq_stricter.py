# -*- coding: utf-8 -*-
"""
# 1Process Fastq to count reads mapped to each sgID-BC region
@author: Chuan Li (chuanli@stanford.edu)
@editor: Myung Chang Lee (Noah Lee, mclee18@stanford.edu)

Each read is expected to contain an 8-nucleotide sgID region followed by a 30-nucleotide barcode (BC) region (GCNNNNNTANNNNNGCNNNNNTANNNNNGC),
and each of the 20 Ns represent random nucleotides with roughly equal representation of A, T, G and C.
The sgID region identifies the putative tumour suppressor gene being targeted,
for which we require perfect match between the sequence in the forward read and one of the 102 forward sgIDs with known sequences.
We require the forward and reverse read to agree completely within the 30 nucleotide sequence to be further processed.

"""

###
# 1. Load packages
import regex as re
import gzip
import sys


def main(File1='S1_R1.fastq.gz'):
    '''
    This function reads the merged reads, merge reads together and get high quality ones

    File1: Forward read file
    return the a dictionary of barcodes
    '''
    # 1. Read in files, get the second line of sequence,
    # the sequencing quality is so high that there is no need to consider if the forward and reverse matches
    File2 = File1.replace('_R1', '_R2')

    # Open the two fastq files
    f1 = gzip.open(File1, 'rt')
    f2 = gzip.open(File2, 'rt')

    sgIDBCdict = {}

    line_seq1 = 1

    print("Processing", File1)

    while (line_seq1):
        # read in the files line by line
        # remove the first and third lines
        # store the second line into line_seq, this is the sequence
        # store the fourth line into line, this is the quality score

        # for R1 file
        line1 = f1.readline().rstrip()  # skip the first line
        line_seq1 = f1.readline().rstrip()  # sequence line
        line1 = f1.readline().rstrip()  # skip the third line
        line_qua1 = f1.readline().rstrip()  # get the sequencing quality

        # for R2 file
        line2 = f2.readline().rstrip()
        line_seq2 = f2.readline().rstrip()
        line2 = f2.readline().rstrip()
        line_qua2 = f2.readline().rstrip()

        line_seq2_complement = revcom(line_seq2)

        # find the locator for R1
        # regex = re.compile('GA' + (sgID) + (random_barcode) + 'ATGCCCAAGAAG') # the last piece nt for alignment purpose
        # 20 nt for alignment purpose
        # Original fuzzy criteria: R1 = e<2, R2 = e<3
        regexR1 = re.compile(
            'GA' + '(........)' + '(GC.....TA.....GC.....TA.....GC)' + '(ATGCCCA){e<2}')
        # make criteria looser for R2
        regexR2 = re.compile(
            '(GC.....TA.....GC.....TA.....GC)' + '(ATGCCCA){e<3}')

        if regexR1.search(line_seq1) and regexR2.search(line_seq2_complement):
            # if we can see patterns in both reads
            k1 = regexR1.search(line_seq1)  # align R1
            R1sgID = k1.group(1)  # upstream of R1
            R1BC = k1.group(2)  # downstream of R1

            k2 = regexR2.search(line_seq2_complement)  # align R2
            R2BC = k2.group(1)  # downstream of R2_RC

            # combine reads only if barcode are correct from both ends
            # use the sequence of R1sgID as the sgID
            # sgID = getsgID(sgIDDict, R1sgID=R1sgID)
            # Just use sgID sequence
            sgID = R1sgID

            # if exists barcode and sgID combination
            if R1BC == R2BC and ('N' not in R1BC):
                myKey = sgID + "," + R1BC
                if myKey in sgIDBCdict:
                    sgIDBCdict[myKey] += 1
                else:
                    sgIDBCdict[myKey] = 1
    f1.close()
    return(sgIDBCdict)

# reverse complement


def revcom(text):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement[base] for base in text[::-1])

# the following three function getsgID, compareText and hamming_distance are nested function


def getsgID(sgIDDict, R1sgID="GTAAGGAG"):
    for gene, SGIDSeq in sgIDDict.items():
        if R1sgID == SGIDSeq:
            return(gene)
    return("None")


# MAIN PROGRAM

# read in sgID
sgIDDict = {}
with open('sgID Master File.tsv', 'r') as f:
    for line in f:
        request = line.split("\t")
        sgIDDict[request[1]] = request[4].rstrip()

Files = ["Fastq1"]

try:
    Files = sys.argv[1:]

    if len(Files) < 1:
        print("Please input FastQ files names to process")
        exit()

except Exception:
    print("Please input FastQ file names to process")
    exit()


for File in Files:
    sourceDir = "./All Fastq Files/"
    sgIDBCdict = main(sourceDir + File)
    outputdir = "./Output/"
    FileOut = outputdir + File + '_MergeReadOut.csv'

    with open(FileOut, "wt") as out:
        for w in sorted(sgIDBCdict, key=sgIDBCdict.get, reverse=True):
            out.write(",".join(map(str, [w, sgIDBCdict[w]])) + "\n")
