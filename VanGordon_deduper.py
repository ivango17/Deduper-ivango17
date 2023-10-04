#!/usr/bin/env python

# The purpose of this script is to take a sam file and remove all PCR replicates

import argparse

# def get_args():
#     parser = argparse.ArgumentParser(description="This program generates an output file from a sam file without PCR replicates.")
#     parser.add_argument("-f", "--filename", help="What is the filepath for the sam file to be read", type=str)
#     parser.add_argument("-o", "--outputfile", help="What do you want this file to be called", type=str)
#     parser.add_argument("-u", "--umifile", help="What is the filepath for the file that contains all of the UMIs", type=str)
#     return parser.parse_args()

# args = get_args()

# sam = args.filename
# output = args.outputfile
# umis = args.umifile

def BitwiseInterpreter(val, bitflag):
    '''This function takes a bitwise flag and an int of interest. Returns a true or false using if-else statement comparing integer and bitwiseflag.'''
    if ((bitflag & val) == val):
        return True
    else:
        return False

def StartPosCalc(CIGARstring, startPos):
    '''This function takes a CIGAR string and starting position from a sam file header and determines the start position for that alignment on the genome to compare possible duplicates to each other. The function will parse the CIGAR string and return the actual start int.'''
    if 'S' in CIGARstring:
        for i in range(len(CIGARstring)):
            if CIGARstring[i] == 'S':
                clipped = int(CIGARstring[0:i])
                break
        return startPos - clipped

    else:
        return startPos


#This dictionary will serve to hold all of the UMI seqs as well as their counts, an unknow counter will also be added
UMIcount = {}
UMIcount["Unknown"] = 0


# This code block is to save the UMI seqs into a dictionary
# with open(umis) as umis:
with open("./STL96.txt") as umis:
    while True:
        UMI = umis.readline().strip()

        if UMI == "":
            break
        else:
            UMIcount[UMI] = 0


# sam = open(sam)
sam = open("./test.sam")


sam.close()

