#!/usr/bin/env python


# The purpose of this script is to take a sam file and remove all PCR replicates from a sorted sam file.
# Author: Ian VanGordon
# 10/15/2023


import argparse
import re



#===========================================================================================================================================================================================================================
# Argparse setup
#===========================================================================================================================================================================================================================

def get_args():
    parser = argparse.ArgumentParser(description="This program generates an output file from a sorted sam file without PCR replicates. This script does not take hard clipping into account.")
    parser.add_argument("-f", "--file", help="What is the filepath for the sam file to be read?", type=str)
    parser.add_argument("-o", "--outfile", help="What should this file to be called?", type=str)
    parser.add_argument("-u", "--umi", help="What is the filepath for the file that contains all of the UMIs?", type=str)
    parser.add_argument("-s", "--summaryfile", help="What should the summary file be called?", type=str)
    parser.add_argument("-d", "--duplicatesfile", help="What should the duplicates file be called?", type=str)
    return parser.parse_args()

args = get_args()

sam = args.file
output = args.outfile
umis = args.umi

# Conditional statements making summary and duplicates files optional
if args.summaryfile != None:
    summary = args.summaryfile
else:
    print("Summary file option not selected, no summary will be output.")

if args.duplicatesfile != None:
    duplicates = args.duplicatesfile
else: 
    print("Duplicates file option not selected, no duplicates file will be output.")



#===========================================================================================================================================================================================================================
# Function definitions
#===========================================================================================================================================================================================================================

def BitwiseInterpreter(val, bitflag):
    '''This function takes a bitwise flag and an int of interest. Returns a true or false using if-else statement comparing integer and bitwiseflag.'''
    if ((bitflag & val) == val):
        return True
    else:
        return False

def StartPosCalc(CIGARstring, startPos, minusStrand = False):
    '''This function takes a CIGAR string and starting position from a sam file header and determines the start position for that alignment on the genome to compare possible duplicates to each other. The function will parse the CIGAR string and return the actual start int.'''
    cString = re.findall(r'(\d+)(\w)', CIGARstring)
    
    if minusStrand == False:

        if cString[0][1] == 'S':
            return startPos - int(cString[0][0])
        else:
            return startPos

    else:
        bpCounter = 0
        for i in range(len(cString)):
            if ((i != 0) or (cString[0][1] != 'S')) and (cString[i][1] != "I"):
                bpCounter += int(cString[i][0])
                
        return startPos + bpCounter


    
#===========================================================================================================================================================================================================================
# Reading unique molecular identifier (UMI) file
#===========================================================================================================================================================================================================================

#This dictionary will serve to hold all of the UMI seqs as well as their counts, an unknow counter will also be added
UMIcount = {}
UMIcount["Unknown"] = 0
UMIcount["Duplicate"] = 0


# This code block is to save the UMI seqs into a dictionary
# with open(umis) as umis:
with open(umis) as umis:
    while True:
        UMI = umis.readline().strip()

        if UMI == "":
            break
        else:
            UMIcount[UMI] = 0



#===========================================================================================================================================================================================================================
# PCR Duplicates removal
#===========================================================================================================================================================================================================================

# Here are variables that need to be initialized prior to the loop

curChrom = 0
duplicateHunter = {}
preLineCount = 0
postLineCount = 0

# Opening files

output = open(output, "w")
sortedSam = open(sam)

if args.duplicatesfile != None:
    duplicatesFile = open(duplicates, "w")


for line in sortedSam:

    #First conditional statement is to ensure header inclusion

    if line[0] != "@":

        preLineCount += 1

        preChrom = curChrom

        # This series of variable assignments is for setting up the rest of the loop for comparison and file writing
        line = line.strip()
        curLine = line.split('\t')

        # These variables hold information about the current alignment in the sam file to compare against duplicateHunter dict
        curQName = curLine[0]
        curQName = curQName.split(":")
        curUMI = curQName.pop()
        curBitFlag = BitwiseInterpreter(16, int(curLine[1]))
        curCIGAR = curLine[5]
        curStart = StartPosCalc(curCIGAR, int(curLine[3]), curBitFlag)
        curChrom = curLine[2]


        # The statement below wipes the comparison dictionary if the current chromosome changes

        if curChrom != preChrom:

            duplicateHunter.clear()

        
        # The below tuple will be used as a key for duplicateHunter

        currentSet = (curUMI, curBitFlag, curStart, curChrom)


        # This conditional statement compares currentset to the duplicateHunter dictionary to determine if the current line is a PCR duplicate

        if currentSet not in duplicateHunter:

            duplicateHunter[currentSet] = 1
            UMIcount[curUMI] += 1
            postLineCount += 1

            output.write(f"{line}\n")

        else:

            if args.duplicatesfile != None:
                duplicatesFile.write(f"{line}\n")

            if curUMI in UMIcount:
                UMIcount["Duplicate"] += 1
            else:
                UMIcount["Unknown"] += 1


    elif args.duplicatesfile != None:
        output.write(f"{line}")
        duplicatesFile.write(f"{line}")

    else:
        output.write(line)



#===========================================================================================================================================================================================================================
# Summary file generation
#===========================================================================================================================================================================================================================

if args.summaryfile != None:

    summary = open(summary, "w")

    summary.write(f"Summary information about deduplicated {sam}\n\n")
    summary.write(f"Alignments before processing: {preLineCount}\n")
    summary.write(f"Alignments after processing: {postLineCount}\n")
    summary.write(f"Precent surviving: {100 * (postLineCount/preLineCount)}%\n\n")

    propDup = (UMIcount["Duplicate"]*100)/preLineCount
    numDup = UMIcount["Duplicate"]

    summary.write(f"Number of alignments that were PCR duplicates: {numDup}\n")
    summary.write(f"Duplicate\t{propDup}% of original alignments were duplicates\n\n")
    summary.write(f"UMI\tCount\n")

    del UMIcount["Duplicate"]

    for key in UMIcount:
        summary.write(f"{key}\t{UMIcount[key]}\t{(100*UMIcount[key])/postLineCount}%\n")

    summary.close()



#===========================================================================================================================================================================================================================
# Closing files
#===========================================================================================================================================================================================================================

sortedSam.close()
output.close()

if args.duplicatesfile != None:
    duplicatesFile.close()

print("Deduplication Complete!")