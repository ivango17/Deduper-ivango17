#!/usr/bin/env python

# The purpose of this script is to take a sam file and remove all PCR replicates

import argparse
import re
import pysam

def get_args():
    parser = argparse.ArgumentParser(description="This program generates an output file from a sam file without PCR replicates. This script does not take hard clipping into account.")
    parser.add_argument("-f", "--filename", help="What is the filepath for the sam file to be read", type=str)
    parser.add_argument("-o", "--outputfile", help="What do you want this file to be called", type=str)
    parser.add_argument("-u", "--umifile", help="What is the filepath for the file that contains all of the UMIs", type=str)
    parser.add_argument("-s", "--summaryfile", help="What should the summary file be called?", type=str)
    return parser.parse_args()

args = get_args()

sam = args.filename
output = args.outputfile
umis = args.umifile
summary = args.summaryfile

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
            if (i != 0) or (cString[0][1] != 'S'):
                bpCounter += int(cString[i][0])
                
        return startPos + bpCounter


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


# This loop adjusts the start position for every read in the sam file

sam = open(sam)
samAdj = open("./Adjusted.sam", "w")

for line in sam:

    if line[0] != "@":

        line = line.strip('\n')
        line = line.split()
        adjPos = StartPosCalc(line[5], int(line[3]), BitwiseInterpreter(16, int(line[1])))
        line[3] = adjPos
        samAdj.write(f"{line[0]}\t{line[1]}\t{line[2]}\t{line[3]}\t{line[4]}\t{line[5]}\t{line[6]}\t{line[7]}\t{line[8]}\t{line[9]}\t{line[10]}\n")

    else:
        samAdj.write(f"{line}")

sam.close()
samAdj.close()


# This code block is to sort the sam file with adjusted start positions

sortedSam = pysam.sort("-o", "sorted.sam", "Adjusted.sam")


#Here are variables that need to be initialized prior to the loop

curQName = 0
curBitFlag = 0
curCIGAR = 0
curStart = 0
curUMI = 0
curLine = 0
curChrom = 0

duplicateHunter = {}


output = open(output, "w")
summary = open(summary, "w")
duplicates = open("./Duplicates.sam", "w")
sortedSam = open("./sorted.sam")

preLineCount = -1
postLineCount = -1


for line in sortedSam:

    #First conditional statement is to avoid extra sam data

    if line[0] != "@":

        preLineCount += 1

        #This series of variable assignments is for setting up the rest of the loop for comparison and file writing
        #The first block saves the previous line info while the second overwrites for the new "current" line

        preQName = curQName
        preBitFlag = curBitFlag
        preUMI = curUMI
        preStart = curStart
        preLine = curLine
        preChrom = curChrom


        line = line.strip()
        curLine = line.split('\t')

        curQName = curLine[0]
        curQName = curQName.split(":")
        curUMI = curQName.pop()

        curBitFlag = BitwiseInterpreter(16, int(curLine[1]))
        curCIGAR = curLine[5]
        curStart = curLine[3]
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

            output.write(f"{line}")

        else:
            duplicates.write(f"{line}")

            if curUMI in UMIcount:
                UMIcount["Duplicate"] += 1
            else:
                UMIcount["Unknown"] += 1


    else:
        output.write(line)
        duplicates.write(line)


# Everything below is simply generating the summary file.

summary.write(f"Summary information about deduplicated {sam}\n\n")
summary.write(f"Alignments before processing: {preLineCount}\n")
summary.write(f"Alignments after processing: {postLineCount}\n")
summary.write(f"Precent surviving: {100 * (postLineCount/preLineCount)}%\n\n")

propDup = (UMIcount["Duplicate"]*100)/preLineCount

summary.write(f"Duplicate\t{propDup}% of original alignements were duplicates\n\n")
summary.write(f"UMI\tCount\n")

del UMIcount["Duplicate"]

for key in UMIcount:
    summary.write(f"{key}\t{UMIcount[key]}\t{(100*UMIcount[key])/postLineCount}%\n")


sortedSam.close()
output.close()
summary.close()
duplicates.close()

