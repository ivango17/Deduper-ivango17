#!/usr/bin/env python

# The purpose of this script is to take a sam file and remove all PCR replicates

import argparse
import re

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
UMIcount["duplicate"] = 0


# This code block is to save the UMI seqs into a dictionary
# with open(umis) as umis:
with open("./STL96.txt") as umis:
    while True:
        UMI = umis.readline().strip()

        if UMI == "":
            break
        else:
            UMIcount[UMI] = 0


#Here are variables that need to be initialized prior to the loop

curQName = 0
curBitFlag = 0
curCIGAR = 0
curStart = 0
curUMI = 0
curLine = 0

# sam = open(sam)
sam = open("./newTest.sam")
output = open("./newTestOutput.sam", "w")


for line in sam:

    #First conditional statement is to avoid extra sam data

    if line[0] != "@":

        #This series of variable assignments is for setting up the rest of the loop for comparison and file writing
        #The first block saves the previous line info while the second overwrites for the new "current" line

        preQName = curQName
        preBitFlag = curBitFlag
        preUMI = curUMI
        preStart = curStart
        preLine = curLine


        line = line.strip()
        curLine = line.split('\t')

        curQName = curLine[0]
        curQName = curQName.split(":")
        curUMI = curQName.pop()

        curBitFlag = int(curLine[1])
        curCIGAR = curLine[5]
        curStart = StartPosCalc(curCIGAR, int(curLine[3]), BitwiseInterpreter(16, curBitFlag))

        #The below conditional statement is checking wether or not the alignment is valid. If not, it will not be writted to the new file.
        
        if ((curUMI == preUMI) and (curStart == preStart) and (BitwiseInterpreter(16, curBitFlag) == BitwiseInterpreter(16, preBitFlag))) or (curUMI not in UMIcount):
            if curUMI in UMIcount:
                UMIcount["duplicate"] += 1
            else:
                UMIcount["Unknown"] += 1

        else:
            UMIcount[curUMI] += 1
            output.write(f"{line}\n")

    else:
        output.write(line)

    


sam.close()

