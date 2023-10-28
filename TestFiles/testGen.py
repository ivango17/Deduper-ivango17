#!/usr/bin/env python

file = "./test.sam"

count = 0


with open(file, "r") as rf, open("./newTest.sam", "w") as out:

    for line in rf:

        if line[0] != "@":
            line = line.strip()

            lineList = line.split()

            out.write(f"{line}\n")

            lineList[1] = 16

            for i in lineList:
                out.write(f"{i}\t")
            out.write("\n")

            for i in lineList:
                out.write(f"{i}\t")
            out.write("\n")


with open(file, "r") as rf, open("./newTest.sam", "a") as out:

    for line in rf:
        
        count += 1

        if (line[0] != "@") and (count % 3 == 0):
            
            line = line.strip()

            lineList = line.split()

            lineList[0] = "NS500451:154:HWKTMBGXX:1:11101:25533:AAAAAAAA"

            for i in lineList:
                out.write(f"{i}\t")

            out.write("\n")

            for i in lineList:
                out.write(f"{i}\t")

            out.write("\n")
                    

