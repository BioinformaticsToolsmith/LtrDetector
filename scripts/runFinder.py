import subprocess
import os
import shutil
import time
import sys

def runLTRFinder(chromFile, resultDir, name):

    params = ['../../tools/LTR_Finder/source/ltr_finder', chromFile]

    output = subprocess.check_output(params).split("\n")

    convertFinder(output,resultDir, name)



def convertFinder(lines, bedFile, name):

    repeats = []

    if not os.path.isdir(resultDir):
        os.mkdir(resultDir)

    info = []
    for l in lines:
        fields = l.split()
        print fields

        if len(fields) > 0:
            if fields[0] == "5'-LTR":
                info.extend([fields[2], fields[4]])
            elif fields[0] == "3'-LTR":
                info.extend([info[0], fields[2], fields[4], info[1]])
                repeats.append(info)
                info = []

    start = name.rfind('/')+1
    chromName = name[start:]

    output = resultDir + "/" + chromName+"Finder.bed"

    with open(output, "w") as outFile:
        
        for r in repeats:
            print r
            outFile.write(chromName)
            outFile.write("\t")

            outFile.write("\t".join(r)+"\n")





if __name__ == "__main__":

    chromFile = sys.argv[1]
    resultDir = sys.argv[2]
    name = sys.argv[3]

    print chromFile
    print resultDir
    print name

    runLTRFinder(chromFile, resultDir, name)
