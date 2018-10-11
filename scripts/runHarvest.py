import subprocess
import os
import shutil
import time
import sys



def runLTRHarvest(chromFile, resultDir, name):

    tempDir = name
    temp_name = 'temp'

    if not os.path.isdir(tempDir):
        os.mkdir(tempDir)

    index = tempDir+"/"+temp_name+".index"

    runSuffixerator(chromFile, index)

    params = ['sudo', '../tools/genometools-1.5.9/bin/gt','ltrharvest', '-index', index]
    #,'-minlenltr',args[0],'-maxlenltr',args[1],'-mindistltr',args[2],'-maxdistltr',args[3]
    print params
    output = subprocess.check_output(params).split('\n')

    convertHarvest(output, resultDir, name)
    shutil.rmtree(tempDir)


def convertHarvest(lines, resultDir, name):

    repeats = []

    if not os.path.isdir(resultDir):
        os.mkdir(resultDir)

    for line in lines:

        fields = line.rstrip().split()
        print fields

        if len(fields) >= 1 and fields[0] != "#":
            save = (fields[0], fields[1],fields[3],fields[4],fields[6],fields[7])
            repeats.append(save)

    start = name.rfind('/')+1
    chromName = name[start:]

    output = resultDir +"/"+ chromName+"Harvest.bed"

    with open(output, "w") as outFile:

        for r in repeats:
            outFile.write(chromName)
            outFile.write("\t")

            outFile.write("\t".join(r))
            outFile.write("\n")


def runSuffixerator(chromFile, name):

    params = ['sudo', '../tools/genometools-1.5.9/bin/gt', 'suffixerator', '-db',
              chromFile, '-indexname', name, '-tis', '-suf', '-lcp', '-des', '-ssp', '-sds', '-dna']

    subprocess.call(params)


if __name__ == "__main__":

    chromFile = sys.argv[1]
    resultDir = sys.argv[2]
 

    name = sys.argv[3]

    
    print name
    
    runLTRHarvest(chromFile,resultDir,name)
