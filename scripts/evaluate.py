
import subprocess
import os
import shutil
import time
import sys



def evaluateTool(tool, directory, output_name):

    path_results = directory+"/"+tool+"/"

    if not os.path.isdir(path_results):
        os.mkdir(path_results)

    path_truth = directory + "/"+"Truth"+"/"
    path_output = path_results+output_name

    results = os.listdir(path_results)
    results = [f for f in results if f.endswith(".bed")]
    print results

    truth = os.listdir(path_truth)
    truth = [t for t in truth if t.endswith(".bed")]
    print truth

    pairs = []

    for r in results:
        for o in truth:
            slice = 0
            if tool == "detector":
                slice = -12
            elif tool == "harvest":
                slice = -11
            elif tool == "finder":
                slice = -10
            print slice
            print o[slice:]

            pre_len = len(r[:slice])

            if o.startswith(r[:slice]) and o[pre_len:] == "Truth.bed":
                name = directory+"/"+r[:slice]
                info = (name, path_results+r, path_truth+o)
                pairs.append(info)
                break

    storage = []

    for p in pairs:
        print p
        info = [p[0]]+list(runBedTools(p[2], p[1], '0.95'))
        info = [str(x) for x in info]
        storage.append(info)

    with open(path_output, 'w') as outFile:
        outFile.write("Trial,Found,Guesses,Total \n")
        for s in storage:
            outFile.write(",".join(s)+"\n")


def runBedTools(truthFile, bedFile, similarity):
    ''' Runs bedtools intersect utility to calculate the entries in truthFile that overlap with entries in bedFile above some similarity. Returns number of overlaps
    '''
    params = ['bedtools', 'intersect','-sorted', '-a', truthFile, '-b',bedFile, '-u', '-f', similarity, '-F', similarity]

    #params = ['bedops', "-e", "95%", truthFile, bedFile]
    intersect = subprocess.Popen(params, stdout=subprocess.PIPE)
    line_count = subprocess.Popen(
        ['wc', '-l'], stdin=intersect.stdout, stdout=subprocess.PIPE)
    overlap = int(line_count.communicate()[0])

    guesses = countLines(bedFile)
    total = countLines(truthFile)

    return overlap, guesses, total


def countLines(file):
    pipe = open(file, 'r')
    wc = subprocess.Popen(['wc', '-l'], stdin=pipe, stdout=subprocess.PIPE)
    line_count = int(wc.communicate()[0])
    return line_count


if __name__ == "__main__":

    tool = sys.argv[1]
    directory = sys.argv[2]
    output_name = sys.argv[3]

    evaluateTool(tool,directory,output_name)
