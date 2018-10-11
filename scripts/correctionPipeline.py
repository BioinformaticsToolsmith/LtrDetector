

import os
import subprocess
import shutil
from collections import OrderedDict
import evaluate


def findRepeats(parent,name):
    
    families = {}
    storage = []

    outFile = parent+"/Out/"+name+".fa.out"
    print outFile

    with open(outFile, 'r') as infile:

        lines = infile.readlines()

        for l in lines:

            fields = l.split()

            if len(fields) > 0:

                head = fields[0]

                if head.isdigit():

                    ltr = "LTR" in fields[10] or "LTR" in fields[9]
                    complex = fields[10] != "Simple_repeat" and fields[10] != "Low_complexity"
                    sine_line = "SINE" in fields[10] or "LINE" in fields[10]
                    dna = "DNA" in fields[10]
            
                    
                    if complex and not ltr and sine_line :
                        family = fields[9]
                        start = fields[5]
                        end = fields[6]
                        info = (start, end)

                    
                        storage.append((start,end))
                           
                        
    repeats_dir = parent +"/"+"Repeats"

    if not os.path.isdir(repeats_dir):

        os.mkdir(repeats_dir)

    repeats_file = repeats_dir+"/"+name+"Repeats.bed"
    print repeats_file

    storage.sort(cmp=lambda x, y: compare(x, y))

    with open(repeats_file,'w') as outFile:

        for s in storage:
            
            outFile.write(name+"\t"+"\t".join(s)+"\n") 


def compare(x, y):

    first = int(x[0])
    second = int(y[0])

    diff = first - second

    val = -1 if diff < 0 else 1
    if diff == 0:
        val = 0
    return val




def findPotentialFP(parent,tool,name):
    ''' removes true positives from detection '''

    detection = parent+"/"+tool+"/"+name+getFileSuffix(tool)+".bed"

    groundTruth = parent +"/"+"Truth"+"/"+name+"Truth"+".bed"

    params = ['bedtools', 'intersect', '-a', detection, '-b',groundTruth, '-v', '-f', '0.95', '-F', '0.95']

    #params = params = ['bedops', "-n", "95%", detection, groundTruth]

    result = subprocess.check_output(params)

    fp = parent+"/"+tool+"/"+name+"NoOverlap.bed"

    with open(fp,"w") as outFile:

        outFile.write(result)



def getFileSuffix(tool):

    start =tool[0].upper()

    return start+tool[1:]


def findOverlaps(parent,tool,name):
    ''' Finds overlaps between potential false positives and all repeats '''

    fp = parent+"/"+tool+"/"+name+"NoOverlap.bed"
    all_repeats = parent+"/Repeats/"+name+"Repeats.bed"

    params = ['bedtools', 'intersect',"-sorted", '-a', fp, '-b',all_repeats, '-wo']
    #params = ['bedops', '-e',"1%",fp,all_repeats]

    result = subprocess.check_output(params)

    overlaps = parent+"/"+tool+"/"+name+"Overlaps.bed"

    with open(overlaps,"w") as outFile:

        outFile.write(result)


def processOverlaps(parent,tool,name):

    overlapFile = parent+"/"+tool+"/"+name+"Overlaps.bed"

    print overlapFile

    storage = []

    with open(overlapFile, 'r') as inFile:

        lines = inFile.readlines()

        i = 0

        while i < len(lines)-1:

            curr = lines[i]
            check_for_more = True
            j = i+1

            save = [curr]

            while check_for_more:

                next = lines[j]

                print next
                print curr

                start_curr = curr.split("\t")[1]
                start_next = next.split("\t")[1]

                if start_curr != start_next or j == len(lines)-1:

                    if len(save) >1:
                        storage.append(save)

                    save = []
                    check_for_more = False
                    i = j

                else:
                    save.append(next)
                    j += 1

 
    families = recordFamilies(parent,name)
    false_positives = []

    for s in storage:
        print "PRINTING S"
        print s
        repeat_overlaps =[]

        for el in s:

            fields = el.split("\t")
            repeat_start = fields[8]
            repeat_end = fields[9]

            s1 = fields[3]
            e1 = fields[4]

            s2 = fields[5]
            e2 = fields[6]

            if intersect(repeat_start,repeat_end,s1,e1,0.8):
                repeat_overlaps.append(repeat_start)

            elif intersect(repeat_start,repeat_end,s2,e2,0.8):
                repeat_overlaps.append(repeat_start)
        
        origins = [families[f] for f in repeat_overlaps if f in families]
        print origins

        if len(origins)>=2:

            false_positives.append(fields[:7])

    false_positives_dir = parent+"/"+tool+"/"+"FP"

    if not os.path.isdir(false_positives_dir):
        os.mkdir(false_positives_dir)

    false_positives_file = false_positives_dir+"/"+name+"FP.bed"


    with open(false_positives_file, "w") as outFile:

        for f in false_positives:
            outFile.write("\t".join(f)+"\n")


            
def intersect(start1,end1,start2,end2,pct):

    s1 = int(start1)
    e1 = int(end1)
    len1 = e1-s1

    s2 = int(start2)
    e2 = int(end2)
    len2 = e2-s2

    a_in_b = s2<=s1 and s1<=e2
    b_in_a = s1<=s2 and s2<=e1

    overlap = 0

    if a_in_b:

        overlap = e2 - s1
        overlap = len1 if overlap > len1 else overlap
    elif b_in_a:

        overlap = e1 -s2
        overlap = len2 if overlap > len2 else overlap

    a_overlap = float(overlap)/len1
    b_overlap = float(overlap)/len2

    return a_overlap >= pct and b_overlap >= pct 



def countLines(file):
    pipe = open(file, 'r')
    wc = subprocess.Popen(['wc', '-l'], stdin=pipe, stdout=subprocess.PIPE)
    line_count = int(wc.communicate()[0])
    pipe.close()
    return line_count
    
def recordFamilies(parent,name):
    '''returns dictionary mapping starting indexes of repeats to the family of repeat that it belongs to'''

    families = {}
    storage = []

    outFile = parent+"/Out/"+name+".fa.out"

    with open(outFile, 'r') as infile:

        lines = infile.readlines()

        for l in lines:

            fields = l.split()

            if len(fields) > 0:

                head = fields[0]

                if head.isdigit():

                  
                    ltr = "LTR" in fields[10] or "LTR" in fields[9]
                    complex = fields[10] != "Simple_repeat" and fields[10] != "Low_complexity"
                                  
                    if complex and not ltr:

                        family = fields[9]
                        start = fields[5]
                        end = fields[6]
                        info = (start, end)

                        families[start] = family
    return families
   
def deleteFiles(parent,tool,name):
    no_overlap = parent+"/"+tool+"/"+name+"NoOverlap.bed"
    overlap = parent+"/"+tool+"/"+name+"Overlaps.bed"
    os.remove(no_overlap)
    os.remove(overlap)


def idFalsePositives(parent,tool,name):

    findPotentialFP(parent,tool,name)
    print name
    findRepeats(parent,name)
    findOverlaps(parent,tool,name)
    processOverlaps(parent,tool,name)

def countFalsePositives(parent,tool,summary):

    out_dir = parent+"/Out"
    outs = os.listdir(out_dir)

    outs = [o[:-7] for o in outs if o.endswith("fa.out")]

    storage = []

    for o in outs:
        idFalsePositives(parent,tool,o)

        result = parent+"/"+tool+"/FP/"+o+"FP.bed"

        fp = str(countLines(result))

        deleteFiles(parent,tool,o)

        info = (o,fp)
        storage.append(info)

    summary = parent+"/"+tool+"/"+summary
    
    with open(summary,'w') as outFile:

        for s in storage:

            outFile.write(",".join(s)+"\n")



if __name__ == "__main__":
    #countFalsePositives("sbicolor", "detector", "summaryFP.txt")
    #countFalsePositives("barley","detector","summaryFP_9_13.txt")
    #countFalsePositives("thaliana","detector", "summary_9_12.txt")
    #countFalsePositives("zMays","detector","summaryFP.txt")
    countFalsePositives("irgsp","detector","summaryFP_9_12.txt")
