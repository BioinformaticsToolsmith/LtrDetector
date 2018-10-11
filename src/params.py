import subprocess, os

def runLTRCollector(i,params):
    
    subprocess.call(params)
    
    out = open('output/trialResults.csv', 'a')
    
    cols = []
    
    cols.append(str(i))
    for x in params[3:]:
        cols.append(str(x))
    
    bedFile = params[2]

    pipe =open(bedFile,'r')  
    wc = subprocess.Popen(['wc','-l'],stdin=pipe,stdout=subprocess.PIPE)
    line_count = int(wc.communicate()[0])
    pipe.close()
  
    cols.append(str(line_count))
    
    
    found = runBedTools('data/fifteenTestTruth.bed',bedFile,'0.95')
    cols.append(str(found))
    

    line = ','.join(cols)
    
    out.write(line)
    print line
    out.write("\n")
    out.close()
    os.remove(bedFile)
    

def runBedTools(truthFile,bedFile,similarity):
    ''' Runs bedtools intersect utility to calculate the entries in truthFile that overlap with entries in bedFile above some similarity. Returns number of overlaps
    '''
    params = ['bedtools','intersect','-a',truthFile,'-b',bedFile,'-u','-f',similarity,'-F',similarity]
    intersect = subprocess.Popen(params,stdout=subprocess.PIPE)
    line_count = subprocess.Popen(['wc','-l'], stdin = intersect.stdout,stdout=subprocess.PIPE)
    overlap = line_count.communicate()[0]
    return int(overlap)

if __name__ == "__main__":
    
    out = open('output/trialResults.csv', 'w')
    cols ="trial,k-mer size,minimum distance, maximum distance,minimum plateau length,difference tolerance,gap tolerance,candidates found,overlaps at 95%,total \n"
    out.write(cols) #first row is column names
    out.close()
    #params = ["../bin/TestTr","../src/FlyBaseR5.48/chrX-chromosome.fasta","../src/output/test1.bed",'14','11000','2000','100','10','25']
    #runLTRCollector(params)

    k_val = [10,11,12,13,14]
    min_plateau_len = [2,5,8,10,12,15]
    diff_tolerance = [50,100,150,200,500]
    gap_tolerance =  [40,50,80,90,100,110,120,150,200,300,500]
    
 
    ltrBinary = "../bin/TestTr"
    chromFile = "../src/data/fifteenTest.fasta"
    bedFilePrefix = "../src/output/bed/"
    
    #os.mkdir(bedFilePrefix)
    
    print("Testing....")
    i = 1
    
    for a in k_val:

        for d in min_plateau_len:
            for e in diff_tolerance:
        
                for f in gap_tolerance:
                    bedFile = bedFilePrefix+"Trial"+str(i)+".bed"
                    params = [ltrBinary,chromFile,bedFile,str(a),'2000','8773',str(d),str(e),str(f)]
                    try:
                        
                        runLTRCollector(i,params)
                        i+=1
                    except Exception:
                        continue
                    
    print("Test completed!Results available at output/trialResults.csv")
    
    