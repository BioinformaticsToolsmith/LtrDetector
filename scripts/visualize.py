import os,sys 
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def slice(bedfile,name):

    if not os.path.exists(name):
        os.mkdir(name)

    coords = []

    with open(bedfile,'r') as infile:

        lines = infile.readlines()

        for l in lines:
            fields = l.split()
            info = (fields[1],fields[2],fields[4],fields[5])
            coords.append(info)
    
    scores =[]

    with open("/home/joseph/Projects/ltr/output/cleanedScores.csv",'r') as infile:
        scores = infile.readlines()
        

    for el in coords:



        
        print el
        start = int(el[0])
        end = int(el[1])

       

        

        mid1 = int(el[2]) 
        mid2 = int(el[3]) 
        

        bounds =[start,mid1,mid2,end-1]
        
        first = [start,mid1]
        second = [mid2,end-1]
        
        slice =name+"/"+el[0]+"-"+el[1]+".png"

        curr = scores[start:end]

        x = [int(l.split(',')[0]) for l in curr]
        y = [int(l.split(',')[1]) for l in curr]

        

        

        d = {"index":x, "score":y}
        

        df = pd.DataFrame(d,index = d["index"])
        title = "LTR("+str(start)+":"+str(end)+")"

        graph = df.plot(x = "index", y = "score",legend=False,title=title)

        graph.set_xlabel("Index")
        graph.set_ylabel("Score")

        graph.plot(bounds,[0.0]*len(bounds),"rv")



        fig = graph.get_figure()
        fig.tight_layout()

        fig.savefig(slice,dpi = 200)


if __name__ == "__main__":

    repeats = sys.argv[1]
    outFolder = sys.argv[2]
    #slice("thaliana/Test/chr5Detector.bed", "thaliana/test") 
    slice(repeats,outFolder)

        
        


