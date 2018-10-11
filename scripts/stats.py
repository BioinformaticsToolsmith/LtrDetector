import numpy as np

def getLengthStats(database):


    with open(database,'r') as infile:

        lines = infile.readlines()

        interiors = []
        ltrs= []
        repeats =set()

        i =0

        while i <len(lines)-1:

            curr = lines[i]
            
            
            if curr[0] == ">":

                label = curr.split()[0][1:]


                seq = ""
                next =""
                
                print label
                print seq
        
            while i<len(lines)-1:
                i+=1
                next = lines[i][0]

                if next == ">":
                    break
                seq+= lines[i].rstrip()


            if label[-3:]=="LTR":
                ltrs.append(len(seq))

            elif label[-1]=="I":
                interiors.append(len(seq))

    print " MAX LTR LENGTH= " +str(max(ltrs))
    print "MIN LTR LENGTH= "+str(min(ltrs))
    mean_ltr = sum(ltrs)/len(ltrs)
    print "AVERAGE LTR LENGTH= "+str(mean_ltr)
    stdev_ltr = np.std(ltrs)

    print "StdDeviation LTR = "+str(stdev_ltr)
    normal = [f for f in ltrs if (f >= mean_ltr -2*stdev_ltr and f <= mean_ltr +2*stdev_ltr)]
    print "PCT LTR +-2StdDev" +str(float(len(normal))/len(ltrs))

    print " MAX Interior LENGTH= " + str(max(interiors))
    print "MIN Interior LENGTH= "+str(min(interiors))
    mean_interior = sum(interiors)/len(interiors)

    print "AVERAGE Interior LENGTH= "+str(mean_interior)
    stdev_interior= np.std(interiors)


    print "StdDeviation Interior = "+str(stdev_interior)

    filtered = [f for f in interiors if f>=mean_interior-2*stdev_interior and f <= mean_interior+2*stdev_interior]
    print "PCT Interiors +-2StdDev " + str(float(len(filtered))/len(interiors))


            
if __name__ == "__main__":
    
    #getLengthStats("thaliana/Arabidopsis_thaliana_LTR_Retrotransposon.txt")
    #getLengthStats("sbicolor/Panicoideae_LTR_Retrotransposon.txt")
    #getLengthStats("gmax/Soybean_LTR_Retrotransposon.txt")
    #getLengthStats("irgsp/Oryza_LTR_Retrotransposon.txt")
    getLengthStats("zMays/Panicoideae_LTR_Retrotransposon.txt")
