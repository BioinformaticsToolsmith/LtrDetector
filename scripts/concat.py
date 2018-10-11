

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
    print "AVERAGE LTR LENGTH= "+str(sum(ltrs)/len(ltrs))

    print " MAX Interior LENGTH= " + str(max(interiors))
    print "MIN Interior LENGTH= "+str(min(interiors))
    print "AVERAGE Interior LENGTH= "+str(sum(interiors)/len(interiors))

            
      



            
if __name__ == "__main__":
    
    getLengthStats("thaliana/Arabidopsis_thaliana_LTR_Retrotransposon.txt")