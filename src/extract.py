def getSequence(file, start,end):

    with open(file, 'r') as infile:
        label = infile.readline()
        print label
        lines = [ l.rstrip() for l in infile.readlines()]

        i = 0
        seq =""
        begin = 0

        for j in range(len(lines)):
            l = lines[j]
            length = len(l)

            if i+length>start:
                remainder = length -abs(i+length - start)
                seq+=l[remainder-1:]
                print "i=" +str(i+remainder)
                i+=length
                
                begin = j+1
                break
                
            else:  
                i+=length

        for j in range(begin,len(lines)):
            l = lines[j]
            length = len(l)

            if i+length>end:
                remainder = length -abs(i+length - end)
                seq+=l[:remainder]
                i+=length
                begin = j+1
                break
                
            else:
                seq+=l  
                i+=length
    return seq


def riceChr1FP():

    with open("../data/irgsp/detector/FP/chr1FP.bed", "r") as inFile:

        lines = inFile.readlines()

        
        storage = []

        for l in lines:

            fields = l.split()

            s1 = int(fields[3])
            e1 = int(fields[4])

            s2 = int(fields[5])
            e2 = int(fields[6])

            left = getSequence("../data/irgsp/Fa/chr1.fa",s1,e1)

            right = getSequence("../data/irgsp/Fa/chr1.fa",s2,e2)

            info = (s1,e1,left,s2,e2,right)

            storage.append(info)

    with open("FPSequences.txt", "w") as outFile:

        for s in storage:

            outFile.write(str(s[0])+":"+str(s[1])+"=="+s[2]+"\n")
            outFile.write(str(s[3])+":"+str(s[4])+"=="+s[5]+"\n"+"\n")



if __name__ == "__main__":
    #getSequence("../data/irgsp/Fa/chr1.fa", 9585316,9585550)
    #getSequence("../data/irgsp/Fa/chr1.fa", 4011582,	4011834)
    riceChr1FP()

    #getSequence("data/thaliana/Fa/chr1.fa", 1305117, 13005214)



        
