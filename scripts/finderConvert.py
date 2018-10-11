
def finderConvert(finderFile,bedFile,name):

with open(finderFile, 'r') as inFile:

    lines = inFile.readlines()

    repeats = []
    info = []
    for l in lines:
        fields = l.split(' ')

        if fields[0] == "5'-LTR":
            info.extend([fields[2],fields[4]])
        elif fields[0] == "3'-LTR":
            info.extend([info[0],fields[2],fields[4],info[1]])
            repeats.append(info)
            info = []

with open(bedFile, "w" ) as outFile:
    
    outFile.write()
    for r in repeats:

            

