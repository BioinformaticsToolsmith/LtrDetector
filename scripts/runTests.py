import subprocess,os,shutil,time,sys



def runTool(tool,directory):

    path_fa = directory+"/"+"Fa"+"/"
    path_destination = directory+"/"+tool+"/"

    if not os.path.isdir(path_destination):
        os.mkdir(path_destination)

    storage = []

    fa = os.listdir(path_fa)
    fa = [f for f in fa if f.endswith(".fa")]
    print fa

    instructions =[]

    for fa_file in fa:
        print fa_file
        name = directory+"/" + os.path.splitext(fa_file)[0]
        print "testing "+ name

        if tool == "harvest":
            path = path_fa+fa_file
            commands = " ".join(['python',"runHarvest.py",path,path_destination,name])
            instructions.append (commands)
        elif tool == "finder":
            path = path_fa+fa_file
            commands = " ".join(['python', "runFinder.py",path, path_destination, name])
            instructions.append(commands)
    
    instruction_file = "temp.txt"
    log_file = directory+"/"+tool+"/"+tool+"_"+directory+".txt"

    print instructions

    with open(instruction_file, 'w') as outFile:

        for i in instructions:
            outFile.write(i+"\n")

    bulkRunTool(instruction_file,log_file)

def bulkRunTool(instruction_file,log_file):
    output = subprocess.check_output(['/usr/bin/time','parallel', '-j3','-a', instruction_file])

    with open(log_file, "w") as outFile:
        outFile.write(output)


if __name__=="__main__":
    #runLTRHarvest("synthetic/zeroTest.fasta","zero")
    #runLTRHarvest("synthetic/fiveTest.fasta","five")
    #runLTRHarvest("synthetic/tenTest.fasta","ten")
    #runLTRHarvest("synthetic/fifteenTest.fasta","fifteen")
    #print(runLTRHarvest("synthetic/fiveTest.fasta", "synthetic/fiveTestTruth.bed", "synthetic/harvest"))
    #runLTRHarvest("synthetic/thirtyTest.fasta","thirty")
    #test("synthetic","summary.txt")

    tool = sys.argv[1]
    genome = sys.argv[2]
    print tool
    print genome
    runTool(tool,genome)
    
    #runTool("finder","thaliana")
    #runTool("finder","irgsp")
    #runTool("finder","gmax")
    
    #runTool("finder","synthetic")

    
    #test("finder", "thaliana","summary.txt")

