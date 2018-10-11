import subprocess,os,shutil

def queryDatabase(database,query,name):
    ''' database = FASTA formatted chromosome
        query = FASTA formatted representative sequences
    '''
    dbDirectory = "dbTemp"

    if not os.path.isdir(dbDirectory):
        os.mkdir(dbDirectory)
        
    out = dbDirectory+ "/"+name

    db_params = ['makeblastdb','-in',database,'-dbtype','nucl','-out',out]
    subprocess.call(db_params)
    identity = '70'
    format = "7 qseqid sseqid qstart qend sstart send length pident qcovhsp qcovs"
    blast_params =['blastn','-db',out,'-query',query,'-perc_identity', identity,'-outfmt',format]
    locAlign = subprocess.check_output(blast_params)
    
    output = name+".blast"

    with open(output,'w') as outfile:
        outfile.write(locAlign)

    shutil.rmtree(dbDirectory)

if  __name__ =="__main__":
    queryDatabase("thaliana/Fa/chr2.fa","thaliana/thaliana_totals.fasta", "thaliana/ATChr2")
    queryDatabase("thaliana/Fa/chr2.fa","thaliana/Arabidopsis_thaliana_LTR_Retrotransposon.txt", "thaliana/ATChr2LTR")
