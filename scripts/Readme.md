### Visualization

pip is required to install dependencies. virtualenv is recommended to isolate packages.

`pip install -r requirements.txt`

You will need to run LtrDetector with the -rawScores and/or -cleanedScores flags and pass the name of a bedfile containing repeats, a destination folder, and the location of the scores file produced by LtrDetector.

`python visualize.py bedFile destination scoresFile `


## Ground Truth Generation

1. Download bedtools :
https://bedtools.readthedocs.io/en/latest/content/installation.html

2. Download NCBI BLAST+ :
https://www.ncbi.nlm.nih.gov/books/NBK52640/

3. You will need to produce your own REPEATMASKER files (extension .out) :
http://www.repeatmasker.org/

4. create a folder \<organismName\>

5. Download the LTR-RTs sequences from REPBASE (you will need an account) :
https://www.girinst.org/repbase/update/browse.php?type=LTR+Retrotransposon&format=FASTA

6. Save as \<organismName>/\<repeatFileName>


7. Place all FASTA files for an organism into a folder called \<organismName\>/Fa.

8. Place all REPEATMASKER output files in a folder \<organismName\>/Out. 

9. ` python generateTruthPipeline.py /organismName  repeatFileName `

## False Positives Evaluation

Assumes that LtrDetector results are in a file called organismName/detector, and that REPEATMASKER output has been generated.

` python correctionPipeline.py detector organismName FPSummary.txt `

## True Positives Evaluation

Assumes that LtrDetector results are in a file called organismName/detector, and that ground truth script has been run.

` python evaluate.py detector organismName TPSummary.txt `



