# LtrDetector
 LtrDetector is an accurate and efficient software tool for de-novo detection of Long Terminal Repeat Retro-transposons. It is currently available for Unix/Linux/MacOS. Academic papers should cite 

## Installation

`git clone ourRepo ` 

`cd LtrDetector `

Run start.py to set up the proper directory structure. Requires Python 2.7

` python start.py `

Compilation requires GNU Make.

` make tr `

Compilation is sucessfull if the binary LtrDetector exists in /bin

## Usage

LtrDetector is a command line tool. The only required arguments are a chromosome directory and an output directory. 

```
bin/LtrDetector -chromDir <pathToChromosomeDirectory> -destDir <pathToOutputDirectory> 

```

Other arguments can be invoked in the form: -arg val 

### Optional arguments
| -arg     | Description | Default |
| ---------------- | ----------- | ------- |
| -minLen  | Minimum length of complete LTR-RTs. Constrains scoring system and filters. | 2000 |
| -maxLen  |  Maximum length of complete LTR-RTs. Constrains scoring system and filters. | 18000 |
|-minLenLTR | Minimum length of LTR direct repeat. Affects length filter. | 100 |
|  -maxLenLTR   | Maximum length of LTR direct repeat. Affects length filter. ( \*Note\* Runtime is highly dependent on this parameter, as it provides an upper bound for alignment length in the boundary adjustment step. Any value larger than the default will be set to the default)  | 2000 |                
| -id | Minimum identity [0-100] between 5' and 3' LTRs. | 85 |
| -k  | Length of k-mers to adjust scoring system. Tradeoff between noise and resistance to mutation.| 14 |
| -plateauSeed | Minimum length of plateaus to be initially considered 'Keep' in merging step. | 10 |
| -gapTol  | Number of base pairs that two plateaus can differ by in height/distance. Affects both plateau merging and pairing steps. | 200 |

### Flags
-rawScores prints the raw scores to a file rawScores.txt 
-cleanedScores prints the scores after merging to a file cleanedScores.txt 



## Evaluation scripts

You will need to download bedtools and NCBI blastn to run the ground truth generation pipeline.

pip is required to install dependencies. virtualenv is recommended to isolate packages.

`pip install -r requirements.txt`

###groundTruth