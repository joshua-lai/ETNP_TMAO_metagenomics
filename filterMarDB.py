import re
import sys
import os

def parseBlastHits(blastHitFilepath):
    accessionList = []
    foundOnce = []
   
    with open (blastHitFilepath, 'r') as blastHits:
        # split the line by whitespace
        for entry in blastHits:
            line = blastHits.readline().split("\t")
            if len(line) > 1:
                accession = line[1]
                # name as fourth item, non-alphanumeric replaced w/ '_'
                notableConstant = re.search('\s.+mmp_id=', line[3]).group()
                name = re.sub('\W+','_',notableConstant)
                if name not in foundOnce:
                    foundOnce += [name]
                    accessionList += [accession]
                #can add an elif to add to a foundTwice list or something
    return(accessionList)

def writeDeduplicated(blastHitFilepath, outputFilepath, accessionList):

    with open (blastHitFilepath, 'r') as blastHits:
        with open (outputFilepath, 'w') as outfile:
            for entry in blastHits:
                line = blastHits.readline()
                if len(line) > 1:
                    accession = line.split("\t")[1]     
                    if accession in accessionList:
                        outfile.write(line)
                        


def main(args):
    ARGS_EXPECTED = 3
    name, blastHitFilepath, outputFilepath = args[0:ARGS_EXPECTED]
    for num, otherArgs in enumerate(args[ARGS_EXPECTED:]):
        print("argument",num+ARGS_EXPECTED,": \"",otherArgs,"\" not used")

    # gets a list of accessions with some important bits deduplicated
    accessionList = parseBlastHits(blastHitFilepath)
    writeDeduplicated(blastHitFilepath, outputFilepath, accessionList)

#example input
#testArgs = ['thisy', "C:/Users/joshu/Desktop/TMAO_project_files/raw_blast_hits/torAFocus/taf_marDB_e120_seq25000_blasthits_Copy.txt", "C:/Users/joshu/Desktop/TMAO_project_files/raw_blast_hits/torAFocus/taf_marDB_e120_seq25000_blasthits_shortened.txt"]

if __name__ == "__main__":
   #main(testArgs)
   main(sys.argv)