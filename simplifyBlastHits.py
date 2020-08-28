import re
import sys
import os
from Bio import SeqIO


def crossReferenceSequences (blastHitFilepath, outputFilepath, numCounts):
    """takes in the filepath of initial sequence that were blasted
    takes in a list of accession numbers / sequence id's
    returns a list of the corresponding sequences"""
  
    # opens fasta files to search
    with open (outputFilepath, 'w') as outfile:
        with open (blastHitFilepath, 'r') as blastHits:
            contigList = []
            count = 0
            for entry in blastHits:

                if count % numCounts == 0:
                    # split the line by whitespace
                    line = blastHits.readline().split("\t")
                    query = line[0]
                    splitpoint = query.find("_")
                    contigName = query[splitpoint+1:]
                    contigPartNum = query[:splitpoint]

                    if contigName not in contigList:
                        outfile.write("\n" + contigName + "\n")
                        contigList += [contigName]
                    outfile.write("\t" + contigPartNum + ":\n") 
                
                if (count % numCounts) in [1,2,3,4,5]:
                    line = blastHits.readline().split("\t")
                    name = line[3]
                    outfile.write("\t\t" + line[3] + "\t" + line[2])
                count += 1

                


def main(args):
    ARGS_EXPECTED = 3
    name, blastHitFilepath, outputFilepath  = args[0:ARGS_EXPECTED]
    for num, otherArgs in enumerate(args[ARGS_EXPECTED:]):
        print("argument",num+ARGS_EXPECTED,": \"",otherArgs,"\" not used")
    numCounts = 30
    # cross references the accessions agianst the sequences for the blast database
    crossReferenceSequences (blastHitFilepath, outputFilepath, numCounts)

#crossReferenceSequences ("fakeGrep.py", "/c/Users/joshu/Desktop/TMAO_project_files/scripts/test_fasta.faa", "C:/Users/joshu/Desktop/TMAO_project_files/scripts/test_fasta_results.faa")
#testArgs = ["fakeGrep.py", "/C:/Users/joshu/Desktop/TMAO_project_files/scripts/test_fasta.faa", "C:/Users/joshu/Desktop/TMAO_project_files/scripts/test_fasta_results.faa"]

if __name__ == "__main__":
    #main(testArgs)
    main(sys.argv)