import re
import sys
import os
from Bio import SeqIO

def readPatternList(patternListInput):
    patternList =[]
    try:
        with open(patternListInput, 'r') as infile:
            for line in infile:
                patternList += [line.rstrip()]
    except:
        patternList = [str(patternListInput)]
    return(patternList)

def crossReferenceSequences (blastDatabaseAsFastaFilepath, outputFilepath, patternList):
    """takes in the filepath of initial sequence that were blasted
    takes in a list of accession numbers / sequence id's
    returns a list of the corresponding sequences"""
  
    fastaFilesToLookIn = []
    # searches through a variety of fasta files if needed, kinda assumes exclusive of each other
    # prefers later copies if not
    if blastDatabaseAsFastaFilepath[-1] == "/":
        filenames = os.listdir(blastDatabaseAsFastaFilepath)
        for file in filenames:
            if ((file[-3:] == ".fa") or (file[-4:] == ".faa")) or ((file[-6:] == ".fasta") or (file[-4:] == ".txt")):
                fastaFilesToLookIn += [blastDatabaseAsFastaFilepath + file]
    else:
        fastaFilesToLookIn = [blastDatabaseAsFastaFilepath]

    # opens fasta files to search
    with open (outputFilepath, 'w') as outfile:
        for fastaFile in fastaFilesToLookIn:
            with open (fastaFile, 'r') as bdFasta:
                for record in SeqIO.parse(bdFasta, "fasta"):
                    #format accession of record.id for comparison
                    for num,pattern in enumerate(patternList):                        
                        if pattern in record.description:
                            #outfile.write(">" + str(record.id) + " " + str(record.description) +" "+ str(num) + "\n" + str(record.seq) + "\n")
                            #outfile.write(">athL " + str(record.id) + " " + str(record.description) + "\n" + str(record.seq) + "\n")
                            outfile.write(">" + str(record.id) + " " + str(record.description) +"\n" + str(record.seq) + "\n")


def main(args):
    ARGS_EXPECTED = 4
    name, blastHitFilepath, outputFilepath, patternListInput  = args[0:ARGS_EXPECTED]
    for num, otherArgs in enumerate(args[ARGS_EXPECTED:]):
        print("argument",num+ARGS_EXPECTED,": \"",otherArgs,"\" not used")
    '''
    patternList = ["130034_kc3_L621_2012", "524620_kc3_L551_2012", "331148_kc4_L835_2012", "295231_kc5_L10221_2012", \
                    "352435_kc6_L4116_2012", "90334_kc3_L796_2012", "349037_kc3_L759_2012", "152437_kc6_L923_2012", \
                    "565808_kc28_L2041_2012","609033_kc6_L9343_2012","432009_kc7_L2013_2012","687906_kc4_L1318_2012", \
                    "286980_kc5_L1244_2012","243402_kc6_L3894_2012", "349530_kc5_L2393_2012","803036_kc3_L736_2012",\
                    "410293_kc5_L3498_2012","139786_kc17_L21529_2012","106366_kc4_L578_2012","649432_KC4_L1039_2012"]
    patternList2 = ["kc18_L19220_ETNP_all_2012","kc37_L43361_ETNP_all_2012", \
                    "kc23_L10497_ETNP_all_2012","kc8_L8918_ETNP_all_2012",\
                    "kc106_L34799_ETNP_all_2012"]
    '''
    patternList = readPatternList(patternListInput)
    print(patternList)
    # cross references the accessions agianst the sequences for the blast database
    crossReferenceSequences (blastHitFilepath, outputFilepath, patternList)

#crossReferenceSequences ("fakeGrep.py", "/c/Users/joshu/Desktop/TMAO_project_files/scripts/test_fasta.faa", "C:/Users/joshu/Desktop/TMAO_project_files/scripts/test_fasta_results.faa")
#testArgs = ["fakeGrep.py", "/C:/Users/joshu/Desktop/TMAO_project_files/scripts/test_fasta.faa", "C:/Users/joshu/Desktop/TMAO_project_files/scripts/test_fasta_results.faa"] #obselete

if __name__ == "__main__":
    #main(testArgs)
    main(sys.argv)