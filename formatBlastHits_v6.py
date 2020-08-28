import re
import sys
import os
from Bio import SeqIO


def formatAccession(preAccession):
    """gets accesion from within the |'s or the entirety
    this was made with first seeing clara's prokdb but i think most things are fine w/o it"""
    if '|' in preAccession:
        start = preAccession.find('|')
        end = preAccession.rfind('|')
        postAccession = preAccession[start + 1:end]
    else:
        postAccession = preAccession
    return(postAccession)


def parseBlastHits(blastHitFilepath, dbName):
    """takes in a string of a filepath to a fasta formatted document
    returns three lists, one as a list of unique acession numbers
    a corresponding list of minimum e-values
    and a corresponding list of names formatted to remove non-alphanumeric characters"""
    accessionList = []
    eValList = []
    nameList = []
   
    with open (blastHitFilepath, 'r') as blastHits:
        for entry in blastHits:
            # split the line by whitespace
            line = entry.split("\t")
            if len(line) > 1:
                # sequence id as the second item, possibly surrounded by '|'s
                try:
                    if dbName == "prokdb":
                        accession = formatAccession(line[1])
                    else:
                        accession = line[1]
                except:
                    print("error in accession with line", line)
                    accession = "fake_someErrorHappened"
                    break
                # e-value as the third item
                eVal = float(line[2])
                # name as fourth item, non-alphanumeric replaced w/ '_'
                name = re.sub('[^0-9a-zA-Z]+', '_', line[3])
                
                #adds sequence id, eval, name to a list
                if accession not in accessionList:
                    accessionList += [accession]
                    eValList += [eVal]
                    nameList += [name]
                # if name already included, replaces e-value with the smaller e-value
                else:
                    indexOfMatch = accessionList.index(accession)
                    if eValList[indexOfMatch] > eVal:
                        eValList[indexOfMatch] = min(eValList[indexOfMatch], eVal)
                        nameList[indexOfMatch] = name
    return(accessionList, eValList, nameList)



def copyQuerySequences(querySequenceFilepath, accessionList, sequenceList):
    """takes in the filepath of the query sequences
    or 'na' to not add it.
    returns three lists:     
    a list of accession numbers / sequence ids,
    a list of corresponding names,
    a list of corresponding sequences"""
    accessionListQ = []
    nameListQ = []
    sequenceListQ = []
    if querySequenceFilepath !=  "na":
        with open (querySequenceFilepath, 'r') as querySeq:
            for record in SeqIO.parse(querySeq, "fasta"):
                #format accession
                accession = record.id
                #format name from description
                name = record.description.split(maxsplit = 1)[1]
                name = [re.sub('[^0-9a-zA-Z]+', '_', name)]
                
                if accession not in accessionList:
                    if record.seq not in sequenceList:
                        accessionListQ += [accession]
                        nameListQ += [name]
                        sequenceListQ += [record.seq]

    return(accessionListQ, nameListQ, sequenceListQ)



def crossReferenceSequences (blastDatabaseAsFastaFilepath, accessionList, minSeqLength, errorOutfile):
    """takes in the filepath of initial sequence that were blasted
    takes in a list of accession numbers / sequence id's
    returns a list of the corresponding sequences"""
    sequenceList = ["nothing_found"] * len(accessionList)
    numU = 0
    numDup = 0
    numFrag = 0
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
    print(accessionList)
    for fastaFile in fastaFilesToLookIn:
        with open (fastaFile, 'r') as bdFasta:
            with open (errorOutfile, 'w') as errorOutput:
                #print("opening ",bdFasta)
                for record in SeqIO.parse(bdFasta, "fasta"):
                    #format accession of record.id for comparison
                    accession = record.id
                    #print(accession, accession in accessionList)
                    if accession in accessionList:
                        #if accession in the accession list, copy the corresponding sequence
                        accessionMatchIndex = accessionList.index(accession)
                        sequenceToCopy = str(record.seq)
                        # if sequence already in the list or contains a U, adjusts it, and copies the sequence to list
                        if sequenceToCopy in sequenceList:
                            sequenceMatchIndex = sequenceList.index(sequenceToCopy)
                            sequenceList[accessionMatchIndex] = "!dupOf" + str(accessionList[sequenceMatchIndex])
                            errorOutput.write(f"{accession} is a duplicate of {accessionList[sequenceMatchIndex]}\n")
                            numDup += 1 
                        # trim short sequences 
                        elif len(sequenceToCopy) < minSeqLength:
                            errorOutput.write(f"{accession} was too short with only {len(sequenceToCopy)} aa's\n")
                            sequenceList[accessionMatchIndex] = "!tooShort" + str(len(sequenceToCopy))
                            numFrag += 1
                        elif "U" in sequenceToCopy:
                            errorOutput.write(f"selenocysteine (U) in {record}, not removed\n")
                            #print("AAHHH! why is there a 'U' in " + str(record) + ", U removed")
                            sequenceList[accessionMatchIndex] = sequenceToCopy#.replace('U','')
                            numU += 1
                        # copies sequence into sequence list
                        else:
                            sequenceList[accessionMatchIndex] = sequenceToCopy
                            errorOutput.write(f"{accession} added\n")
                errorOutput.write(f"There were {numDup} duplicates\n")
                errorOutput.write(f"There were {numFrag} sequences shorter than {minSeqLength}\n")
                errorOutput.write(f"There were {numU} sequences with a U in them\n")
    return(sequenceList)



def writeFasta(accessionList, nameList, sequenceList, eValList, outputFilepath, dbName):
    """takes in a list of accession numbers, a corresponding list of names 
    and a corresponding sequence of sequences
    takes in a filepath to write a fasta formatted file to"""
    with open (outputFilepath, 'w') as outfile:
        for ii in range(len(accessionList)):
            if sequenceList[ii][0] != "!":
                outfile.write(f">{accessionList[ii]} {nameList[ii]} {dbName}_{eValList[ii]}\n")
                temp = str(sequenceList[ii])
                outfile.write(f"{temp}\n")



def main(args):
    ARGS_EXPECTED = 8
    name, blastHitFilepath, querySequenceFilepath, blastDatabaseAsFastaFilepath, outputFilepath, errorOutfile, dbName, minSeqLength  = args[0:ARGS_EXPECTED]
    for num, otherArgs in enumerate(args[ARGS_EXPECTED:]):
        print("argument",num+ARGS_EXPECTED,": \"",otherArgs,"\" not used")

    # reads a fasta file with the blast hits
    accessionList, eValList, nameList = parseBlastHits(blastHitFilepath, dbName)
    # cross references the accessions agianst the sequences for the blast database
    sequenceList = crossReferenceSequences (blastDatabaseAsFastaFilepath, accessionList, int(minSeqLength), errorOutfile)
    # could add the original query sequences
    accessionListQ, nameListQ, sequenceListQ = copyQuerySequences(querySequenceFilepath, accessionList, sequenceList)
    accessionList += accessionListQ
    eValList += ["0.0"]*len(accessionList)
    nameList += nameListQ
    sequenceList += sequenceListQ
    # writes results and sequences to a fasta file
    writeFasta(accessionList, nameList, sequenceList, eValList, outputFilepath, dbName)
    

if __name__ == "__main__":
   main(sys.argv)