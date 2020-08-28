import matplotlib.pyplot as plt
import sys
import numpy as np
from Bio import SeqIO
from collections import Counter

def alignmentFastaParser(alignmentFastaFilepath, minLength, errorOutputFilepath):
    """takes in the filepath of a fasta file
    returns three lists with id's descriptions and sequences"""
    recordIDList = []
    recordDesList = []
    recordSeqList = []
    nonDashCountList = []
    with open (alignmentFastaFilepath, 'r') as fastaInput:
        with open (errorOutputFilepath, 'a+') as errorOutput:
            for record in SeqIO.parse(fastaInput, "fasta"):
                
                nonDashCount = len(record.seq) - str(record.seq).count('-') 
                nonDashCountList += [nonDashCount]
                if nonDashCount >= minLength:
                    recordIDList += [record.id]
                    recordDesList += [record.description]
                    recordSeqList += [record.seq]
                else:
                    errorOutput.write(record.id + " was removed for only having " +str(nonDashCount) + " characters\n")
    plt.hist(nonDashCountList, bins = 100)
    plt.show()
    return(recordIDList,recordDesList, recordSeqList, nonDashCountList)

def maskSequenceByDashes(recordSeqList, errorOutputFilepath, minColAgreement, maxSeqDashContent):
    """takes in a list of aligned sequences, and thresholds for the maximum number of dashes
    in the sequences to be kept. lower numbers are more restrictive.
    returns two numpy arrays with sequences and sequence positions to keep 
    """
    #gets the number of sequences (numRows) and length of sequence (numCols)
    numRows = len(recordSeqList)
    numCols = len(recordSeqList[0])
    

    #populates a numpy array with the sequences
    testStringArray = np.empty([numRows,numCols], dtype='str')
    #initialize a dictionary of protein names, with values of lists of positions
    barChartDict = {}
    for aaLetter in ["Q","W","E","R","T","Y","U","I","O","P", \
                     "A","S","D","F","G","H","J","K","L", \
                     "Z","X","C","V","B","N","M","-"]:
        barChartDict[aaLetter] = [0]*numCols

    for rowNum, string in enumerate(recordSeqList):
        for colNum, char in enumerate(string):
            testStringArray[rowNum, colNum] = char
            # for the aaLetter of the char, adds one to the col num position
            try: 
                barChartDict[char][colNum] += 1
            except:
                barChartDict[char.upper()][colNum] += 1
    #creates a bool array by if a certain sequence position has fewer dashes than the threshold percent
    colsToKeep = np.full(numCols, False, dtype='bool')
    with open (errorOutputFilepath, 'a+') as errorOutput:
        for colNum in range(numCols):
            maxKey = ""
            maxKeyNum = 0
            for key in barChartDict.keys():
                if key != '-':
                    if barChartDict[key][colNum] >= minColAgreement:
                        colsToKeep[colNum] = True
                    if barChartDict[key][colNum] >= maxKeyNum:
                        maxKeyNum = barChartDict[key][colNum]
                        maxKey = key
            errorOutput.write(str(colNum) + ", " + maxKey + "(" + str(maxKeyNum) + ")\n")
    #creates a bool array by if a certain sequence has fewer dashes among kept positions than a threshold percent
    seqsToKeep = np.empty(numRows, dtype = 'bool')
    for rowNum in range(numRows):
        numDashes = np.count_nonzero(testStringArray[rowNum,colsToKeep] == '-')
        seqsToKeep[rowNum] = (numDashes <= maxSeqDashContent * sum(colsToKeep) / 100)

    '''
    #make a barChart
    noneBelow = True
    below = ''
    r = range(numCols)
    for key in barChartDict.keys():
        randomColor = tuple(np.random.uniform(0,1, size=3))
        if noneBelow:
            plt.bar(r, barChartDict[key], color = randomColor)
            below = barChartDict[key]
        else:
            plt.bar(r, barChartDict[key], color = randomColor, bottom = below)
            below += barChartDict[key]
    plt.legend(loc="upper right")
    plt.show()
    '''


    return(colsToKeep,seqsToKeep)


def writeFasta(recordIDList, recordDesList, recordSeqList, colsToKeep, seqsToKeep, outputFilepath, errorOutputFilepath):
    """takes in info on the record id, description, sequences and which part of it to keep
    writes a fasta file to the output filepath given"""
    with open (outputFilepath, 'w') as outfile:
        with open (errorOutputFilepath, 'a+') as errorOutput:
            for seqNum, seq in enumerate(recordSeqList):
                seq = np.array(seq)
                if seqsToKeep[seqNum]:
                    outfile.write(f">{recordIDList[seqNum]} {recordDesList[seqNum].replace(' ','_')}\n")
                    #outfile.write(">" + str(recordDesList[seqNum]).replace(" ","_") + "\n")
                    seqToWrite =''.join(seq[colsToKeep])
                    seqToWrite = seqToWrite.replace("U","-")
                    outfile.write(seqToWrite + "\n\n")    
                else:
                    numEmpty = np.count_nonzero(seq[colsToKeep] == '-')
                    totCols = sum(colsToKeep)
                    percentEmpty = round((numEmpty * 100 / totCols),2)
                    errorOutput.write("removed " + str(recordIDList[seqNum]) + " for having " +  \
                        str(numEmpty) + " empty positions out of " +              \
                        str(totCols) + "  (" + str(percentEmpty)+ "%)")
                
def main(args):
    #takes in arguments from command line
    ARGS_EXPECTED = 7
    name, alignmentFastaFilepath, outputFilepath, errorOutputFilepath, minColAgreement, maxSeqDashContent, minLength  = args[0:ARGS_EXPECTED]
    for num, otherArgs in enumerate(args[ARGS_EXPECTED:]):
        print("argument",num+ARGS_EXPECTED,": \"",otherArgs,"\" not used")
    with open (errorOutputFilepath, 'a+') as errorOutput:
        errorOutput.write("\n"+str(args)+"\n")
    #reads the original filepath
    recordIDList, recordDesList, recordSeqList, nonDashCountList = alignmentFastaParser(alignmentFastaFilepath, int(minLength), errorOutputFilepath)
    #figures out which sequences and parts of sequences to keep based on dashes
    colsToKeep,seqsToKeep = maskSequenceByDashes(recordSeqList, errorOutputFilepath, int(minColAgreement), float(maxSeqDashContent))
    #writes a fasta file to a filepath
    writeFasta(recordIDList, recordDesList, recordSeqList, colsToKeep, seqsToKeep, outputFilepath, errorOutputFilepath)


if __name__ == "__main__":
   main(sys.argv)
     