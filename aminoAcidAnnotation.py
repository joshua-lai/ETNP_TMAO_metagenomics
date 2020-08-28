#aminoAcidAnnotation.py
import sys
import re
import numpy as np
from Bio import SeqIO

def positionFinder(alignmentFastaFilepath, referenceRecord, segmentPattern, chars):
    """takes in the filepath of a fasta file
    and gets a particular record, to find a particular segment
    and where it is in an alignment"""
    with open (alignmentFastaFilepath, 'r') as fastaInput:
        records = SeqIO.parse(fastaInput, "fasta")
        rCapDorASeq = str(SeqIO.to_dict(records)[referenceRecord].seq)
        #the alignment has dashes but this pattern should search for it regardless of intermittant dashes
        #doesn't work great is part fo this seqeuence is missing but hopefully it shouldn't go away
        dashedSegment = re.search(segmentPattern, rCapDorASeq).group()
        #gets the position of the pattern
        segPos = rCapDorASeq.find(dashedSegment)
        #finds where certain postions are in the sequence to tell how much they are offset even with the dashes
        posList = []
        for char in chars:
            charOffset = dashedSegment.find(char)
            charPosition = segPos + charOffset
            posList += [charPosition]

        #things for testing matches
        """print(rCapDorASeq)
        print(dashedSegment)
        print(segPos)
        print(rCapDorASeq[segPos-12:segPos+12])"""
        return(posList)            

def writeAaAnnotation(alignmentFastaFilepath, yPos, wPos, uPos, periPos, annotationFilepath):
    """records to an tab-delimited file what is in particular columns"""
    with open (alignmentFastaFilepath, 'r') as fastaInput:
        with open (annotationFilepath, 'w') as outputFilepath:            
            outputFilepath.write("taxa\tisY\tisW\tisYW\tuPos\tnumc1\tnumc2\tperiSecMotifMatch\n")
            for record in SeqIO.parse(fastaInput, "fasta"):
                
                recString = str(record.seq)
                isY = record.seq[yPos] == 'Y'
                isW = record.seq[wPos] == 'W'
                isYW = isY and isW
                uPosAA = record.seq[uPos]
                # certain motifs i found relating to fe-s clusters
                numClusterOnes = len(re.findall("[HC](-*[^-C]){2}C(-*[^-C]){3}C(-*[^-C]){15,45}C", recString))
                numClusterTwos = len(re.findall("[HC](-*[^-C]){3}C(-*[^-C]){3}C(-*[^-C]){15,45}C", recString))
                # matches in the twin arg pattern
                numPeriMatch = 0
                numPeriMatch += recString[periPos + 0] == 'S' or recString[periPos] == 'T'
                numPeriMatch += recString[periPos + 1] == 'R'
                numPeriMatch += recString[periPos + 2] == 'R'
                numPeriMatch += recString[periPos + 4] == 'F'
                numPeriMatch += recString[periPos + 5] == 'L'
                numPeriMatch += recString[periPos + 6] == 'K'
                outputFilepath.write(f"{record.id}\t{isY}\t{isW}\t{isYW}\t{uPosAA}\t{numClusterOnes}\t{numClusterTwos}\t{numPeriMatch}\n")

                
def main(args):
    #takes in arguments from command line
    ARGS_EXPECTED = 3
    name, alignmentFastaFilepath, annotationFilepath  = args[0:ARGS_EXPECTED]
    for num, otherArgs in enumerate(args[ARGS_EXPECTED:]):
        print("argument",num+ARGS_EXPECTED,": \"",otherArgs,"\" not used")

    #reads the original filepath to get the position of y114 and w116 and the twin arginine motif
    yPos, wPos = positionFinder(alignmentFastaFilepath, "sp|Q52675|DSTOR_RHOCA", "G-*G-*S-*Y-*G-*W-*K-*S-*P-*G", ('Y','W'))
    uPos = positionFinder(alignmentFastaFilepath, "sp|P32176|FDOG_ECOLI", "D-*N-*Q-*A-*R-*V-*[XU]-*H-*G-*P-*T", ('U'))[0]
    # i probably need a mroe elegant way to do this because right now it basically assumes no dash gaps in the twin-arg pattern
    periPos = positionFinder(alignmentFastaFilepath, "sp|P81186.2|NAPA_DESDA", "[ST]-*R-*R-*.-*F-*L-*K", ('S'))[0]

    #writes a fasta file to a filepath
    writeAaAnnotation(alignmentFastaFilepath, yPos, wPos, uPos, periPos, annotationFilepath)
 
if __name__ == "__main__":
   main(sys.argv)
     

