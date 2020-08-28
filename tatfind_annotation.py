#aminoAcidAnnotation.py
import sys
import re
import numpy as np

def tatfindReader(input, output):
    with open(input, 'r') as infile:
        with open(output, 'w') as outfile:
            outfile.write(f"taxa\thasTat\n")
            for row in infile:
                if row[0:7] == "Results":
                    taxa, hasTat = re.findall("Results for (.*?): (.*?)\n", row)[0]
                    hasTat = hasTat == "TRUE"
                    outfile.write(f"{taxa}\t{hasTat}\n")

def main(args):
    #takes in arguments from command line
    ARGS_EXPECTED = 3
    name, tatfindDocumentFilepath, tatfindOutputFilepath  = args[0:ARGS_EXPECTED]
    for num, otherArgs in enumerate(args[ARGS_EXPECTED:]):
        print("argument",num+ARGS_EXPECTED,": \"",otherArgs,"\" not used")
    tatfindReader(tatfindDocumentFilepath, tatfindOutputFilepath)
    


if __name__ == "__main__":
   main(sys.argv)
     

