# adapted from
# https://github.com/edgraham/GhostKoalaParser/blob/master/GhostKOALA-taxonomy-to-anvio
# ty anvi'o and meren lab thing
#!/usr/bin/env python
import pandas as pd
import sys

def keggTaxaToTxt(inputTopFile, outputFile):
    """takes in the input for the kegg taxa .gz file thing and a filepath for an output file
    writes to the file the data in tab-delimited format"""
    with open (inputTopFile, 'r') as infile:
        with open (outputFile, 'w') as outfile:
            #.tz file to pandas i think
            GK_taxonomy = pd.read_table(infile,header=None).replace({'genecall_': ''}, regex=True)
            for index, item in enumerate(GK_taxonomy[0]):
                if item[:5] == "user:":
                    GK_taxonomy[0][index] = item[5:]
                    print(item[0:])
            #saw somewhere that to_dict then iterating is faster and i might as well get used to it
            keggRecords = GK_taxonomy.to_dict("records")
            outfile.write('taxa\tkegg_domain\tkegg_class\tkegg_species\ttaxaScore\n')
            #skip the first row, this is inelegant oh well
            #counter = False
            for row in keggRecords:
                #if counter:
                outfile.write(f"{row[0]}\t{row[2]}\t{row[3]}\t{row[4]}\t{row[6]}\n")
                #counter = True


def main(args):
    ARGS_EXPECTED = 3
    name, inputTopFilepath, outputFilepath = args[0:ARGS_EXPECTED]
    for num, otherArgs in enumerate(args[ARGS_EXPECTED:]):
        print("argument",num+ARGS_EXPECTED,": \"",otherArgs,"\" not used")
    keggTaxaToTxt(inputTopFilepath, outputFilepath)


if __name__ == "__main__":
   main(sys.argv)