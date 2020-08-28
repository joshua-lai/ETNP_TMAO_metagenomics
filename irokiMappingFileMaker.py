import sys
import pandas

    # branch_width, branch_color \
    # , bar2_height, bar2_color \
    # , bar
    # , bar3_height, bar3_color
anoDict = \
    {"etnpall":"\t3\tgreen",\
     "satl":"\t3\tdarkcyan",\
     "query":"\t5\tblack",\
     "marref":"\t1\tsaddlebrown",\
     "mardb":"\t1\tbrown",\
     "prokdb":"\t1\tsandybrown",\
     "from_ETNP_all":"\t1\tgreen",\
     "from_satl":"\t1\tdarkcyan",\
     "Nitrate_reductase":"\t1\tlavender",\
     "unknown/other/recheck":"\t1\tgrey",\
     "DMSO_reductase":"\t1\tcrimson",\
     "TMAO_reductase(torA)":"\t1\tpurple",\
     "Something_Molybdopterin":"\t1\tdarkblue",\
     "Something_Molydopterin":"\t1\tdarkblue",\
     "Something_BiotinSulfoxide":"\t1\torange",\
     "Formate_reductase":"\t1\tblue",\
     "Pyrogallol_hydroxytransferase":"\t1\tyellow",\
     "torA; trimethylamine-N-oxide reductase (cytochrome c) [EC:1.7.2.3]":"\t1\tpurple",\
     "torZ; trimethylamine-N-oxide reductase (cytochrome c) [EC:1.7.2.3]":"\t1\tdarkviolet",\
     "bisC; biotin/methionine sulfoxide reductase [EC:1.-.-.-]":"\t1\torange",\
     "dmsA; anaerobic dimethyl sulfoxide reductase subunit A [EC:1.8.5.3]":"\t1\tdarkgoldenrod",\
     "narG; nitrate reductase / nitrite oxidoreductase, alpha subunit [EC:1.7.5.1 1.7.99.-]":"\t1\tlavender",\
     "nasA; assimilatory nitrate reductase catalytic subunit [EC:1.7.99.-]":"\t1\tpink",\
     "napA; nitrate reductase (cytochrome) [EC:1.9.6.1]":"\t1\tcoral",\
     "narB; ferredoxin-nitrate reductase [EC:1.7.7.2]":"\t1\tpeachpuff",\
     "fdoG; formate dehydrogenase major subunit [EC:1.17.1.9]":"\t1\taquamarine",\
     "fdhF; formate dehydrogenase (acceptor) [EC:1.17.99.7]":"\t1\tsteelblue",\
     "fdnG; formate dehydrogenase-N, alpha subunit [EC:1.17.5.3]":"\t1\troyalblue",\
     "phsA; thiosulfate reductase / polysulfide reductase chain A [EC:1.8.5.5]":"\t1\tolive",\
     "aoxB; arsenite oxidase large subunit [EC:1.20.2.1 1.20.9.1]":"\t1\tdarkred",\
     "cysJ; sulfite reductase (NADPH) flavoprotein alpha-component [EC:1.8.1.2]":"\t1\tsalmon",\
     "-":"\t0\twhite",\
     "True":"\t1\tblack",\
     "False":"\t0\twhite"\
    }




def writeMappingFile(totalDF, irokiAnnotationOutputPathway):
    with open (irokiAnnotationOutputPathway, 'w') as outputFile:
        outputFile.write("name\tbranch_width\tbranch_color\tbar1_height\tbar1_color\tbar2_height\tbar2_color\tbar3_height\tbar3_color\tbar4_height\tbar4_color\n")
        for index, row in totalDF.iterrows():
            for num, columnSet in enumerate(row):
                if num == 0:
                    outputFile.write(columnSet)
                else:
                    outputFile.write(f"{anoDict[str(columnSet)]}")
            outputFile.write("\n")

def main(args):
    #takes in arguments from command line
    ARGS_EXPECTED = 5
    #name, keggAnnotationPathway  = args[0:ARGS_EXPECTED]
    name, keggAnnotationPathway, wordSearchAnnotationPathway, aaAnnotationPathway, irokiAnnotationOutputPathway  = args[0:ARGS_EXPECTED]
    for num, otherArgs in enumerate(args[ARGS_EXPECTED:]):
        print("argument",num+ARGS_EXPECTED,": \"",otherArgs,"\" not used")
    keggDF = pandas.read_csv(keggAnnotationPathway, sep='\t', usecols =["taxa","k_number_annotation"])
    wordSearchDF = pandas.read_csv(wordSearchAnnotationPathway, sep='\t', usecols =["taxa","function", "database"])
    aaDF = pandas.read_csv(aaAnnotationPathway, sep='\t', usecols =["taxa","isY","isW"])
    
    partialDF = pandas.merge(keggDF, wordSearchDF, on=["taxa"])
    totalDF = pandas.merge(partialDF, aaDF, on=["taxa"])
    
    writeMappingFile(totalDF, irokiAnnotationOutputPathway)



if __name__ == "__main__":
   main(sys.argv)