import sys
import pandas

anoColorDict = \
    {"etnpall":"#a1f0c6",\
     "etnppool":"#166924",\
     "satl":"#1fddff",\
     "query":"#000000",\
     "marref":"#cda441",\
     "mardb":"#a4947b",\
     "prokdb":"#806b76",\
     "from_ETNP_all":"#a1f0c6",\
     "from_ETNP_pool":"#166924",\
     "from_satl":"#1fddff",\
     "Nitrate_reductase":"#2250cc",\
     "unknown/other/recheck":"#ffffff",\
     "DMSO_reductase":"#421682",\
     "TMAO_reductase(torA)":"#421682",\
     "Something_Molybdopterin":"#78828b",\
     "Something_Molybdopterin":"#78828b",\
     "Something_Dehydrogenase":"#7a576c",\
     "Something_BiotinSulfoxide":"#ff7900",\
     "Something_amidotransferase":"#4A4A4A",\
     "Formate_reductase":"#698701",\
     "Pyrogallol_hydroxytransferase":"#edcaef",\
     "torA; trimethylamine-N-oxide reductase (cytochrome c) [EC:1.7.2.3]":"#421682",\
     "torZ; trimethylamine-N-oxide reductase (cytochrome c) [EC:1.7.2.3]":"#620062",\
     "bisC; biotin/methionine sulfoxide reductase [EC:1.-.-.-]":"#ff7900",\
     "dmsA; anaerobic dimethyl sulfoxide reductase subunit A [EC:1.8.5.3]":"#ed64ea",\
     "narG; nitrate reductase / nitrite oxidoreductase, alpha subunit [EC:1.7.5.1 1.7.99.-]":"#2250cc",\
     "narG, narZ, nxrA; nitrate reductase / nitrite oxidoreductase, alpha subunit [EC:1.7.5.1 1.7.99.-]":"#2250cc",\
     "nasA; assimilatory nitrate reductase catalytic subunit [EC:1.7.99.-]":"#80e1f2",\
     "napA; nitrate reductase (cytochrome) [EC:1.9.6.1]":"#36f2ea",\
     "narB; ferredoxin-nitrate reductase [EC:1.7.7.2]":"#4e79d4",\
     "fdoG; formate dehydrogenase major subunit [EC:1.17.1.9]":"#4a7038",\
     "fdhF; formate dehydrogenase (acceptor) [EC:1.17.99.7]":"#507932",\
     "fdoG, fdhF, fdwA; formate dehydrogenase major subunit [EC:1.17.1.9]":"#307843",\
     "fdnG; formate dehydrogenase-N, alpha subunit [EC:1.17.5.3]":"#669150",\
     "phsA; thiosulfate reductase / polysulfide reductase chain A [EC:1.8.5.5]":"#64a0bc",\
     "phsA, psrA; thiosulfate reductase / polysulfide reductase chain A [EC:1.8.5.5]":"#64a0bc",\
     "aoxB; arsenite oxidase large subunit [EC:1.20.2.1 1.20.9.1]":"#b6ad11",\
     "cysJ; sulfite reductase (NADPH) flavoprotein alpha-component [EC:1.8.1.2]":"#d9f598",\
     "ttrA; tetrathionate reductase subunit A":"#8f92bd",\
     "soeA; sulfite dehydrogenase (quinone) subunit SoeA [EC:1.8.5.6]":"#bda38f",\
     "iscS, NFS1; cysteine desulfurase [EC:2.8.1.7]":"#111111",\
     "dmsD; putative dimethyl sulfoxide reductase chaperone":"#111111",\
     "ynfE; Tat-targeted selenate reductase subunit YnfE [EC:1.97.1.9]":"#222222",\
     "A": "#b4c2cf","B": "#b4c2cf","C": "#36f2ea","D": "#cc9761","E": "#b4c2cf","F": "#b4c2cf","G": "#b4c2cf","H": "#b4c2cf","I": "#b4c2cf",\
     "J": "#b4c2cf","K": "#b4c2cf","L": "#b4c2cf","M": "#b4c2cf","N": "#b4c2cf","O": "#b4c2cf","P": "#b4c2cf","Q": "#b4c2cf","R": "#b4c2cf",\
     "S": "#421682","T": "#b4c2cf","U": "#307843","W": "#b4c2cf","X": "#b4c2cf","Y": "#b4c2cf","Z": "#b4c2cf",\
     "-":"#ffffff",\
     "True":"#3c1611",\
     "False":"#f7e5e3"\
    }
anoWidthDict = \
    {"etnpall":"1.2",\
     "etnppool":"1.2",\
     "satl":"1.2",\
     "query":"2",\
     "marref":"1",\
     "mardb":"1",\
     "prokdb":"1"}
anoStyleDict = \
    {"etnpall":"normal",\
     "etnppool":"normal",\
     "satl":"normal",\
     "query":"bold",\
     "marref":"normal",\
     "mardb":"normal",\
     "prokdb":"normal"}


def writeColorStringFile(totalDF, irokiAnnotationOutputPathway, datasetLabel):
    filepath = irokiAnnotationOutputPathway + "_" + datasetLabel + "_colorString.txt"
    with open (filepath, 'w') as csOutfile:
        csOutfile.write("DATASET_COLORSTRIP\nSEPARATOR TAB\nDATASET_LABEL\t" + datasetLabel + "\n")
        csOutfile.write("COLOR\t#ff0000\nSTRIP_WIDTH\t25\nMARGIN\t0\nBORDER_WIDTH\t1\nBORDER_COLOR\t#000\nSHOW_INTERNAL\t0\nDATA\n")
        for index, row in totalDF.iterrows():
            taxa = row["taxa"]
            csOutfile.write(f"{taxa}\t{anoColorDict[row[datasetLabel]]}\t{row[datasetLabel]}\n")

def writeBinaryFile(totalDF, irokiAnnotationOutputPathway):
    filepath = irokiAnnotationOutputPathway + "_binary.txt"
    with open (filepath, 'w') as bOutfile:    
        bOutfile.write("DATASET_BINARY\nSEPARATOR TAB\nDATASET_LABEL\taaPos\nCOLOR\t#00ff00\n")
        bOutfile.write("FIELD_LABELS\tisW\tisY\ttype1motif\ttype2motif\tperiSecMotifMatchOver3\ttatfind\n")
        bOutfile.write("FIELD_COLORS\t#e6b800\t#806600\t#99ccff\t#333399\t#006699\t#2288bb\n")
        bOutfile.write("FIELD_SHAPES\t2\t2\t1\t1\t3\t3\nSHOW_INTERNAL\t1\nMARGIN\t10\nHEIGHT_FACTOR\t1\nSYMBOL_SPACING\t10\nDATA\n")
        for index, row in totalDF.iterrows():
            taxa = row["taxa"]
            isW = int(row["isW"])
            isY = int(row["isY"])
            type1motif = int(row["numc1"] >= 1)
            type2motif = int(row["numc2"] >= 1)
            periSecMatch = int(row["periSecMotifMatch"] >= 3)
            tatfind = int(row["hasTat"])
            bOutfile.write(f"{taxa}\t{isW}\t{isY}\t{type1motif}\t{type2motif}\t{periSecMatch}\t{tatfind}\n")

def writeDatasetStyleFile(totalDF, irokiAnnotationOutputPathway):
    filepath = irokiAnnotationOutputPathway + "_datasetStyle.txt"
    with open (filepath, 'w') as dsOutfile:    
        dsOutfile.write("DATASET_STYLE\nSEPARATOR TAB\nDATASET_LABEL\tdatabaseDatasetFormat\nCOLOR\t#000000\nDATA\n")
        for index, row in totalDF.iterrows():
            taxa = row["taxa"]
            database = row["database"]
            dsOutfile.write(f"{taxa}\tbranch\tnode\t{anoColorDict[database]}\t{anoWidthDict[database]}\tnormal\n")
            dsOutfile.write(f"{taxa}\tlabel\tnode\t{anoColorDict[database]}\t1\t{anoStyleDict[database]}\n")

def writeLeafWordingFile(totalDF, irokiAnnotationOutputPathway, leafWording):
    filepath = irokiAnnotationOutputPathway + "_" + leafWording + "_relabel.txt"
    with open (filepath, 'w') as lOutfile:    
        lOutfile.write("LABELS\nSEPARATOR TAB\nDATA\n")
        for index, row in totalDF.iterrows():
            taxa = row["taxa"]
            newLabel = row[leafWording]
            lOutfile.write(f"{taxa}\t{newLabel}\n")

def writeTextLabelFile(totalDF, irokiAnnotationOutputPathway, textAnnotation):
    filepath = irokiAnnotationOutputPathway + "_" +textAnnotation + "_textlabel.txt"
    with open (filepath, 'w') as tlOutfile:    
        tlOutfile.write(f"DATASET_TEXT\nSEPARATOR TAB\nDATASET_LABEL\t{textAnnotation}_textLabel\nCOLOR\t#0000ff\
                        \nMARGIN\t25\nSHOW_INTERNAL\t0\nALIGN_TO_TREE\t1\nDATA\n")
        for index, row in totalDF.iterrows():
            taxa = row["taxa"]
            textLabel = row[textAnnotation]
            tlOutfile.write(f"{taxa}\t{textLabel}\t1\t#909090\tnormal\t2\t0\n")


def main(args):
    """takes in sets of tab delimited annotations with "taxa" as a column and labled headers
    creates a set of annotation datasets for iToL"""
    ARGS_EXPECTED = 7
    name, wordSearchAnnotationPathway, keggGeneAnnotationPathway, keggTaxaAnnotationPathway, aaAnnotationPathway, tatfindAnnotationPathway, irokiAnnotationOutputPathway  = args[0:ARGS_EXPECTED]
    for num, otherArgs in enumerate(args[ARGS_EXPECTED:]):
        print("argument",num+ARGS_EXPECTED,": \"",otherArgs,"\" not used")
    #opens annotation files to pandas dataframes
    keggGeneDF = pandas.read_csv(keggGeneAnnotationPathway, sep='\t', usecols =["taxa","k_number_annotation", "k_score"])
    keggTaxaDF = pandas.read_csv(keggTaxaAnnotationPathway, sep='\t', usecols =["taxa","kegg_domain","kegg_class","kegg_species"])
    wordSearchDF = pandas.read_csv(wordSearchAnnotationPathway, sep='\t', usecols =["taxa","function", "description", "database"])
    aaDF = pandas.read_csv(aaAnnotationPathway, sep='\t', usecols =["taxa","isY","isW","uPos","numc1","numc2","periSecMotifMatch"])
    tatfindDF = pandas.read_csv(tatfindAnnotationPathway, sep='\t', usecols =["taxa","hasTat"])
    #merges things together; this is messy but oh well.
    df1 = pandas.merge(keggGeneDF, wordSearchDF, on=["taxa"])
    df2 = pandas.merge(df1, keggTaxaDF, on=["taxa"])
    df3 = pandas.merge(df2, aaDF, on=["taxa"])
    totalDF = pandas.merge(df3, tatfindDF, on=["taxa"])

    # makes a couple color string datasets
    
    for colorStringCat in ("k_number_annotation","function","database","uPos"):
        writeColorStringFile(totalDF, irokiAnnotationOutputPathway, colorStringCat)
    #writes a binary file dataset
    writeBinaryFile(totalDF, irokiAnnotationOutputPathway)
    #writes a dataset changing how the branches and leaves are stylized
    writeDatasetStyleFile(totalDF, irokiAnnotationOutputPathway)
    #writes a dataset for the original leaf text and one based on the description
    for leafWording in ("taxa", "description"):
        writeLeafWordingFile(totalDF, irokiAnnotationOutputPathway, leafWording)
    for textAnnotation in ("kegg_domain","kegg_class","kegg_species"):
        writeTextLabelFile(totalDF, irokiAnnotationOutputPathway, textAnnotation)
    
if __name__ == "__main__":
   main(sys.argv)