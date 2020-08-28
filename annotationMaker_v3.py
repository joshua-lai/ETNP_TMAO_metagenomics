import sys
import re
from Bio import SeqIO

def readWriteAnnotation(inputFolderFilepath, outputFolderFilepath):
    """takes in a file path of a fastafile
    outputs a tab delimited file with some traits
    determined by the presence of a phrase within the description"""
    with open(inputFolderFilepath, 'r') as infile:
        with open(outputFolderFilepath, 'w') as outfile:

            generalAnnotationsList = ["taxa","description","function","database","species"]

            functionNameList = ["TMAO_reductase(torA)", "DMSO_reductase", "Nitrate_reductase", "Formate_reductase", "Pyrogallol_hydroxytransferase", \
                            "Something_BiotinSulfoxide", "Something_amidotransferase", "Something_Dehydrogenase", "Something_Molybdopterin"]
            functionTermsList = [["tora", "trimethyl","trimethylamine"], ["dora","dmso","dimethyl"], ["nitrate"], ["formate"], \
                                ["pyrogallol"], ["biotin","bisc"], ["amidotransferase"], ["dehydrogenase"], ["molybdopterin","molydopterin"]]
            """
            speciesList = ["photobacterium", "ruegeria","roseobacter","kitasatospora","puniceispirillum","alkaliphilus",
                            "scalindua","chloroflexus","Methylobacterium","Methylococcus","Syntrophobacter","Methylophaga","Desulfovibrio",
                            "Desulfocapsa","Desulfuromonas","Alkalilimnicola","Lutibacter","Oligotropha","Puniceispirillum","Marichromatium",
                            "Alcanivorax","Sulfitobacter","Saccharophagus","Thiobacillus"]
            """

            sourceNameList = ["query","etnppool","etnpall","prokdb","marref","mardb"]

            evalThresholdList = ["e-15","e-30","e-50","e-80","e-100"]
            evalThreshold = [1e-15, 1e-30, 1e-50, 1e-80, 1e-100]


            #writing the header
            for item in (generalAnnotationsList + functionNameList + sourceNameList + evalThresholdList):
                 outfile.write(item + "\t")
            outfile.write("\n")

            # per line
            for record in SeqIO.parse(infile, "fasta"):
                taxa = record.id
                description = record.description.replace(" ","_")
                
                #functionList and sourceList should be a lot of dashes and one name
                functionList = ['-'] * len(functionNameList)
                for num, functionName in enumerate(functionNameList):
                    #print(num,functionName)
                    #print(functionName, functionTermsList[num], recLow)
                    functionList[num] = searchXinY(functionName, functionTermsList[num], description)
                sourceList = ['-'] * len(sourceNameList)
                for num, sourceName in enumerate(sourceNameList):
                    sourceList[num] = searchXinY(sourceName, [sourceName], description)

                #for the more universal groups
                database = "unknown/other/recheck"
                for jj in range(len(sourceList)):
                    if sourceList[jj] != '-':
                        database = sourceNameList[jj]
                function = "unknown/other/recheck"
                species = "unknown/other/recheck"
                if database == "etnppool":
                    function = "from_ETNP_pool"
                    species = "from_ETNP_pool"
                elif database == "etnpall":
                    function = "from_ETNP_all"
                    species  = "from_ETNP_all"
                else:
                    #get species things
                    if database == "prokdb":
                        try:
                            species = str(re.search("(_)([A-Za-z]+_*proteobacter[a-z]+)(_)", description).group(2))
                        except:
                            try:
                                species = str(re.search('(_)([A-Z][a-z]{3,})(_)', description).group(2))
                            except:
                                try:
                                    species = str(re.search('([A-Za-z]+_[A-Za-z]+_[A-Za-z0-9]*_[A-Za-z0-9]*)(_prokdb_)', description).group(1))
                                except:
                                    species = "this_seems_to_have_weird_regex_for_prokdb"
                    elif database == "marref" or database == "mardb":
                        try:
                                species = str(re.search("(_)([A-Z][a-z]{3,})(_*[a-z]*)(_*?[A-Za-z0-9]*?)(_mmp_id)", description).group(2))
                        except:
                            try:
                                species = str(re.search('(\w+_)([A-Z][a-z]{3,})(\w+?)(_mmp_id)', description).group(2))
                            except:
                                species = "this_seems_to_have_weirder_regex_for_marref/mardb"
                    elif database == "query":
                        species = "i'm_too_lazy_to_code_this_for_the_queries_rn"
                    for ii in range(len(functionList)):
                        if functionList[ii] != '-':
                            function = functionNameList[ii]
                            break
                        



                #getting eval things
                evalList = ['-'] * len(evalThresholdList)
                try:
                    evalue = float(re.search('\d+\.\d+e-\d+$', description).group())
                except:
                    try:
                        evalue = float(re.search('0.0$', description).group())
                    except:
                        evalue = 99999
                for num,threshold in enumerate(evalThreshold):
                    if evalue < threshold:
                        evalList[num] = "includedAtThisEval"
                    else:
                        evalList[num] = '-'
                
                
                
                #writes to the file
                for entry in ([taxa,description,function,database,species] + functionList + sourceList + evalList):
                    outfile.write(entry + "\t")
                outfile.write("\n")

def searchXinY(name, terms, phrase):
    """takes in a name as a string, a set of matching terms as a list, and a phrase or description as a string,
    outputs the name if any of the terms are in the phrase, or a '-' otherwise"""
    output = '-'
    for term in terms:
        if term.lower() in phrase.lower():
            output = name
    return(output)



def main(args):
    """takes in a filepath for a fasta file,
    outputs a fasta file with the accession number, description, and a guessed function from a list"""
    # allows function to takes in input from the console
    ARGS_EXPECTED = 3
    name, inputFolderFilepath, outputFolderFilepath = args[0:ARGS_EXPECTED]
    for num, otherArgs in enumerate(args[ARGS_EXPECTED:]):
        print("argument",num+ARGS_EXPECTED,": \"",otherArgs,"\" not used")
    readWriteAnnotation(inputFolderFilepath, outputFolderFilepath)
    


if __name__ == "__main__":
    #testArgs = ["python annotationMaker_v2.py",  "C:/Users/joshu/Desktop/TMAO_project_files/alignments/clusterSet/cs_ePmRpDB_clusterSet_mafftAlignment.faa",  "C:/Users/joshu/Desktop/TMAO_project_files/annotations/cs_mafft_ePmRpDB_annotations_v2.txt"]
    #main(testArgs)
    main(sys.argv)
                
                