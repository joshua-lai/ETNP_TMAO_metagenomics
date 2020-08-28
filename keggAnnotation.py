import sys 

#basically removes the blank tab delimited bits to have dashes or something
with open (sys.argv[1], 'r') as keggList:
    with open (sys.argv[2], 'w') as keggOutput:
        counter = 0
        max = 10000
        moreToGo = True
        keggOutput.write("taxa\tk_number\tk_number_annotation\tk_score\n")
        while moreToGo:
            counter += 1
            # split the line by whitespace
            line = keggList.readline()
            #print(counter, line)
            line = line.split("\t")
            if len(line) > 1:
                try:
                    taxa, k_number, k_number_annotation, k_score =  line
                except:
                    taxa, k_number, k_number_annotation, k_score, k_2nd_best_annotation, k_score2 = line
                if k_number == "":
                    k_number = '-'
                if k_number_annotation == "":
                    k_number_annotation = '-'
                if k_score2 == "":
                    k_score2 = '-'
                if k_2nd_best_annotation == "":
                    k_2nd_best_annotation = '-'
                if k_score2 == "":
                    k_score2 = '-'
            else:
                moreToGo = False
            if counter > max:
                moreToGo = False
            keggOutput.write(f"{taxa}\t{k_number}\t{k_number_annotation}\t{k_score}\n")


