"""
15-110 Hw6 - Protein Sequencing Project
Name: M.Bharadwaj
AndrewID: 2023501046
"""

import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):                                           # This function is used to retrieve the NIH data.
    dna = ""
    f = open(filename,"r")                                        # Reading the text file and adding each string.
    reader = f.read()
    reader = reader.strip().split()
    for i in reader:
        dna += i
    f.close()
    return dna


'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):                                                 # this function is used to convert DNA to RNA by reading the DNA string.
    codons = []
    for i in range(startIndex,len(dna),3):
        codon = dna[i:i+3]
        if codon in ["TAA", "TGA", "TAG"]:                                     # ["TAA","TAG","TGA"] are the codon signals for start and end
            codons.append(codon.replace("T","U"))
            break
        codons.append(codon.replace("T","U"))
    return codons


'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):                     # This function is used to map each codon in codons list to amino acids to turn RNA into protein. This is done by reading the json file.
    import json
    codon = {}
    f = open(filename,"r")
    reader = json.load(f) 
    for aminoAcid, codons in reader.items():
        for i in codons:
            j = i.replace("T","U")
            codon[j] = aminoAcid
    return codon


'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):                     # This function is used to turn RNA into protein by identifing each RNA sequence and its associated amino acid to chain and that chain becomes our protein.
    protein = []
    for i in range(len(codons)):
        if i == 0 and codons[i] == "AUG":
            protein.append("Start")
        elif codonD[codons[i]] == "Stop":
            protein.append(codonD[codons[i]])
            break
        else:
            protein.append(codonD[codons[i]])

    return protein


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):              # This function is used to synthesize proteins from data file.
    dna = readFile(dnaFilename)                       # Reading the DNA from file.
    CodonDictionary = makeCodonDictionary(codonFilename)       # Producing the codons dictionary.
    proteins = []
    count = 0
    i = 0
    while i < len(dna):
        if dna[i:i+3] == "ATG":
            rna = dnaToRna(dna,i)
            prt = generateProtein(rna,CodonDictionary)
            proteins.append(prt)
            i += 3 * len(rna)
        else:
            i += 1
            count += 1
    print("Total Bases:", len(dna))
    print("Unused bases:", count)
    print("Total synthesized:", len(proteins))

    return proteins


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):          # This function is used to determine the common proteins that occur in both genes for checking the similarities in both the genes.
    uniqueProteins = []                                  # This list contains unique proteins that occur in both protein lists.
    for i in proteinList1:
        for j in proteinList2:
            if i == j and i not in uniqueProteins:
                # print(i,j)
                uniqueProteins.append(i)

    return uniqueProteins


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):                    # This function is used for comparing the amino acids produced for two genes by combining them.
    aminoAcids = []                                  # This list contains amino acids that occur in protein list.
    for i in proteinList:
        aminoAcids.extend(i)
    return aminoAcids


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):                     #This function is used to determine the frequences of each amino acid in each gene
    aminoAciddict = {}                               # Amino acid dictionary for frequency
    for aminoAcid in aaList:
        if aminoAcid in aminoAciddict:
            aminoAciddict[aminoAcid] += 1
        else:
            aminoAciddict[aminoAcid] = 1
    return aminoAciddict


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):                            # this function is used to compare amino acids frequencies between two genes of different lengths
    combinedProteinList1 = combineProteins(proteinList1)     
    combinedProteinList2 = combineProteins(proteinList2)
    amino_acid_dict1 = aminoAcidDictionary(combinedProteinList1)
    amino_acid_dict2 = aminoAcidDictionary(combinedProteinList2)
    elementslist = []                           # this list contains three elements first is amino acid and second and third is frequency of that amino acid in protein list1 and protein list2 and that frequency must be greater than given cutoff
    for i in amino_acid_dict2:
        if i not in amino_acid_dict1:
            amino_acid_dict1[i] = 0
    for j in amino_acid_dict1:
        if j not in amino_acid_dict2:
            amino_acid_dict2[j] = 0
    for k in amino_acid_dict1:
        if (k != "Start") and (k != "Stop") and (abs((amino_acid_dict1[k]/len(combinedProteinList1))-(amino_acid_dict2[k]/len(combinedProteinList2)))>cutoff):
            elementslist.append([k,amino_acid_dict1[k]/len(combinedProteinList1),amino_acid_dict2[k]/len(combinedProteinList2)])
    return elementslist


'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):                    # This function takes commons proteins for two genes and big differences for two genes and displays the result      
    for i in commonalities:
        print(" ".join(i))
        print()
    for j in differences:
        print(f"Amino acids are : {j[0]}")
        print(f"Percentage in set1 : {j[1]:.3f}")
        print(f"Percentage in set2 : {j[2]:.3f}")
        print()
    return


def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):                       # this takes the amino acids from both the protein list and returns the amino acid labels
    combinedProteins = combineProteins(proteinList1 + proteinList2)
    AminoAcidDictionary = aminoAcidDictionary(combinedProteins)
    AminoAcidLabels = sorted(AminoAcidDictionary.keys()) 
    return AminoAcidLabels


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):                           # this function takes labels list and protein list and returns frequency list which contains frequencies of each gene in protein list
    combinedProteinlist = combineProteins(proteinList)
    dictonary = aminoAcidDictionary(combinedProteinlist)
    frequency = []
    for i in labels:
        if i in combinedProteinlist:
            frequency.append(dictonary[i]/len(combinedProteinlist))
        else:
            frequency.append(0.0)
    return frequency


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):                           # this function is used to create and to plot the graph of frequences and labels in sequences
    import matplotlib.pyplot as plt

    a = 0.35
    plt.bar(xLabels, freqList1, width = -a,align = "edge", label = label1, color = edgeList)
    plt.bar(xLabels, freqList2, width = a, align = "edge", label = label2, color = edgeList)
    plt.legend()
    plt.show()

    return


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):                  # this function will return the 1D list in which it contains each element of black color if that element is in biggestDiffs else it other elements in white color
    colorList = []
    # for items in biggestDiffs
    aminoAcid = [items[0] for items in biggestDiffs]
    for i in labels:
        if i in aminoAcid:                               # if it is in amino acid list black color is appended
            colorList.append("black")
        else:
            colorList.append("white")                    # else white color is appended
    return colorList


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():                                                   #This is the main function to generate graph and the flow of prcess takes place to generate the analysis 
    Human_protein = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")       #Generated the human protein sequences
    Elephant_protein = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json") #Generated the elephant protein sequences
    commons = commonProteins(Human_protein,Elephant_protein)                                # Took the common sequences from human and elephant protein sequences
    Differencs = findAminoAcidDifferences(Human_protein,Elephant_protein,0.005)             # Took the differences from human and elephant protein sequences
    displayTextResults(commons,Differencs)                                                  # this function is used to display the commons and differences in protein sequences
    xlables = makeAminoAcidLabels(Human_protein,Elephant_protein)                           # amino acid labels for protein sequences
    First_frequency = setupChartData(xlables,Human_protein)                                 # frequences in labels and huamn protein sequences
    Second_frequency = setupChartData(xlables,Elephant_protein)                             # frequences in labels and elephant protein sequences
    Edge_list = makeEdgeList(xlables,Differencs)                                            # to plot the graph
    createChart(xlables,First_frequency,"Human", Second_frequency,"Elephant", Edge_list)    # to create chart for the graph
    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    test.week1Tests()
    test.testReadFile()
    test.testDnaToRna()
    test.testMakeCodonDictionary()
    test.testGenerateProtein()
    test.testSynthesizeProteins()
    print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    runWeek1()

    ## Uncomment these for Week 2 ##
    
    print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    test.week2Tests()
    test.testCommonProteins()
    test.testCombineProteins()
    test.testAminoAcidDictionary()
    test.testFindAminoAcidDifferences()
    print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    runWeek2()
    
    ## Uncomment these for Week 3 ##
    
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.testMakeAminoAcidLabels()
    test.testSetupChartData()
    test.testCreateChart()
    test.testMakeEdgeList()
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
    
