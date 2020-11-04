from itertools import combinations
import os
import pandas as pd

def test_function(a,b):

    print(a,b)


def multi_RBH(*SP):
    """
    Arguments: can take N number of arguments
    Returns: all the possible non redondant combinations for each especies to have a RBH
    """

    #turns all the arguments into a list
    species = list(SP)

    #creates a non redondant pair of two combination
    comb = combinations(species,2)

    #loops in the list of combination
    for i in comb:

        #launches the RBH for that combination
        test_function(*i)

def RBH_comparator():

    os.system("ls RBH* > filename.txt")

    f = open("filename.txt", "r")

    filenames = f.readlines()

    dict_combination = {}
    for filename in filenames:
        filename = filename.replace("\n","")
        cleaned_prot_paired = []
        with open(filename,"r") as file_content:
            raw_prot_pairs = file_content.readlines()

            for prot_pairs in raw_prot_pairs:
                prot = prot_pairs.replace("\n","")
                prot = prot.split("\t")
                cleaned_prot_paired.append(prot)
        dict_combination[filename]=cleaned_prot_paired
    print(dict_combination)
    return dict_combination

def combinator(list1, list2):
    liste = list1 + list2
    comb = combinations(liste, 2)

    return comb

def RBH_analysor(dict):

    total_species = []
    RBH_filename_dict = {}

    for RBH_filename, total_prot_pairs in dict.items():
        print("The file name is: ",RBH_filename)
        RBH_species_filename = RBH_filename.split("_")

        print("Has these pairs: ",total_prot_pairs)
        RBH_filename_species = RBH_species_filename[1:len(RBH_species_filename)-1]

        for RBH_species_finder in RBH_filename_species:

            RBH_filename_dict[RBH_filename]=RBH_filename_species

            if RBH_species_finder not in total_species:
                total_species.append(RBH_species_finder)

    for RBH_filename, RBH_filename_species in RBH_filename_dict.items():
        print("first species is: ",RBH_filename_species[0])
        print("second species is: ",RBH_filename_species[1])
        to_do_combination = []
        to_combine_with = RBH_filename_species
        for species in total_species:

            if RBH_filename_species[0] != species and RBH_filename_species[1] != species:
                print("found new link:",species)
                to_do_combination.append(species)

        print("to do combinations: ",to_do_combination)
        print("to combine with: ", to_combine_with)
        possible_combinations = combinator(to_do_combination, to_combine_with)
        for true_combinations in possible_combinations:
            print([*true_combinations])
            if true_combinations == RBH_filename_species:
                print("this is possible:", true_combinations)


    print("Total species: ", total_species)
    print("These are all the files: ", RBH_filename_dict)



SP4 = "bact1-prot.faa"
SP5 = "bact2-prot.faa"
SP6 = "E-coli-prot.faa"

"""multi_RBH(SP1,SP2,SP3)"""
RBH_analysor(RBH_comparator())
