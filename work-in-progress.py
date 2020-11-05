from itertools import combinations
import os
import pandas as pd
import glob

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
                prot = prot.split()
                cleaned_prot_paired.append(prot)
        dict_combination[filename]=cleaned_prot_paired
    print(dict_combination)
    return dict_combination

def RBH_analysor(dict):

    cluster = []
    for RBH_filename1, prot_pairs1 in dict.items():
        for RBH_filename2, prot_pairs2 in dict.items():

            print(RBH_filename2)
            print(prot_pairs2)
            mini_cluster = []
            if RBH_filename1 != RBH_filename2 and prot_pairs1 != prot_pairs2:


                for unique_prot_pair2 in prot_pairs2:
                    print(unique_prot_pair2)

                    for unique_prot_pair1 in prot_pairs1:



                        if unique_prot_pair1[0] in unique_prot_pair2 or unique_prot_pair1[1] in unique_prot_pair2:
                            print(unique_prot_pair1[1],unique_prot_pair2)
                            print(unique_prot_pair1[0],unique_prot_pair2)
                            temporary=[]
                            for i in range(len(unique_prot_pair2)):
                                temporary.append(unique_prot_pair2[i])

                            temporary.append(unique_prot_pair1[1])
                            mini_cluster.append(temporary)



            print(mini_cluster)



RBH_analysor(RBH_comparator())
