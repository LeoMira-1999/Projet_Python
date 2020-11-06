from itertools import combinations
import os
import pandas as pd
import glob
from operator import itemgetter

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

    return dict_combination

def RBH_analysor(dict):

    cluster = []
    for RBH_filename1, prot_pairs1 in dict.items():

        for unique_prot_pair1 in prot_pairs1:
            temporary = []

            for unique_prot1 in unique_prot_pair1:

                for RBH_filename2, prot_pairs2 in dict.items():


                    if RBH_filename1 != RBH_filename2:

                        for unique_prot_pair2 in prot_pairs2:


                            if unique_prot1 == unique_prot_pair2[0]:
                                temporary.append(unique_prot_pair2[1])

                            elif unique_prot1 == unique_prot_pair2[1]:
                                temporary.append(unique_prot_pair2[0])

                            elif unique_prot1 not in temporary:
                                temporary.append(unique_prot1)

                            temporary = list(dict.fromkeys(temporary))

                            temporary = sorted(temporary)

            cluster.append(temporary)
    print(cluster)

test = [["a","b"],["a","b"],["a","b","c"],["d","c"]]

test1 = list(set(zip(*test)))
print(test1)

RBH_analysor(RBH_comparator())
