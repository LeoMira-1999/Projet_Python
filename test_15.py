

import os #importing os to let us use bash
import itertools #import itertools

def cluster_species_finder():
    """
    Arguments: takes a list
    Returns:
    Author: Duplan Alexandre
    """

    os.system("ls *-protein.faa > filename-protein.txt")
    file_prot = open("filename-protein.txt", "r")
    file_name= file_prot.readlines()
    for line in file_name:
        line_2=line.split('.faa')
        filename=str(line_2[0])+'.faa'
#boucle qui modifie le nom des prot dans les clusters
        for cluster1 in list:
            for prot in cluster1:
                with open(filename) as file:
                    if prot in file.read():
                        organism_name=filename.split("-protein")
                        for index, item in enumerate(cluster1):
                            if prot in item:
                                cluster1[index] = organism_name[0]

    os.system("rm filename-protein.txt")

    return list
cluster_species_finder()
