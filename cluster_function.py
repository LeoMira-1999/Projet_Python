#!/usr/bin/python3
#@ Authors: Mirandola Leonardo and Duplan Alexandre

import os #importing os to let us use bash
import itertools #import itertools


def RBH_comparator():
    """
    Arguments: None
    Returns: creates a dictionary that will contain as a key the filename and the value will be a list containing a list of each pairs
    Author: Mirandola Leonardo
    """
    #read each RBH filename and cconcatenate them
    os.system("ls reciprocal* > filename.txt")

    #open the file containing the names of the files
    f = open("filename.txt", "r")

    #read the lines
    filenames = f.readlines()

    #create a dictionnary
    dict_combination = {}

    #read each filename
    for filename in filenames:

        #remove each \n in the file
        filename = filename.replace("\n","")

        #create an empty list to store the prot pairs
        cleaned_prot_paired = []

        #open the filename to read what's inside
        with open(filename,"r") as file_content:

            #read the lines
            raw_prot_pairs = file_content.readlines()

            #loop inside the lines
            for prot_pairs in raw_prot_pairs:

                #remove \n
                prot = prot_pairs.replace("\n","")

                #create a list that is split by a space in the line
                prot = prot.split()

                #append the protein pair to a list of prot pairs
                cleaned_prot_paired.append(prot)
        #for each filename, set the key as the filename and the value will be a list in a list of prot pairs
        dict_combination[filename]=cleaned_prot_paired

    #remove the file that reads the reciprocal filenames
    os.system("rm filename.txt")

    #return the dictionary
    return dict_combination

def RBH_analysor(dict):
    """
    Arguments: takes the dictionnary created by RBH_comparator
    Returns: a list of all the non redundant clusters from the dictionnary storing the reciprocals
    Author: Mirandola Leonardo
    """

    #create an empty list to store the clusters
    cluster = []

    #read the key (filename), value(lists of prot pair inside a list) of the dictionnary function argument
    for RBH_filename1, prot_pairs1 in dict.items():

        #call each protein pair inside the list containing the pairs
        for unique_prot_pair1 in prot_pairs1:

            #create an empty temporary storage
            temporary = []

            #read each prot inside the prot pair list
            for unique_prot1 in unique_prot_pair1:

                #cycle a second time in the dictionnary in order to compare our first protein to others
                for RBH_filename2, prot_pairs2 in dict.items():

                    #if the filenames are the same there's no need to analyse
                    if RBH_filename1 != RBH_filename2:

                        #call each protein pair inside the list containing the pairs
                        for unique_prot_pair2 in prot_pairs2:

                            #if the protein is found in the other file
                            if unique_prot1 == unique_prot_pair2[0]:

                                #append the other protein to temporary
                                temporary.append(unique_prot_pair2[1])

                            #if the protein is found in the other file but here it's the second one that matches
                            elif unique_prot1 == unique_prot_pair2[1]:

                                #append the other protein to temporary
                                temporary.append(unique_prot_pair2[0])

                            #if there are no proteins that are found in the file
                            elif unique_prot1 not in temporary:

                                #append the couple to temporary
                                temporary.append(unique_prot1)

                            #remove duplicates
                            temporary = list(dict.fromkeys(temporary))

                            #sort the list
                            temporary = sorted(temporary)

            #once a prot has looped all the way, append the results to cluster
            cluster.append(temporary)

    #sort the list containing the clusters
    cluster_sorted = sorted(cluster)

    #remove duplicates of clusters
    cleaned_cluster = list(duplicates for duplicates,_ in itertools.groupby(cluster))

    #read each clusters
    for cluster1 in cleaned_cluster:

        #cycle a second time
        for cluster2 in cleaned_cluster:

            #if two clusters are not different
            if cluster1 != cluster2:

                #read each prot
                for prot in cluster1:

                    #if the prot from the first cluster is in the second cluster
                    if prot in cluster2:

                        #calculate the length of the first and second clusters
                        longueur_cluster1 = len(cluster1)
                        longueur_cluster2 = len(cluster2)

                        #if the length of cluster 1 is bigger than the length of cluster 2
                        if longueur_cluster1 > longueur_cluster2:

                            #add each cluster together
                            final = cluster1+cluster2

                            #remove duplicates
                            final_cluster = list(dict.fromkeys(final))

                            #remove each clusters
                            cluster1.clear()
                            cluster2.clear()

                            #append the result to the final cluster
                            cleaned_cluster.append(final_cluster)

                        #if the length of cluster 1 is smaller than the length of cluster 2
                        elif longueur_cluster1 < longueur_cluster2:

                            #add each cluster together
                            final = cluster1+cluster2

                            #remove duplicates
                            final_cluster = list(dict.fromkeys(final))

                            #remove each clusters
                            cluster1.clear()
                            cluster2.clear()

                            #append the result to the final cluster
                            cleaned_cluster.append(final_cluster)

                        #if cluster 1 and cluster 2 have the same size
                        else:

                            #add each cluster together
                            final = cluster1+cluster2

                            #remove duplicates
                            final_cluster = list(dict.fromkeys(final))

                            #remove each clusters
                            cluster1.clear()
                            cluster2.clear()

                            #append the result to the final cluster
                            cleaned_cluster.append(final_cluster)

    #remove empty clusters
    final_cluster = [empty for empty in cleaned_cluster if empty != []]

    #remove duplicates
    cleaned_final_cluster = list(duplicate for duplicate,_ in itertools.groupby(sorted(final_cluster)))

    #return the final list containing all the clusters clustered
    return cleaned_final_cluster

def cluster_species_finder(list):
    """
    Arguments: takes a list of protein clusters  (list of lists)
    Returns: list of organism clusters (list of lists)
    Author: Duplan Alexandre
    """

    #os command to get all file named " *-protein.faa " in the current directory and save them in "filename-protein.txt"
    os.system("ls *-protein.faa > filename-protein.txt")

    #read the file create
    file_prot = open("filename-protein.txt", "r")
    file_name= file_prot.readlines()

    #for each line, we just get the file name (we have a little problem because spaces can be add at the end of the line)
    #so we split the line and add the end ourself
    for line in file_name:
        line_2=line.split('.faa')
        filename=str(line_2[0])+'.faa'

        #loop which modifies the name of prot in clusters by the organism name:
        #we need to modified prot by prot in all clusters
        for cluster1 in list:
            for prot in cluster1:

                #open proteomes by proteomes (file.faa)
                with open(filename) as file:

                    #if a prot is in the proteomes
                    if prot in file.read():

                        organism_name=filename.split("-protein")#we just want the name before "-protein" => organism name

                        #for all proteins in the cluster, if it's in a proteomes file, change the protein name by the organism name
                        for index, item in enumerate(cluster1):
                            if prot in item:
                                cluster1[index] = organism_name[0]

    #remove the file create at the beginning
    os.system("rm filename-protein.txt")

    return list

def cluster_species_redundance_remover(cluster_AC, cluster_SP):
    """
    Arguments: takes the clusters with accession code and the same list of clusters but with species instead
    Returns: a non redundant list of clusters AC and SP but redundant species from the species cluster have been removed from both lists
    Author: Mirandola Leonardo
    """

    #set 2 empty lists
    nr_AC = []
    nr_SP = []

    #cycle in species cluster and with the index
    for counter_cluster, cluster in enumerate(cluster_SP):

        #set 2 temporary lists
        temporary_cluster_SP = []
        temporary_cluster_AC = []

        #cycle in each clusters and with the index
        for counter_SP, SP in enumerate(cluster):

            #if the species is not in the temporary list
            if SP not in temporary_cluster_SP:

                #append SP
                temporary_cluster_SP.append(SP)
                #append the AC from the cluster index and the general cluster index
                temporary_cluster_AC.append(cluster_AC[counter_cluster][counter_SP])

        #append the non redundant clusters
        nr_AC.append(temporary_cluster_AC)
        nr_SP.append(temporary_cluster_SP)

    #return the non redundant list of accession and species
    return nr_AC, nr_SP
