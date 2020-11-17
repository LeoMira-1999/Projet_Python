#!/usr/bin/python3
import os #import os to use bash
import glob #import glob to read files


def RBH_DB_creator(list):
    """
    Arguments: takes in a list of the proteome files
    Result: create databases for each proteome
    """

    #for each filename in the list
    for filename in list:

        #remove the -protein.faa
        SP = filename.split("-protein.faa")

        #create a database that will have the name of the proteome
        os.system("makeblastdb -in "+filename+" -parse_seqids -blastdb_version 5 -dbtype prot -out "+SP[0]+"/"+SP[0]+"")


def cluster_alignment(non_redundant_AC_list,non_redundant_SP_list):
    """
    Arguments: list of non redundant AC and list of non redundant species from the cluster AC
    Returns:
    """

    #create a directory to store the clusters
    os.system("mkdir clusters")

    #initialise 2 counters
    j = 1
    i = 0

    #for each clusters in the non redundant AC
    for cluster in non_redundant_AC_list:

        #for each prot in the cluster
        for prot in cluster:

            #retreive the index of the cluster inside the non redundant AC list
            index1 = non_redundant_AC_list.index(cluster)

            #retreive the index of the cprot inside the cluster
            index2 = cluster.index(prot)

            #retreive the species by using the index of the non redundant AC onto the non redundant SP
            SP = non_redundant_SP_list[index1][index2]

            #retrieve sequence of the protein and store it in a file
            os.system("blastdbcmd -entry "+prot+" -db "+SP+"/"+SP+" -dbtype prot -out clusters/protein"+str(i)+".fa")

            #increment i by 1
            i += 1

        #concatenate each file created into 1 to recreate the original cluster
        os.system("cat clusters/protein* > clusters/raw_cluster"+str(j)+".fa")

        #increment j by 1
        j += 1

        #remove the already concatenated old protein files
        os.system("rm clusters/protein*")

    #create empty list
    raw_clusters =[]

    #call each file that start with raw in clusters directory
    for file in glob.glob("clusters/raw*"):

        #remove clusters/ from the filename
        filename = file.split("clusters/")

        #append each filename to raw_clusters
        raw_clusters.append(filename[1])

    #create firectory aligned clusters
    os.system("mkdir aligned_clusters")

    #cycle for each cluster name in raw cluster
    for cluster in raw_clusters:

        #remove the raw_ from cluster
        cleaned_cluster = cluster.split("raw_")

        #align each cluster with MAFFT and output result to aligned clusters directory
        os.system("mafft --anysymbol --auto --quiet clusters/"+cluster+" > aligned_clusters/"+cleaned_cluster[1]+"")



def RBH_DB_remover(list):
    """
    Arguments: takes in a list of the proteome files
    Returns: remove databases for each proteome
    """

    #for each file in list
    for file in list:

        #remove -protein.faa from file
        filename = file.split("-protein.faa")

        #remove the directory that is named like file
        os.system("rm -r "+filename[0]+"")
