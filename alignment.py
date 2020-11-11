from main import *
from cluster_function import *
import os
import glob

def RBH_DB_creator(list):

    for filename in list:
        print(filename)
        SP = filename.split("-protein.faa")

        os.system("makeblastdb -in "+filename+" -parse_seqids -blastdb_version 5 -dbtype prot -out "+SP[0]+"/"+SP[0]+"")


def cluster_alignment(non_redundant_AC_list,non_redundant_SP_list):

    os.system("mkdir clusters")


    j = 1
    i = 0
    for cluster in non_redundant_AC_list:

        for prot in cluster:

            index1 = non_redundant_AC_list.index(cluster)
            index2 = cluster.index(prot)

            SP = non_redundant_SP_list[index1][index2]

            os.system("blastdbcmd -entry "+prot+" -db "+SP+"/"+SP+" -dbtype prot -out clusters/cluster"+str(i)+".fa")

            i += 1

        os.system("cat clusters/cluster* > clusters/raw_cluster"+str(j)+".fa")

        j += 1

        os.system("rm clusters/cluster*")

    raw_clusters =[]
    for file in glob.glob("clusters/raw*"):
        filename = file.split("clusters/")
        raw_clusters.append(filename[1])

    os.system("mkdir aligned_clusters")
    for cluster in raw_clusters:
        cleaned_cluster = cluster.split("raw_")
        os.system("mafft --auto --quiet clusters/"+cluster+" > aligned_clusters/"+cleaned_cluster[1]+"")



def RBH_DB_remover(list):

    for file in list:
        filename = file.split("-protein.faa")
        os.system("rm -r "+filename[0]+"")
