from main import *
from cluster_function import *
import os

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

def RBH_DB_remover(list):

    for file in list:
        os.system("rm -r "+file+"")


proteomes = proteome_file_finder()

RBH_DB_creator(proteomes)

cluster_AC = RBH_analysor(RBH_comparator())

cluster_SP_nr = cluster_species_redundance_remover(cluster_AC)

"""print(cluster_SP_nr)
print(len(cluster_SP_nr))"""

cluster_alignment(RBH_analysor(RBH_comparator()),cluster_SP_nr)

RBH_DB_remover(proteomes)
