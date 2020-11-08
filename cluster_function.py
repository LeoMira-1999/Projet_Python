import os
import itertools


def RBH_comparator():

    os.system("ls reciprocal* > filename.txt")

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
    os.system("rm filename.txt")
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

    cluster.sort()

    cleaned_cluster = list(k for k,_ in itertools.groupby(cluster))

    for cluster1 in cleaned_cluster:
        for cluster2 in cleaned_cluster:
            if cluster1 != cluster2:
                for prot in cluster1:
                    if prot in cluster2:
                        longueur_cluster1 = len(cluster1)
                        longueur_cluster2 = len(cluster2)
                        if longueur_cluster1 > longueur_cluster2:
                            final = cluster1+cluster2
                            final_cluster = list(dict.fromkeys(final))

                            cluster1.clear()
                            cluster2.clear()
                            cleaned_cluster.append(final_cluster)

                        elif longueur_cluster1 < longueur_cluster2:
                            final = cluster1+cluster2
                            final_cluster = list(dict.fromkeys(final))

                            cluster1.clear()
                            cluster2.clear()
                            cleaned_cluster.append(final_cluster)

                        else:
                            final = cluster1+cluster2
                            final_cluster = list(dict.fromkeys(final))

                            cluster1.clear()
                            cluster2.clear()
                            cleaned_cluster.append(final_cluster)


    final_cluster = [x for x in cleaned_cluster if x != []]

    return final_cluster


print(RBH_comparator())
print(RBH_analysor(RBH_comparator()))