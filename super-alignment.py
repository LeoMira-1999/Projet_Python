from alignment import *
import glob
from itertools import chain
import reprlib
def cluster_reader():

    clusters_dict = {}
    species = []
    total_sp = []
    for file in glob.glob("aligned_clusters/*"):
        cleaned_file = file.split("aligned_clusters/")
        with open("aligned_clusters/"+cleaned_file[1], "r") as cluster_file:
            lines = cluster_file.readlines()
            lines = [w.replace('\n', '') for w in lines]
            temporary_dict = {}
            temporary_sp = []
            for line in lines:
                if line.startswith(">"):
                    raw_species = line.split("[")
                    species = raw_species[-1].split("]")
                    first_name = species[0].split(" ")

                    temporary_sp.append(first_name[0])
                    total_sp.append(first_name[0])

                if line.startswith(">") == False:

                    if temporary_sp[-1] not in temporary_dict:
                        temporary_dict[temporary_sp[-1]]= [line]
                    else:
                        temporary_dict[temporary_sp[-1]].append(line)

            for SP, seq in temporary_dict.items():
                if SP not in clusters_dict:
                    clusters_dict[SP] = [''.join(seq)]
                else:
                    clusters_dict[SP].append(''.join(seq))

    proteome_sp_real = []
    for raw_names in glob.glob("*.faa"):
        cleaned_names = raw_names.split("-protein.faa")
        names = cleaned_names[0].split("-")
        proteome_sp_real.append(names[0])

    proteome_sp = []
    for sp in total_sp:
        if sp not in proteome_sp:
            proteome_sp.append(sp)

    odd_sp =list(set(proteome_sp)-set(proteome_sp_real))


    for odd in odd_sp:
        del clusters_dict[odd]


    cleaned_clusters_dict = {}

    seqs = []
    for sp in proteome_sp_real:
        seqs.clear()
        seqs.append("".join(clusters_dict[sp]))
        cleaned_clusters_dict[sp]=seqs
    print(reprlib.repr(cleaned_clusters_dict))

cluster_reader()
