from alignment import *
import glob
from itertools import chain
from Bio import Phylo

def gap_filer(proteomes):

    proteome_sp_real = []
    for raw_names in proteomes:
        cleaned_names = raw_names.split("-protein.faa")
        names = cleaned_names[0].split("-")
        proteome_sp_real.append(names[0])


    for file in glob.glob("aligned_clusters/*"):
        cleaned_file = file.split("aligned_clusters/")

        with open("aligned_clusters/"+cleaned_file[1], "r") as cluster_file:
            lines = cluster_file.readlines()
            lines = [w.replace('\n', '') for w in lines]

            list_sp_file = []
            seq = []
            i = 0
            while i < len(lines):
                if lines[i].startswith(">"):
                    raw_species = lines[i].split("[")
                    species = raw_species[-1].split("]")
                    first_name = species[0].split(" ")
                    list_sp_file.append(first_name[0])
                i+=1
                one_seq= []

                while i < len(lines):

                    if lines[i].startswith(">"):
                        break
                    else:
                        one_seq.append(lines[i])

                    i+=1
                seq.append("".join(one_seq))


            with open("aligned_clusters/processed_"+cleaned_file[1], "a") as processed_cluster:
                missing_sp = list(set(proteome_sp_real)-set(list_sp_file))

                if len(missing_sp) == 0:

                    for sp, sequence in zip(list_sp_file,seq):

                        processed_cluster.write(">"+sp+"\n"+sequence+"\n")
                if len(missing_sp) > 0:
                    seq_len = len(seq[0])
                    for sp, sequence in zip(list_sp_file,seq):

                        processed_cluster.write(">"+sp+"\n"+sequence+"\n")

                    for sp_missing in missing_sp:

                        processed_cluster.write(">"+sp_missing+"\n"+seq_len*"-"+"\n")





def cluster_reader(proteomes):

    clusters_dict = {}
    species = []
    total_sp = []
    for file in glob.glob("aligned_clusters/processed_*"):
        cleaned_file = file.split("aligned_clusters/")
        with open("aligned_clusters/"+cleaned_file[1], "r") as cluster_file:
            lines = cluster_file.readlines()
            lines = [w.replace('\n', '') for w in lines]
            temporary_dict = {}
            temporary_sp = []
            for line in lines:
                if line.startswith(">"):
                    first_name = line.split(">")

                    temporary_sp.append(first_name[1])
                    total_sp.append(first_name[1])

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
    for raw_names in proteomes:
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

    for sp, seqs in clusters_dict.items():
        cleaned_clusters_dict[sp]="".join(seqs)

    with open("final_super_alignment.afa" , "a") as file:

        """file.write(""+str(len(proteome_sp_real))+" "+str(len("".join(cleaned_clusters_dict[proteome_sp[0]])))+"\n")"""
        for SP, seq in cleaned_clusters_dict.items():

            file.write(">"+SP+"\n"+seq+"\n")

    os.system("raxmlHPC -s final_super_alignment.afa -n final_tree -f a -m PROTCATBLOSUM62 -x 2 -N 10 -p 1")

    trees = Phylo.read('RAxML_bipartitions.final_tree', 'newick')
    Phylo.draw(trees)
