#!/usr/bin/python3
#@ Author: Mirandola Leonardo


import glob #import glob to read files
from Bio import Phylo #import phylo from bio to draw a phylogenetic tree
import os #import os to use bash

def gap_filer(proteomes):
    """
    Arguments: Takes the list of proteomes
    Results: fills the aligned clusters missing proteomes with gaps
    Author: Mirandola Leonardo
    """

    #create an empty list that stores the family names of the proteomes
    proteome_sp_real = []

    #retreive each proteome from proteomes
    for raw_names in proteomes:

        #remove -protein.faa
        cleaned_names = raw_names.split("-protein.faa")

        #take only the family name
        names = cleaned_names[0].split("-")

        #append the family name to the list of family names
        proteome_sp_real.append(names[0])

    #for each file in the directory aligned_clusters
    for file in glob.glob("aligned_clusters/*"):

        #remove aligned_clusters from the name file
        cleaned_file = file.split("aligned_clusters/")

        #open the edited filenames
        with open("aligned_clusters/"+cleaned_file[1], "r") as cluster_file:

            #read the lines
            lines = cluster_file.readlines()

            #replace each \n in the lines
            lines = [word.replace('\n', '') for word in lines]

            #create empty lists
            list_sp_file = []
            seq = []

            #initialize counter
            i = 0

            #begin loop until i has reached the length of lines
            while i < len(lines):

                #if the line startswith a ">"
                if lines[i].startswith(">"):

                    #split by the the first "["
                    raw_species = lines[i].split("[")

                    #split by the complementary "]"
                    species = raw_species[-1].split("]")

                    #retreive only the family name
                    first_name = species[0].split(" ")

                    #append the family name to the list of species
                    list_sp_file.append(first_name[0])

                #increment i by 1
                i+=1

                #create empty list
                one_seq= []

                #begin loop until i has reached the length of lines
                while i < len(lines):

                    #if the line startswith a ">"
                    if lines[i].startswith(">"):

                        #do nothing
                        break

                    #if the line doesn't start with a ">"
                    else:

                        #append the line of the sequence to the list of sequences of the family
                        one_seq.append(lines[i])

                    #increment i by 1
                    i+=1

                #concatenate a list of seqs in a list and append them to the sequence list
                seq.append("".join(one_seq))

            #create new file called processed_name_of_the_cluster_file
            with open("aligned_clusters/processed_"+cleaned_file[1], "a") as processed_cluster:

                #find out if there are missing species in the file original file
                missing_sp = list(set(proteome_sp_real)-set(list_sp_file))

                #if there are no missing species
                if len(missing_sp) == 0:

                    #retreive the species name and its sequence
                    for sp, sequence in zip(list_sp_file,seq):

                        #duplicate file and rename it as appropriate processed_cluster_name
                        processed_cluster.write(">"+sp+"\n"+sequence+"\n")

                #if there are 1 or more missing species
                if len(missing_sp) > 0:

                    #retreive the length of the first species sequence
                    seq_len = len(seq[0])

                    #retreive the species name and its sequence
                    for sp, sequence in zip(list_sp_file,seq):

                        #write down the species that are present in the actual file
                        processed_cluster.write(">"+sp+"\n"+sequence+"\n")

                    #for each missing species
                    for sp_missing in missing_sp:

                        #write down the species name and the sequence will be gaps equal to the length of the first species sequence
                        processed_cluster.write(">"+sp_missing+"\n"+seq_len*"-"+"\n")


def cluster_reader(proteomes):
    """
    Arguments: Takes the list of proteomes
    Result: creates the super_alignment file and phylogenetic tree
    Author: Mirandola Leonardo
    """

    #create empty dictionary
    clusters_dict = {}

    #create empty lists
    species = []
    total_sp = []

    #retreive from aligned clusters directory every file that startswith processed_
    for file in glob.glob("aligned_clusters/processed_*"):

        #remove aligned_clusters/ fromt the filename
        cleaned_file = file.split("aligned_clusters/")

        #open the filenames
        with open("aligned_clusters/"+cleaned_file[1], "r") as cluster_file:

            #read the lines
            lines = cluster_file.readlines()

            #replace each \n in the lines
            lines = [word.replace('\n', '') for word in lines]

            #create a temporary dictionnary
            temporary_dict = {}

            #create a temporary list
            temporary_sp = []

            #for every line
            for line in lines:

                #if the line begins with ">"
                if line.startswith(">"):

                    #split by ">"
                    first_name = line.split(">")

                    #retreive the family name of the species and append to temporary species list and total species
                    temporary_sp.append(first_name[1])
                    total_sp.append(first_name[1])

                #if the line doesn't start with ">"
                if line.startswith(">") == False:

                    #if the last item in the list of temporary species is not in the temporary_dict
                    if temporary_sp[-1] not in temporary_dict:

                        #set up the first key (species family name) and value (first line of the sequence)
                        temporary_dict[temporary_sp[-1]]= [line]

                    #if the species is already present in the dictionary
                    else:

                        #append the sequence to the already existing sequence
                        temporary_dict[temporary_sp[-1]].append(line)

            #for each species and sequences in the temporary dictionnary
            for SP, seq in temporary_dict.items():

                #if the species is not in the dictionnary
                if SP not in clusters_dict:

                    #set up the first key (species family name) and value (concatenation of each list of seqs inside the list)
                    clusters_dict[SP] = [''.join(seq)]

                #if the species is already present in the dictionary
                else:

                    #append the sequence to the already existing sequence
                    clusters_dict[SP].append(''.join(seq))

    #create empty list
    proteome_sp_real = []

    #for each protemones
    for raw_names in proteomes:

        #remove the -protein.faa
        cleaned_names = raw_names.split("-protein.faa")

        #split by each "-"
        names = cleaned_names[0].split("-")

        #store in the proteome_sp_real the family name of the proteome
        proteome_sp_real.append(names[0])

    #create empty list
    proteome_sp = []

    #retreive each species in the total species from the previous cluster analysis
    for sp in total_sp:
        if sp not in proteome_sp:
            proteome_sp.append(sp)

    #show the species that are not found from the proteome files
    odd_sp =list(set(proteome_sp)-set(proteome_sp_real))

    #for each of these odd species
    for odd in odd_sp:

        #remove the extra species that are not present in the species that were initialy selected by the user
        del clusters_dict[odd]

    #create an empty dictionary
    cleaned_clusters_dict = {}

    #for each species and sequence from the cluster dictionary
    for sp, seqs in clusters_dict.items():

        #join each list of sequence together for each species
        cleaned_clusters_dict[sp]="".join(seqs)

    #open a super_alignment file
    with open("final_super_alignment.afa" , "a") as file:

        #for each species and sequence from the cluster dictionary
        for SP, seq in cleaned_clusters_dict.items():

            #write down in a fasta format the species family name and its sequence
            file.write(">"+SP+"\n"+seq+"\n")

    #create a newick file by using RAxML with the super_alignment, which will be called final_tree, using fast bootstrap method (-f a), blosum62 matrix, seed 2, 10 bootstraps (1000 would be better but takes too long), -p 1 to combine the fast bootstrap and the seed 2
    os.system("raxmlHPC -s final_super_alignment.afa -n final_tree -f a -m PROTCATBLOSUM62 -x 2 -N 10 -p 1")

    #once the newick file is created retreive RAxML_bipartitions.name to draw the final tree and it will automatically let the user save or discard the tree drawn
    trees = Phylo.read('RAxML_bipartitions.final_tree', 'newick')
    Phylo.draw(trees)

    #create a directory that will have the name of each proteome family name analysed
    os.system("mkdir "+"-".join(proteome_sp_real))

    #move every file created for the analysis into the directory created above
    os.system("mv aligned_clusters clusters final_super_alignment.afa RAxML* reciprocal* "+"-".join(proteome_sp_real))
