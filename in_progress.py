from itertools import combinations
import os

def test_function(a, b):

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


def reciprocal_species_file_counter():
    os.system("ls reciprocal* > filename.txt")

    f = open("filename.txt", "r")

    lines = f.readlines()

    dict_combination = {}

    for line in lines:
        line = line.split("_")
        key = '_'.join(line[1:len(line)-1])
        full_key= 'reciprocal-hits_'+key+'_.ids'
        dict_combination[full_key]=line[1:len(line)-1]

    os.system("rm filename.txt")
    return dict_combination

def retreive_RBH_species_sequence(dict):

    for RBH_filename, species_filename in dict.items():

        first_species = ''.join(species_filename[0])
        first_species_name = first_species.split(".")[0]


        os.system("makeblastdb -in "+first_species+" -parse_seqids -blastdb_version 5 -dbtype prot -out "+first_species_name+"/"+first_species_name+"")


        RBH_filename_sequences = "sequence-"+RBH_filename+".fa"
        os.system("blastdbcmd -entry_batch "+RBH_filename+" -db "+first_species_name+"/"+first_species_name+" -dbtype prot -out "+RBH_filename_sequences+"")
        os.system("rm -r "+first_species_name+"")


def reciprocal_sequence_species_file_counter():
    os.system("ls sequence* > filename.txt")

    f = open("filename.txt", "r")

    lines = f.readlines()

    dict_combination = {}

    for line in lines:
        line = line.split("_")
        key = '_'.join(line[1:len(line)-1])
        full_key= 'sequence-reciprocal-hits_'+key+'_.ids.fa'
        dict_combination[full_key]=line[1:len(line)-1]

    os.system("rm filename.txt")
    return dict_combination


def comparing_RBH_to_diff_species(dict):

    species = []

    for species_filename in dict.values():

        for called_name in species_filename:

            if called_name not in species:
                species.append(called_name)

    for RBH_sequence_filename, species_filename in dict.items():
         RBH_species = []

         RBH_sequence_filename.split("_").append(RBH_species)
         del RBH_species[1:len(RBH_species)]
         print(RBH_species)


comparing_RBH_to_diff_species(reciprocal_sequence_species_file_counter())


"""#multi_RBH("a","b","c")
retreive_RBH_species_sequence(reciprocal_species_file_counter())"""
