def reciprocal_species_file_counter():
    os.system("ls reciprocal* > filename.txt")

    f = open("filename.txt", "r")

    lines = f.readlines()

    dict_combination = {}

    for line in lines:
        line = line.split("_")
        if line[1] == line[0]:
            line.remove(line[0])
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



def reciprocal_sequence_species_file_counter():
    os.system("ls sequence* > filename.txt")

    f = open("filename.txt", "r")

    lines = f.readlines()

    dict_combination = {}

    for line in lines:
        line = line.split("_")
        if line[1] == line[0]:
            line.remove(line[0])
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
        print(RBH_sequence_filename)

        RBH_species = RBH_sequence_filename.split("_")

        RBH_species.remove(RBH_species[0])
        RBH_species.remove(RBH_species[len(RBH_species)-1])
        print(species)
        print(RBH_species)
        for excluded_species in species:
            if excluded_species not in RBH_species:
                print(excluded_species)

                bidirectional_blast(RBH_sequence_filename, excluded_species)


retreive_RBH_species_sequence(reciprocal_species_file_counter())
comparing_RBH_to_diff_species(reciprocal_sequence_species_file_counter())
