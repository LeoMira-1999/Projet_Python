#!/usr/bin/python3
#@ Author: Mirandola Leonardo and Duplan Alexandre


import os #importing os to let us use bash
import multiprocessing #importing multiprocessing will let us count the amount of cpu cores for an optimised use and faster processing of blastp
import pandas as pd #importing pandas to read our blasted files and retreive the information needed
import math #import math to use math.floor in mean_prot_length_evalue function
from itertools import combinations # imorting combinations from itertools to generate a combination of non-redundant pairs for the multi_RBH function
import glob #import glob to read files


def mean_prot_length_evalue(SP):
    """
    Arguments: a species file name
    Returns: calculated evalue it should use in a later blast program
    Author: Mirandola Leonardo
    """
    #use stats from bash seqkit program to get information about our species and store it in a file
    os.system("./seqkit stats "+SP+" > SP_mean_prot_length.txt")

    #open the stored file
    f= open("SP_mean_prot_length.txt","r")

    #read the opened file
    lines = f.readlines()

    #get second line
    line = lines[1]

    #replace each \n with nothing
    line = line.replace("\n","")

    #create a list of our second line where each space is a separator
    line = line.split()

    #use bash to remove the previously created file
    os.system("rm SP_mean_prot_length.txt")

    #remove the anglosaxon thousand coma
    mean_aa = "".join(line[6].split(","))

    #calculate our evalue from the mean protein length located on the 7th separator and divide it by 100
    #multiply this ration by 20 in order to have a dynamic evalue where every 100 amino acids we add 20 to the evalue
    evalue = (float(mean_aa)/100)*20

    #remove the decimals and transform it into a string
    evalue = int(math.floor(evalue))

    #return the evalue
    return str(evalue)

def bidirectional_blast(SP1, SP2):
    """
    Arguments: species 1 file name and species 2 file name
    Returns: a document where you hava the AC of each reciprocal hits
    Authors: Mirandola Leonardo and Duplan Alexandre
    """
    # retreive species 1 name
    species_name_1 = SP1.split("-protein.faa")

    # retreive species 1 name
    species_name_2 = SP2.split("-protein.faa")

    #use bash seqkit program to separate species 1 into 100 subfiles, in order to handle proteomes and minimise RAM usage
    #these files will be stored in a directory called subset_SP1
    os.system("./seqkit split "+SP1+" --by-part 100 --out-dir subset_SP1/")

    #create a file that has the names of each files from subset_SP1
    os.system("ls subset_SP1/ > files_SP1.txt")

    #create a directory to store first blasted results
    os.system("mkdir result_blast_1vs2")

    #open the file containing the names of each of the 100 files from subset_SP1
    f = open("files_SP1.txt","r")

    #read the lines
    lines = f.readlines()

    #calculate the evalue for the first species
    evalue = mean_prot_length_evalue(SP1)

    #setting a counter to 0
    i = 0

    #for loop to access each of the 100 files one by one
    for line in lines:

        #remove the \n form it file
        line = line.replace("\n","")

        #launching our blastp with each file name, using our database of our second species, format with tabular and comment lines, count the number of cpu cores the computer has to maximise efficiency
        #each of the results will be stored in the directory result_blast_1vs2 where the will be labled with the same number of the query file
        blastp_SP1_vs_SP2 = "blastp -query subset_SP1/"+str(line)+" -db "+species_name_2[0]+"/"+species_name_2[0]+" -outfmt 7 -evalue 1e-"+evalue+" -num_threads "+str(multiprocessing.cpu_count())+" > result_blast_1vs2/blast_raw_1vs2_0"+str(i+1)+".fa"

        #launch command for ourblastp
        os.system(blastp_SP1_vs_SP2)

        #increment counter by 1
        i+=1

    #concatenate each blasted results in result_blast_1vs2 as one 1
    os.system("cat result_blast_1vs2/* > blast_raw_1vs2.fna")

    #remove temporary files and folders
    os.system("rm files_SP1.txt")
    os.system("rm -r result_blast_1vs2 subset_SP1")

    #storing in a file the header to be used by our processed blast after having selected only the first hits (best hits) of each protein
    processing_1 = "echo \"query acc.ver\tsubject acc.ver\t% identity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. start\ts. end\tevalue\tbit score\" > temp_file.txt | cat blast_raw_1vs2.fna |awk '/hits found/{getline;print}' | grep -v \"#\" > blast_best_hits_1vs2.fa"
    os.system(processing_1)

    #if the size of the file is 0 then create an empty file with the appropriate name
    if os.stat("blast_best_hits_1vs2.fa").st_size == 0:
        return os.system("touch reciprocal-hits_"+SP1+"_"+SP2+"_.ids")

    #concatenate the header and the best hits
    processing_2 = "cat temp_file.txt blast_best_hits_1vs2.fa > header_blast_1vs2.fa"
    os.system(processing_2)

    #read the file with pandas
    df = pd.read_csv("header_blast_1vs2.fa", sep='\t')

    #retreive the accession of the subject
    column_ids_seq2 = df['subject acc.ver']

    #remove the index
    column_ids_seq2.reset_index(drop=True, inplace=True)

    #create a file that stores the accession of the subject best hits
    with open('SP2.ids', 'x') as f:
        #remove a index a second time (required)
        column_ids_seq2.to_csv(f, sep='\t',header=False, index=False)

    #retreive accession sequences from the SP2.ids from database ouf our second species to later be blasted
    SP2_blast_hits_seq_query = "blastdbcmd -entry_batch SP2.ids -db "+species_name_2[0]+"/"+species_name_2[0]+" -dbtype prot -out SP2_seq_best_hits_blast_1vs2.fa"

    #launch command
    os.system(SP2_blast_hits_seq_query)

    #separe our file in 100 subfiles in a directory called subset_SP2_Blast_reciprocal
    os.system("./seqkit split SP2_seq_best_hits_blast_1vs2.fa --by-part 100 --out-dir subset_SP2_Blast_reciprocal/")

    #create a file that has the names of each file from subset_SP2_Blast_reciprocal
    os.system("ls subset_SP2_Blast_reciprocal/ > files_SP2_blast_reciprocal.txt")

    #create a directory to store second blast
    os.system("mkdir result_blast_2vs1_reciprocal")

    #open the file containing the name of all the files from subset_SP2_Blast_reciprocal
    f = open("files_SP2_blast_reciprocal.txt","r")

    #read file
    lines = f.readlines()

    #calculate evalue for the mean length of our best hits file
    evalue = mean_prot_length_evalue("SP2_seq_best_hits_blast_1vs2.fa")

    #start counter
    i = 0

    #for loop to access each of the 100 files one by one
    for line in lines:

        #remove the \n form it file
        line = line.replace("\n","")

        #launching our blastp with each file name, using our database of our first species, format with tabular and comment lines, count the number of cpu cores the computer has to maximise efficiency
        #each of the results will be stored in the directory result_blast_2vs1 where the will be labled with the same number of the query file
        blastp_SP1_vs_SP2 = "blastp -query subset_SP2_Blast_reciprocal/"+str(line)+" -db "+species_name_1[0]+"/"+species_name_1[0]+" -outfmt 7 -evalue 1e-"+evalue+" -num_threads "+str(multiprocessing.cpu_count())+" > result_blast_2vs1_reciprocal/blast_raw_2vs1_reciprocal_0"+str(i+1)+".fa"

        #launch command
        os.system(blastp_SP1_vs_SP2)

        #increment counter
        i+=1

    #concatenate each blasted results in result_blast_1vs2 as one 1
    os.system("cat result_blast_2vs1_reciprocal/* > blast_raw_2vs1_reciprocal.fna")

    #storing in a file the header to be used by our processed blast after having selected only the first hits (best hits) of each protein
    processing_1 = "echo \"query acc.ver\tsubject acc.ver\t% identity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. start\ts. end\tevalue\tbit score\" > temp_file.txt | cat blast_raw_2vs1_reciprocal.fna |awk '/hits found/{getline;print}' | grep -v \"#\" > blast_reciprocal_2vs1.fa"
    os.system(processing_1)

    #concatenate the header and the best hits
    processing_2 = "cat temp_file.txt blast_reciprocal_2vs1.fa > header_blast_2vs1.fa"
    os.system(processing_2)

    #read the file with pandas
    df = pd.read_csv("header_blast_2vs1.fa", sep='\t')

    #retreive the accession of the subject
    column_ids_seq2 = df[['query acc.ver','subject acc.ver']]

    #remove the index
    column_ids_seq2.reset_index(drop=True, inplace=True)

    #create a file that stores the accession of the subject best hits
    with open('reciprocal_hits.ids', 'x') as f:
        #remove a index a second time (required)
        column_ids_seq2.to_csv(f, sep='\t',header=False, index=False)

    #rename file with appropriate to the species analysis and remove duplicates
    os.system("cat reciprocal_hits.ids > reciprocal-hits_"+SP1+"_"+SP2+"_.ids")


    #remove temporary files and folders
    os.system("rm blast_raw_1vs2.fna blast_raw_2vs1_reciprocal.fna reciprocal_hits.ids header_blast_2vs1.fa header_blast_1vs2.fa SP2.ids files_SP2_blast_reciprocal.txt  blast_best_hits_1vs2.fa temp_file.txt SP2_seq_best_hits_blast_1vs2.fa blast_reciprocal_2vs1.fa ")
    os.system("rm -r subset_SP2_Blast_reciprocal result_blast_2vs1_reciprocal")

def multi_RBH(*SP):
    """
    Arguments: can take N number of arguments
    Returns: all the possible non redondant combinations for each especies to have a RBH
    Author: Mirandola Leonardo
    """

    #turns all the arguments into a list
    species = list(SP)

    #creates a non redondant pair of two combination
    comb = combinations(species,2)

    #loops in the list of combination
    for i in comb:

        #launches the RBH for that combination
        bidirectional_blast(*i)

def proteome_file_finder():
    """
    Arguments: none
    Returns: a list of all the proteome files in the directory of this file
    Author: Mirandola Leonardo
    """
    #create empty list
    faa_files =[]
    #cycle for each file that ends with .faa
    for file in glob.glob("*.faa"):
        #append in the list the file
        faa_files.append(file)
    #return the list once finished
    return faa_files
