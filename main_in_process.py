import os
import multiprocessing
import pandas as pd

SP1 = "arc1.fna"
SP2 = "arc2.fna"

os.system("./seqkit split "+SP1+" --by-part 100 --out-dir subset_SP1/")
os.system("./seqkit split "+SP2+" --by-part 100 --out-dir subset_SP2/")
os.system("makeblastdb -in "+SP2+" -parse_seqids -blastdb_version 5 -dbtype nucl -out SP2/SP2")
os.system("makeblastdb -in "+SP1+" -parse_seqids -blastdb_version 5 -dbtype nucl -out SP1/SP1")
os.system("ls subset_SP1/ > files_SP1.txt")
os.system("ls subset_SP2/ > files_SP2.txt")
os.system("mkdir result_blast_1vs2")

f = open("files_SP1.txt","r")

lines = f.readlines()

i = 0

for line in lines:
    line = line.replace("\n","")
    blastn_SP1_vs_SP2 = "blastn -query subset_SP1/"+str(line)+" -db SP2/SP2 -outfmt 7 -subject_besthit -num_threads "+str(multiprocessing.cpu_count())+" > result_blast_1vs2/blast_raw_1vs2_0"+str(i+1)+".fa"
    os.system(blastn_SP1_vs_SP2)
    i+=1

os.system("cat result_blast_1vs2/* > blast_raw_1vs2.fna")
os.system("rm files_SP1.txt files_SP2.txt")
os.system("rm -r result_blast_1vs2 subset_SP2 subset_SP1")


processing_1 = "echo \"query acc.ver\tsubject acc.ver\t% identity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. start\ts. end\tevalue\tbit score\" > temp_file.txt | cat blast_raw_1vs2.fna |awk '/hits found/{getline;print}' | grep -v \"#\" > blast_best_hits_1vs2.fa"
os.system(processing_1)

processing_2 = "cat temp_file.txt blast_best_hits_1vs2.fa > header_blast_1vs2.fa"
os.system(processing_2)

df = pd.read_csv("header_blast_1vs2.fa", sep='\t')

print(df)
column_ids_seq2 = df['subject acc.ver']
column_ids_seq2.reset_index(drop=True, inplace=True)

print(column_ids_seq2)

if os.path.isfile("SP2.ids"):
    pass
else:
    with open('SP2.ids', 'x') as f:
        column_ids_seq2.to_csv(f, sep='\t',header=False, index=False)

SP2_blast_hits_seq_query = "blastdbcmd -entry_batch SP2.ids -db SP2/SP2 -out SP2_seq_best_hits_blast_1vs2.fa"
os.system(SP2_blast_hits_seq_query)

os.system("./seqkit split SP2_seq_best_hits_blast_1vs2.fa --by-part 100 --out-dir subset_SP2_Blast_reciprocal/")

os.system("ls subset_SP2_Blast_reciprocal/ > files_SP2_blast_reciprocal.txt")

os.system("mkdir result_blast_2vs1_reciprocal")

f = open("files_SP2_blast_reciprocal.txt","r")

lines = f.readlines()

i = 0

for line in lines:
    line = line.replace("\n","")
    blastn_SP1_vs_SP2 = "blastn -query subset_SP2_Blast_reciprocal/"+str(line)+" -db SP1/SP1 -outfmt 7 -subject_besthit -num_threads "+str(multiprocessing.cpu_count())+" > result_blast_2vs1_reciprocal/blast_raw_2vs1_reciprocal_0"+str(i+1)+".fa"
    os.system(blastn_SP1_vs_SP2)
    i+=1

os.system("cat result_blast_2vs1_reciprocal/* > blast_raw_2vs1_reciprocal.fna")


processing_1 = "echo \"query acc.ver\tsubject acc.ver\t% identity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. start\ts. end\tevalue\tbit score\" > temp_file.txt | cat blast_raw_2vs1_reciprocal.fna |awk '/hits found/{getline;print}' | grep -v \"#\" > blast_reciprocal_2vs1.fa"
os.system(processing_1)

processing_2 = "cat temp_file.txt blast_reciprocal_2vs1.fa > header_blast_2vs1.fa"
os.system(processing_2)

df = pd.read_csv("header_blast_2vs1.fa", sep='\t')

print(df)
column_ids_seq2 = df['subject acc.ver']
column_ids_seq2.reset_index(drop=True, inplace=True)

print(column_ids_seq2)

if os.path.isfile("reciprocal.ids"):
    pass
else:
    with open('reciprocal.ids', 'x') as f:
        column_ids_seq2.to_csv(f, sep='\t',header=False, index=False)

os.system("rm files_SP2_blast_reciprocal.txt blast_raw_2vs1_reciprocal.fna blast_best_hits_1vs2.fa temp_file.txt SP2_seq_best_hits_blast_1vs2.fa header_blast_1vs2.fa header_blast_2vs1.fa blast_reciprocal_2vs1.fa blast_raw_1vs2.fna")
os.system("rm -r subset_SP2_Blast_reciprocal result_blast_2vs1_reciprocal")
