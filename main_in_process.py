import os
import multiprocessing


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
