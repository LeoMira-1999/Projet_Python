import pandas as pd
import os

SP1 = ""

SP2 = ""

blastn_SP1_vs_SP2 = "blastn -query "+SP1+" -subject "+SP2+" -outfmt 7 -subject_besthit > blast_raw_1vs2.fa"
os.system(blastn_SP1_vs_SP2)

processing_1 = "echo \"query acc.ver\tsubject acc.ver\t% identity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. start\ts. end\tevalue\tbit score\" > temp_file.txt | cat blast_raw_1vs2.fa |awk '/hits found/{getline;print}' | grep -v \"#\" > blast_best_hits_1vs2.fa"
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

SP2_DB_creation = "makeblastdb -in "+SP2+" -parse_seqids -blastdb_version 5 -dbtype nucl -out SP2/SP2"
os.system(SP2_DB_creation)

SP2_blast_hits_seq_query = "blastdbcmd -entry_batch SP2.ids -db SP2/SP2 -out SP2_seq_best_hits_blast_1vs2.fa"
os.system(SP2_blast_hits_seq_query)

blastn_SP2_best_hits_vs_SP1 = "blastn -query SP2_seq_best_hits_blast_1vs2.fa -subject GCA_004765815.2_ASM476581v2_genomic.fna -outfmt 7 -subject_besthit > blast_raw_2BH_vs_1.fa"
os.system(blastn_SP2_best_hits_vs_SP1)

processing_1 = "echo \"query acc.ver\tsubject acc.ver\t% identity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. start\ts. end\tevalue\tbit score\" > temp_file.txt | cat blast_raw_2BH_vs_1.fa |awk '/hits found/{getline;print}' | grep -v \"#\" > blast_best_hits_2vs1.fa"
os.system(processing_1)

processing_2 = "cat temp_file.txt blast_best_hits_2vs1.fa > header_blast_2vs1.fa | rm temp_file.txt"
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
