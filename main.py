import pandas as pd
import os

blastn_SP1_vs_SP2 = "blastn -query GCF_001402945.1_ASM140294v1_genomic.fna -subject GCF_001402935.1_ASM140293v1_genomic.fna -outfmt 7 -subject_besthit > blast_raw_1vs2.fa"
os.system(blastn_SP1_vs_SP2)

processing_1 = "touch temp_file.txt | echo \"query acc.ver\tsubject acc.ver\t% identity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. start\ts. end\tevalue\tbit score\" > temp_file.txt | cat blast_raw_1vs2.fa |awk '/hits found/{getline;print}' | grep -v \"#\" > blast_best_hits_1vs2.fa"
os.system(processing_1)

processing_2 = "cat temp_file.txt blast_best_hits_1vs2.fa > header_blast.fa"
os.system(processing_2)

df = pd.read_csv("header_blast.fa", sep='\t')

print(df)
column_ids_seq2 = df['subject acc.ver']
column_ids_seq2.reset_index(drop=True, inplace=True)

print(column_ids_seq2)

if os.path.isfile("SP2.ids"):
    pass
else:
    with open('SP2.ids', 'x') as f:
        column_ids_seq2.to_csv(f, sep='\t',header=False, index=False)

SP2_DB_creation = "makeblastdb -in GCF_001402935.1_ASM140293v1_genomic.fna -parse_seqids -blastdb_version 5 -dbtype nucl -out SP2/SP2"
os.system(SP2_DB_creation)

SP2_blast_hits_seq_query = "blastdbcmd -entry_batch SP2.ids -db SP2/SP2 -out SP2_seq_best_hits_blast_1vs2.fa"
os.system(SP2_blast_hits_seq_query)

blastn_SP2_best_hits_vs_SP1 = "blastn -query SP2_seq_best_hits_blast_1vs2.fa -subject GCF_001402945.1_ASM140294v1_genomic.fna -outfmt 7 -subject_besthit > blast_raw_2BH_vs_1.fa"
os.system(blastn_SP2_best_hits_vs_SP1)

processing_1 = "rm temp_file.txt | touch temp_file.txt | echo \"query acc.ver\tsubject acc.ver\t% identity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. start\ts. end\tevalue\tbit score\" > temp_file.txt | cat blast_raw_1vs2.fa |awk '/hits found/{getline;print}' | grep -v \"#\" > blast_best_hits_1vs2.fa"
os.system(processing_1)

processing_2 = "cat temp_file.txt blast_best_hits_1vs2.fa > header_blast.fa"
os.system(processing_2)

df = pd.read_csv("header_blast.fa", sep='\t')

print(df)
column_ids_seq2 = df['subject acc.ver']
column_ids_seq2.reset_index(drop=True, inplace=True)

print(column_ids_seq2)

if os.path.isfile("SP2.ids"):
    pass
else:
    with open('SP2.ids', 'x') as f:
        column_ids_seq2.to_csv(f, sep='\t',header=False, index=False)
