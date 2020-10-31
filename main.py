import pandas as pd
import os


firstCline = "blastn -query GCF_001402945.1_ASM140294v1_genomic.fna -subject GCF_001402935.1_ASM140293v1_genomic.fna -outfmt 7 -subject_besthit > blast_raw_1vs2.fa"
execute = os.system(firstCline)

secondCline = "touch temp_file.txt | echo \"query acc.ver\tsubject acc.ver\t% identity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. start\ts. end\tevalue\tbit score\" > temp_file.txt | cat blast_raw_1vs2.fa |awk '/hits found/{getline;print}' | grep -v \"#\" > blast_best_hits_1vs2.fa"
execute = os.system(secondCline)

thirdCline = "cat temp_file.txt blast_best_hits_1vs2.fa > header_blast.fa | rm temp_file.txt "
execute = os.system(thirdCline)

df = pd.read_csv("header_blast.fa", sep='\t')

print(df)
column_ids_seq2 = df['subject acc.ver']

print(column_ids_seq2)

with open('SP2.ids', 'x') as f:
column_ids_seq2.to_csv(f, sep='\t',header=False)
