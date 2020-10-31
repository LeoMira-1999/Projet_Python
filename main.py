from Bio.Blast import NCBIWWW
help(NCBIWWW.qblast)
fasta_string1 = open("GCF_001402945.1_ASM140294v1_genomic.fna").read()
fasta_string2 = open("GCF_001402935.1_ASM140293v1_genomic.fna").read()
result_handle = NCBIWWW.qblast("blastn", database = "nr"  ,sequence=fasta_string1)

atom://teletype/portal/d111786f-a868-4a3e-a5f9-cd9eafb2f764
