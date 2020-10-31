import os

S = "blastn -query GCF_001402945.1_ASM140294v1_genomic.fna -subject GCF_001402935.1_ASM140293v1_genomic.fna -outfmt 6 -subject_besthit"

y = os.system(S)
print("--------------------------")
print(type(y))

for x in y:
    print(x)

atom://teletype/portal/87baca77-8b58-4363-8370-460f61c226cd
