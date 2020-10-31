import os

S = "blastn -query GCF_001402945.1_ASM140294v1_genomic.fna -subject GCF_001402935.1_ASM140293v1_genomic.fna -outfmt 6 -subject_besthit"

y = os.system(S)
print("--------------------------")
print(type(y))

for x in y:
    print(x)

atom://teletype/portal/f45ff42d-d108-443b-be8e-20faa482e980
