import pandas as pd
import os

if os.path.isfile("test.fa"):
    pass
else:
    S = "blastn -query GCF_001402945.1_ASM140294v1_genomic.fna -subject GCF_001402935.1_ASM140293v1_genomic.fna -outfmt 6 -subject_besthit > test.fa"
    y = os.system(S)
    print("--------------------------")
    print(type(y))

df = pd.read_csv("test.fa", sep='\t')

print(df)

atom://teletype/portal/d111786f-a868-4a3e-a5f9-cd9eafb2f764
