import ftplib
import os
from datetime import datetime
import pandas as pd
import numpy as np

os.system("rm list_genome_add.txt")

FTP_HOST = "ftp.ncbi.nlm.nih.gov"
FTP_USER = "anonymous"
FTP_PASS = ""

# some utility functions that we gonna need
def get_size_format(n, suffix="B"):
    # converts bytes to scaled format (e.g KB, MB, etc.)
    for unit in ["", "K", "M", "G", "T", "P"]:
        if n < 1024:
            return f"{n:.2f}{unit}{suffix}"
        n /= 1024

def get_datetime_format(date_time):
    # convert to datetime object
    date_time = datetime.strptime(date_time, "%Y%m%d%H%M%S")
    # convert to human readable date time string
    return date_time.strftime("%Y/%m/%d %H:%M:%S")

# initialize FTP session
ftp = ftplib.FTP(FTP_HOST, FTP_USER, FTP_PASS)
# force UTF-8 encoding
ftp.encoding = "utf-8"
# print the welcome message
print(ftp.getwelcome())
# change the current working directory to 'pub' folder and 'maps' subfolder
ftp.cwd('genomes/refseq/bacteria/Acetanaerobacterium_elongatum/all_assembly_versions/GCF_900103835.1_IMG-taxon_2667527408_annotated_assembly')
# LIST a directory
print("*"*50, "LIST", "*"*50)
data = []
ftp.dir((data.append))
print(data)
#__________________________tentative____________________________
from tkinter import * # création d'un menu qui va se remplir avec toute les sequences
fen =Tk()
fen.geometry("200x200")
mb =Menubutton(fen, text="Bacteria", relief=RAISED )
mb.grid()
mb.menu=Menu ( mb, tearoff = 0)
mb["menu"]=mb.menu

df = pd.DataFrame() #création d'un dataframe qui remplira avec toutes les sequences

def OS_Wget():
    os.system("wget --timestamping ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Acetobacter_cerevisiae/latest_assembly_versions/GCF_001580535.1_ASM158053v1/GCF_001580535.1_ASM158053v1_genomic.fna.gz")
    os.system("gunzip GCF_001580535.1_ASM158053v1_genomic.fna.gz ")
with open('list_genome_add.txt', 'w') as f: #fichier qui se remplira avec toutes les sequences
#la boucle qui remplie
    for y in data:
        f.write("%s\n" % str(y))
        df = df.append({'gene': y}, ignore_index=True)
        mb.menu.add_checkbutton ( label=str(y), command= OS_Wget)
f.close()

#print(df)

mb.pack()
fen.mainloop()
