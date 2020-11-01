import ftplib
import os
from datetime import datetime
import pandas as pd
import numpy as np

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


df = pd.DataFrame(columns=['A'])
f = open("list_genome_add2.txt", "x")
with open('list_genome_add2.txt', 'w') as f:
    for y in data:
        f.write("%s\n" % str(y))
        df = df.append({'A': y}, ignore_index=True)
f.close()

print("*"*50, "df", "*"*50)

# dropping null value columns to avoid errors
df.dropna(inplace = True)
# new data frame with split value columns
new = df["A"].str.split("%\n", expand = True)
# making separate name column from new data frame
print(df['A'])

print("*"*50, "on est la", "*"*50)
print(new)
