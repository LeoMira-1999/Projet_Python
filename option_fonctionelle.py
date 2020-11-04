import ftplib #SI USER A DEJA UN FICHIER GENOME
import os
from datetime import datetime
import pandas as pd
import numpy as np
import random


FTP_HOST = "ftp.ncbi.nlm.nih.gov"
FTP_USER = "anonymous"
FTP_PASS = ""

# initialize FTP session
ftp = ftplib.FTP(FTP_HOST, FTP_USER, FTP_PASS)
# force UTF-8 encoding
ftp.encoding = "utf-8"
# print the welcome message
print(ftp.getwelcome())
# change the current working directory to refseq
ftp.cwd('genomes/refseq/')

# LIST1 a directory of bacteria
ftp.cwd('bacteria/')
print("*"*50, "LIST_1", "*"*50)
data_bac = []
ftp.dir(data_bac.append)
print(data_bac)
print("data_bac printed")

organism_bac=[]
for line in data_bac:
    organism_str_bacteria="".join(line)
    organism_line_bac=organism_str_bacteria.split(' ')
    organism_bac.append(str(organism_line_bac[-1]))
print(organism_bac)
print("organisme bac au dessus")

ftp.cwd('..')

# LIST2 a directory of archae
ftp.cwd('archaea/')
print("*"*50, "LIST_2", "*"*50)
data_archaea = []
ftp.dir(data_archaea.append)
print(data_archaea)
print("data_archaea printed")
organism_arc=[]
for line in data_archaea:
    organism_str_archaea="".join(line)
    organism_line_arc=organism_str_archaea.split(' ')
    organism_arc.append(str(organism_line_arc[-1]))
print(organism_arc)
print("organisme arc au dessus")

ftp.cwd('..')

# LIST3 a directory of vertebrate_mammalian
ftp.cwd('vertebrate_mammalian/')
print("*"*50, "LIST_3", "*"*50)
data_vertebrate_mammalian = []
ftp.dir(data_vertebrate_mammalian.append)
print(data_vertebrate_mammalian)
print("data_vertebrate_mammalian printed")

organism_vertebrate_mammalian=[]
for line in data_vertebrate_mammalian:
    organism_str_vertebrate_mammalian="".join(line)
    organism_line_vertebrate_mammalian=organism_str_vertebrate_mammalian.split(' ')
    organism_vertebrate_mammalian.append(str(organism_line_vertebrate_mammalian[-1]))
print(organism_vertebrate_mammalian)
print("organisme vertebrate_mammalian au dessus")

ftp.cwd('..')

#Fonction choix changement de direction, selection d'oganisme, +taxon
def chang_direction():
    global valeur
    valeur = entreeH.get()
    global organism_selected
    global taxon
    if (var_bacteria.get() == 1):
        ftp.cwd('bacteria/')
        organism_selected=organism_bac[organism_bac.index(valeur)]
        taxon='bacteria'
    elif (var_vertebrate_mammalian.get() == 1):
        ftp.cwd('vertebrate_mammalian/')
        organism_selected=organism_vertebrate_mammalian[organism_vertebrate_mammalian.index(valeur)]
        taxon='vertebrate_mammalian'
    else:
        ftp.cwd('archaea/')
        organism_selected=organism_arc[organism_arc.index(valeur)]
        taxon='archaea'


#fonction qui recupere la sequences sur le ncbi
data_2=[]
def get_seq():
    chang_direction()
    print("------------------")
    labelRes.configure(text=" nom de l'organisme : " + str(valeur))
    cwd_dir=str(organism_selected)+'/all_assembly_versions/'
    print(ftp.cwd(cwd_dir)) #cherche l'accession lors que l'animal est selectionné
    print(ftp.dir(data_2.append))
    print(data_2)
    print("data_2 printed")
    str_data="".join(data_2)
    print(str_data)
    accession=str_data.split('/')[-1] #prend la derniere column = code accession
    print(accession)
#    command os.system
    os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/"+taxon+"/"+organism_selected+"/all_assembly_versions/"+accession+"/"+accession+"_protein.faa.gz ")
    os.system("gunzip "+accession+"_protein.faa.gz ")
    os.system("mv "+accession+"_protein.faa "+organism_selected+"-protein.faa")
    fenetre.destroy()

from tkinter import * # création d'un menu qui va se remplir avec toute les sequences

# creation de la fenetre et configuration de son titre
fenetre = Tk()
fenetre.title (" Add genomes")

# création d’un champ text
labelH = Label(fenetre, text=" veuillez CTRL+C/CTRL+V un organisme la liste")
labelRes = Label(fenetre, text = " en attente")

# creation d’un champ de saisie
value1 = StringVar()
entreeH = Entry(fenetre, textvariable= value1, width =30)

# creation de boutons à cliquer
boutonResearch = Button(fenetre, text = "download", command = get_seq)
boutonQuitter = Button(fenetre, text = " quitter", command = fenetre.quit)

#création checkbuttun + label qui affiche la selectionné
l = Label(fenetre, bg='white', width=20, text='empty')
l.pack()

def print_selection():
    var=var_bacteria.get()+var_archae.get()+var_vertebrate_mammalian.get()
    if (var>1):
        l.config(text='trop de case coché')
    elif (var_bacteria.get() == 1):
        l.config(text='Bacteria selected ')
    elif (var_archae.get() == 1):
        l.config(text='archae selected')
    elif (var_vertebrate_mammalian.get() == 1):
        l.config(text='vertebrate_mammalian selected')
    else:
        l.config(text='nothing selected')

var_bacteria = IntVar()
var_archae = IntVar()
var_vertebrate_mammalian = IntVar()
c1 = Checkbutton(fenetre, text='bacteria',variable=var_bacteria, onvalue=1, offvalue=0, command=print_selection)
c1.pack(side=TOP, anchor=W)
c2 = Checkbutton(fenetre, text='archae',variable=var_archae, onvalue=1, offvalue=0, command=print_selection)
c2.pack(side=TOP, anchor=W)
c3 = Checkbutton(fenetre, text='vertebrate_mammalian',variable=var_vertebrate_mammalian, onvalue=1, offvalue=0, command=print_selection)
c3.pack(side=TOP, anchor=W)

# affichage des objets crées
labelH.pack()
entreeH.pack()
labelRes.pack()
boutonResearch.pack()
boutonQuitter.pack()

#création checkbox bacteria
scrollbar = Scrollbar(fenetre)
scrollbar.pack( side = RIGHT, fill = Y )

#création list avec all organisme proposition
mylist = Listbox(fenetre,width=30, height=20, yscrollcommand = scrollbar.set )
for line in range(0,len(organism_bac)):
   mylist.insert(END, "bac:"+str(organism_bac[line]))
for line in range(0,len(organism_arc)):
   mylist.insert(END, "arc:"+str(organism_arc[line]))
for line in range(0,len(organism_vertebrate_mammalian)):
   mylist.insert(END, "ve_m:"+str(organism_vertebrate_mammalian[line]))

mylist.pack( side = LEFT, fill = BOTH )
scrollbar.config( command = mylist.yview )

# lancement de la fenetre et attente des clics
fenetre.mainloop()
# destruction de la fenêtre après l’avoir quittée
fenetre.destroy()
