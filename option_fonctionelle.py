import ftplib
import os
from datetime import datetime
import pandas as pd
import numpy as np
import random
from tkinter import * # création d'un menu qui va se remplir avec toute les sequences

def ftp_refseq_proteome_finder():
    FTP_HOST = "ftp.ncbi.nlm.nih.gov"
    FTP_USER = "anonymous"
    FTP_PASS = ""

    # initialize FTP session
    ftp = ftplib.FTP(FTP_HOST, FTP_USER, FTP_PASS)
    # force UTF-8 encoding
    ftp.encoding = "utf-8"
    # print the welcome message of NCBI
    print(ftp.getwelcome())
    # change the current working directory to refseq
    ftp.cwd('genomes/refseq/')

    # directory of bacteria to
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
        global valeur_temp
        valeur_temp = entreeH.get()
        valeur = valeur_temp.split(":", 1)
        global organism_selected
        global taxon
        if (valeur[0] == 'bac'):
            l.config(text='Bacteria selected ')
            ftp.cwd('bacteria/')
            organism_selected=organism_bac[organism_bac.index(valeur[-1])]
            print()
            taxon='bacteria'
        elif (valeur[0] == 've_m'):
            l.config(text='vertebrate_mammalian selected')
            ftp.cwd('vertebrate_mammalian/')
            organism_selected=organism_vertebrate_mammalian[organism_vertebrate_mammalian.index(valeur[-1])]
            taxon='vertebrate_mammalian'
        else:
            l.config(text='archae selected')
            ftp.cwd('archaea/')
            organism_selected=organism_arc[organism_arc.index(valeur[-1])]
            taxon='archaea'


    #fonction qui recupere la sequences sur le ncbi
    data_2=[]
    def get_seq():
        chang_direction()
        print("------------------")
        labelRes.configure(text=" nom de l'organisme : " + str(valeur[-1]))
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
        os.system("mv "+accession+"_protein.faa "+(organism_selected.replace("_", "-"))+"-protein.faa")
        print(ftp.cwd('..'))# on ressors du directory pour revenir a refseq
        print(ftp.cwd('..'))
        print(ftp.cwd('..'))
        print(ftp.pwd())




    # creation de la fenetre et configuration de son titre
    fenetre = Tk()
    fenetre.title (" Add genomes")


    # création d’un champ text
    labelH = Label(fenetre, text=" veuillez CTRL+C/CTRL+V un organisme de la liste")
    labelRes = Label(fenetre, text = " en attente")

    # creation d’un champ de saisie
    value1 = StringVar()
    entreeH = Entry(fenetre, textvariable= value1, width =30)

    # creation de boutons à cliquer
    boutonResearch = Button(fenetre, text = "download", command = get_seq)
    boutonQuitter = Button(fenetre, text = " quitter", command = fenetre.quit)

    #création label qui affiche la selection
    l = Label(fenetre, bg='white', width=20, text='empty')

    # affichage des objets crées
    l.pack()
    labelRes.pack()
    entreeH.pack()
    boutonResearch.pack()
    boutonQuitter.pack()
    labelH.pack()

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
