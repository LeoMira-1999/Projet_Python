import os #importing os to let us use bash for a later use in blast+ programs
import pandas as pd #importing pandas to read our blasted files and retreive the information needed
import glob
import ftplib
import numpy as np
import tkinter as tk
from tkinter import *
from tkinter import ttk
from main import *
from alignment import *
from option_fonctionelle import *

#__________________
def proteome_file_finder(): #function to create a list of files.faa
        faa_files =[]
        for file in glob.glob("*.faa"):
            faa_files.append(file)
        return faa_files
    #___________________

def select_genome(): #function to select_genome for blast
    faa_files_selected=[]
    global organism_dico

    for name, checkbutton in organism_dico.items():
        if checkbutton.var.get()==1:
            if checkbutton['text'] not in faa_files_selected:
                faa_files_selected.append(checkbutton['text'])

        if checkbutton.var.get()==0:
            if checkbutton['text'] in faa_files_selected:
                faa_files_selected.remove(checkbutton['text'])

    print("list of selected files.faa =>",faa_files_selected)
    return faa_files_selected

def launch():

    multi_RBH(*select_genome())

    proteomes = proteome_file_finder()

    RBH_DB_creator(proteomes)

    cluster_AC_nr, cluster_SP_nr = cluster_species_redundance_remover(RBH_analysor(RBH_comparator()), cluster_species_finder(RBH_analysor(RBH_comparator())))

    cluster_alignment(cluster_AC_nr,cluster_SP_nr)

    RBH_DB_remover(proteomes)

#_____Tkinter_____________
window = Tk()

onglet_system = ttk.Notebook(window)   # creation of tab system
onglet_system.pack()

onglet1 = ttk.Frame(onglet_system)       #add the first tab
onglet1.pack()
onglet_system.add(onglet1, text='launch')      # name of the first tab

onglet2 = ttk.Frame(onglet_system)       # add the second tab
onglet2.pack()
onglet_system.add(onglet2, text='add genome')      # name of the first tab

label = Label( onglet1, text='CHOOSE YOUR GENOMES', relief=RAISED ) #label in the first tab
label.pack()

#creation of checkbutton on the first tab
list_already=[]
organism_dico = dict()
def list_file_faa():
    for i, value in enumerate(proteome_file_finder()):
        if value not in list_already:
            organism_dico[value] = Checkbutton(onglet1, text=value, onvalue=True, offvalue=False, command=select_genome)
            organism_dico[value].var = BooleanVar(onglet1, value=False)
            organism_dico[value]['variable'] = organism_dico[value].var
            organism_dico[value].pack(padx=100, pady=1)
            list_already.append(value)

#creation of buttons the first tab
BoutonEntre = Button(onglet1, text = 'Launch', command = launch)
BoutonEntre.pack(padx=100, pady=10)
boutonRefresh = Button(onglet1, text = "Refresh", command = list_file_faa)
boutonRefresh.pack()
BoutonQuitter = Button(onglet1, text = 'Quit', command = window.destroy)
BoutonQuitter.pack(padx=100, pady=10)
list_file_faa()

#option_1_:__ADD_Genome_From_The_Web
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

# directory of bacteria________________________________________________________
ftp.cwd('bacteria/')
print("*"*50, "addition of the list of bacteria from refseq", "*"*50)
data_bac = []
ftp.dir(data_bac.append)

#creation of name list of bacteria
organism_bac=[]
for line in data_bac:
    organism_str_bacteria="".join(line)
    organism_line_bac=organism_str_bacteria.split(' ')
    organism_bac.append(str(organism_line_bac[-1]))

ftp.cwd('..') # return to directory of refseq

# directory of archae
ftp.cwd('archaea/')
print("*"*50, "addition of the list of archaea from refseq", "*"*50)
data_archaea = []
ftp.dir(data_archaea.append)

#creation of name list of archaea
organism_arc=[]
for line in data_archaea:
    organism_str_archaea="".join(line)
    organism_line_arc=organism_str_archaea.split(' ')
    organism_arc.append(str(organism_line_arc[-1]))

ftp.cwd('..')# return to directory of refseq

# directory of vertebrate_mammalian
ftp.cwd('vertebrate_mammalian/')
print("*"*50, "addition of the list of vertebrate_mammalian from refseq", "*"*50)
data_vertebrate_mammalian = []
ftp.dir(data_vertebrate_mammalian.append)

#creation of name list of vertebrate_mammalian
organism_vertebrate_mammalian=[]
for line in data_vertebrate_mammalian:
    organism_str_vertebrate_mammalian="".join(line)
    organism_line_vertebrate_mammalian=organism_str_vertebrate_mammalian.split(' ')
    organism_vertebrate_mammalian.append(str(organism_line_vertebrate_mammalian[-1]))

ftp.cwd('..')# return to directory of refseq


#Fonction choix changement de direction, selection d'oganisme, +taxon
def directory_choice():
    global valeur
    global valeur_temp
    valeur_temp = entree_organismes.get()# (exemple = bac:salmonella_bongori)
    valeur = valeur_temp.split(":", 1)

    #the word after ":" is used to change the name of organism_selected
    global organism_selected
    organism_selected=valeur[-1]
    tax_id=valeur[0]
    global taxon

    #the word before ":" is used to change the label text, the directory and the taxon
    if (tax_id == 'bac'):
        l.config(text='Bacteria selected ')
        ftp.cwd('bacteria/')
        taxon='bacteria'

    elif (tax_id == 've_m'):
        l.config(text='vertebrate_mammalian selected')
        ftp.cwd('vertebrate_mammalian/')
        taxon='vertebrate_mammalian'

    else:
        l.config(text='archae selected')
        ftp.cwd('archaea/')
        taxon='archaea'

#function find and DDL the protein.faa from refseq
data_2=[]
def get_seq():
    directory_choice()
    labelRes.configure(text=" nom de l'organisme : " + str(organism_selected))
    cwd_dir=str(organism_selected)+'/all_assembly_versions/'
    print(ftp.cwd(cwd_dir)) #when the organism is selected, find the accession code
    print(ftp.dir(data_2.append))
    str_data="".join(data_2)
    accession=str_data.split('/')[-1] #last column = accession code
#    command os.system
    os.system("wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/"+taxon+"/"+organism_selected+"/all_assembly_versions/"+accession+"/"+accession+"_protein.faa.gz ")
    os.system("gunzip "+accession+"_protein.faa.gz ")
    os.system("mv "+accession+"_protein.faa "+(organism_selected.replace("_", "-"))+"-protein.faa")
    print(ftp.cwd('..'))# get out from the directory
    print(ftp.cwd('..'))
    print(ftp.cwd('..'))


# creation of a text field
labelnote = Label(onglet2, text="note: please CTRL+C/CTRL+V a line from the list to download \n \n to search :")
labelRes = Label(onglet2, text = " waiting")

# creation of an input field
value_organism = StringVar()
entree_organismes = Entry(onglet2, textvariable= value_organism, width =30)

# creation of the DDL and quit buttons
boutonResearch = Button(onglet2, text = "download", command = get_seq)
boutonQuitter = Button(onglet2, text = " Quit", command = window.quit)

#label that displays the tax_id => ex: bacteria is selected
l = Label(onglet2, bg='white', width=20, text='empty')

# display of created objects
l.pack()
labelRes.pack()
entree_organismes.pack()
boutonResearch.pack()
boutonQuitter.pack()
labelnote.pack()

#boxlist of all selectable organism:
class Organism_refseq(Frame):
    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.pack()
        self.create_widgets()

# Create main GUI window
    def create_widgets(self):
        self.search_var = StringVar()
        self.search_var.trace("w", self.update_list)
        self.entry = Entry(onglet2, textvariable=self.search_var, width=30)
        self.lbox = Listbox(onglet2, width=45, height=15)

        self.entry.pack(padx=100, pady=1)
        self.lbox.pack(padx=100, pady=1)
    # Function for updating the list/doing the search.
    # It needs to be called here to populate the listbox.
        self.update_list()
    def update_list(self, *args):
        search_term = self.search_var.get()

    # Just a generic list to populate the listbox
        lbox_list = []
        organism_bac_head = ["bac:" + item for item in organism_bac]
        lbox_list.extend(organism_bac_head)

        organism_arc_head = ["arc:" + item for item in organism_arc]
        lbox_list.extend(organism_arc_head)

        organism_vertebrate_mammalian_head = ["ve_m:" + item for item in organism_vertebrate_mammalian]
        lbox_list.extend(organism_vertebrate_mammalian_head)

        self.lbox.delete(0, END)

        for item in lbox_list:
            if search_term.lower() in item.lower():
                self.lbox.insert(END, item)



Organism_refseq() #call the class

window.mainloop() #call the window Tkinter
