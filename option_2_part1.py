import os #importing os to let us use bash for a later use in blast+ programs
import pandas as pd #importing pandas to read our blasted files and retreive the information needed
import glob
import ftplib
from datetime import datetime
import numpy as np
import random
from tkinter import *
#__________________
def proteome_file_finder():
        faa_files =[]
        for file in glob.glob("*.faa"):
            faa_files.append(file)
            print(faa_files)
        return faa_files
    #___________________

def select_genome():
    faa_files_selected=[]
    global cbs
    for name, checkbutton in cbs.items():
        if checkbutton.var.get()==1:
            if checkbutton['text'] not in faa_files_selected:
                faa_files_selected.append(checkbutton['text'])
        if checkbutton.var.get()==0:
            if checkbutton['text'] in faa_files_selected:
                faa_files_selected.remove(checkbutton['text'])
    print(faa_files_selected)
    return faa_files_selected

fenetre = Tk()
text = Text(fenetre)
label = Label( fenetre, text='CHOICE YOUR GENOMES', relief=RAISED )
label.grid()


cbs = dict()
for i, value in enumerate(proteome_file_finder()):
    cbs[value] = Checkbutton(fenetre, text=value, onvalue=True, offvalue=False, command=select_genome)
    cbs[value].var = BooleanVar(fenetre, value=False)
    cbs[value]['variable'] = cbs[value].var

    cbs[value].grid(row=i+2, column=0)

BoutonEntre = Button(fenetre, text = 'LANCER', command = fenetre.destroy)
BoutonEntre.grid(row=0, column=1)
BoutonQuitter = Button(fenetre, text = 'QUITTER', command = fenetre.destroy)
BoutonQuitter.grid(row=1, column=1)

text.grid()
fenetre.mainloop()
