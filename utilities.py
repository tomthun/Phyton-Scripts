# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 11:44:43 2018

@author: heinzinger
"""
from perseuspy import pd
from tkinter import Tk
from tkinter import filedialog

def openfile():
    root = Tk()
    root.withdraw()
    a = filedialog.askopenfilename(initialdir = "D:\Thomas\PhytonSCripts\MQOutput" \
           ,title = "Choose a MQ .deppep file or .tsv or .txt to plot",\
           filetypes = (("all files","*.*"),("deppep files","*.deppep"),("tsv files",".tsv"),\
                        ("MQ Output",".txt")))
   
    table = pd.read_table(a,low_memory=False)
    if ".deppep" in a:
        table = table.drop(0).reset_index(drop = True)
    if 'Intensity' in table.columns:
        table['Intensity'] = table['Intensity'].astype(float)
    return table