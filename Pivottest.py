# -*- coding: utf-8 -*-
"""
Created on Fri Apr  6 13:56:14 2018

@author: heinzinger
"""

import time as time 
from perseuspy import pd
from tkinter import Tk
from tkinter import filedialog


if __name__ == "__main__":    
    
    start_time = time.time();
    root = Tk()
    root.withdraw()
    a = filedialog.askopenfilename(initialdir = "D:\PhytonSCripts\MQOutput" \
           ,title = "Choose a MQ DP_peptides.deppep file to plot",\
           filetypes = (("deppep files","*.deppep"),("all files","*.*")))
    a = pd.read_table(a,low_memory=False)
    a = a.drop(0).reset_index(drop = True)

    a['Intensity'] = a['Intensity'].astype(float)

    table = pd.pivot_table(a, values = 'Intensity',\
    index = ['DP Cluster Index','DP Base Sequence'], columns = 'Raw file')
    
    table = table.reset_index()

    print("---Runtime =  %s seconds ---" % (time.time() - start_time))   