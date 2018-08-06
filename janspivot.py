# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 14:01:19 2018

@author: heinzinger
"""
import time as time 
from perseuspy import pd
from tkinter import Tk
from tkinter import filedialog
import numpy as np
def janspivot(a,reg):
    if reg:
        table = pd.pivot_table(a, values = 'Intensity',\
        index = ['DP Cluster Index','DP AA','DP Base Sequence','DP Probabilities',\
                 'DP Probabilities after Regression'], columns = 'Raw file')
        table = table.reset_index()
        grouped = table.groupby(['DP Cluster Index', 'DP Base Sequence'])
        probs = grouped['DP Probabilities'].apply(lambda df: ';'.join(df))
        means = grouped.apply(lambda df: df[table.columns[4:]].mean())
        probsreg = grouped['DP Probabilities after Regression'].apply(lambda df: ';'.join(df))  
        pivot = pd.concat([probs,probsreg,means], axis=1)
        pivot = pivot.reset_index()
    else:
        table = pd.pivot_table(a, values = 'Intensity',\
        index = ['DP Cluster Index','DP AA','DP Base Sequence','DP Probabilities'], columns = 'Raw file')
        table = table.reset_index()
        grouped = table.groupby(['DP Cluster Index', 'DP Base Sequence'])
        probs = grouped['DP Probabilities'].apply(lambda df: ';'.join(df))
        means = grouped.apply(lambda df: df[table.columns[4:]].mean())
        pivot = pd.concat([probs,means], axis=1)
        pivot = pivot.reset_index()
    return pivot

def MSFrgPivot(a):
    a['Raw file'] = a['Spectrum'].apply(lambda df: df.split('.')[0])
    a[a['Observed Modifications']!='Unknown']
    table = pd.pivot_table(a, values = 'Calculated M/Z',index = ['Peptide','Observed Modifications'], columns = 'Raw file')
    pivot = table.reset_index()    
    return pivot

def comparePivot(fileendingdeppep):
    a = openfile()
    b = openfile()
    x = 3
    if fileendingdeppep:
        piva = janspivot(a,False)
        pivb = janspivot(b,False)
    else:
        piva = janspivot(a,False)
        pivb = MSFrgPivot(b)
        x - 1
    countsa = []
    countsb = []
    for row in piva[piva.columns[3:]].itertuples():
        countsa.append( np.count_nonzero(np.nan_to_num(row[1:])))
    for row in pivb[pivb.columns[x:]].itertuples():
        countsb.append( np.count_nonzero(np.nan_to_num(row[1:])))
    
    countsa = pd.Series(countsa).value_counts()
    countsb = pd.Series(countsb).value_counts()
    countsall = pd.concat([countsa,countsb],axis=1)
    countsall = countsall[countsall[0]>10]
    countsall[0] = (countsall[0]/sum(countsall[0]))*100
    countsall[1] = (countsall[1]/sum(countsall[1]))*100
    if fileendingdeppep:
        countsall = countsall.rename(index = str, columns = {0:'No matching',\
                                 1:'Complete matching'})
    else:       
        countsall = countsall.rename(index = str, columns = {0:'MQ Complete matching',\
                                 1:'MSFragger open search'})
    ax = countsall.plot.bar(rot=0,fontsize = 15)
    ax.legend(fontsize = 15)
    ax.set_xlabel('Raw file',fontsize = 15)
    ax.set_ylabel('Percentage share',fontsize = 15)
    return countsall

def openfile():
    root = Tk()
    root.withdraw()
    a = filedialog.askopenfilename(initialdir = "D:\Thomas\PhytonSCripts\MQOutput" \
           ,title = "Choose a MQ DP_peptides.deppep file or .tsv or .txt to plot",\
           filetypes = (("deppep files","*.deppep"),("tsv files",".tsv"),\
                        ("MQ Output",".txt"),("all files","*.*")))
   
    table = pd.read_table(a,low_memory=False)
    if ".deppep" in a:
        table = table.drop(0).reset_index(drop = True)
    if 'Intensity' in table.columns:
        table['Intensity'] = table['Intensity'].astype(float)
    return table

if __name__ == "__main__":
    start_time = time.time();
    counts = comparePivot(False)
    print("---Runtime =  %s seconds ---" % (time.time() - start_time))   
