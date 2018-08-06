# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 16:07:07 2018

@author: heinzinger
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from perseuspy import pd
from tkinter import Tk
from tkinter import filedialog
from scipy.stats import kde

# =============================================================================
# Tools to extract the rows which to puttogether
# =============================================================================

def extractRows(table):
    dupl = table['DP Cluster Index'].value_counts()>1
    for x in range(len(dupl)):
        dupl = dupl[dupl == True]
        snip = table.loc[(table['DP Cluster Index'] == dupl.index[x])]
        for y in range(len(snip)):
            duplsnip = snip['DP Base Sequence'].value_counts()>1
            duplsnip = duplsnip[duplsnip == True]
            seqs = snip.loc[(snip['DP Base Sequence'] == duplsnip.index[y])]
            # do work with sequences
            table.drop(seqs.index)
            newindex = seqs['DP Probabilities'].str.len().sort_values(ascending=False).index
            seqs.reindex(newindex)
            extractInformation(seqs)
            
def extractInformation(seqs):
    
               
            
# =============================================================================
#     puttogether(inserttostr,str1,str2):
#         Tool to put a String between two Chars together
#         e.g 'AB(0.04)CD(0.03)E(0.2)FG(0.2)A(0.005)' &
#             'AB(0.94)CDE(0.5)FG(0.5)A'
#           = 'AB(0.04,0.94)CD(0.03)E(0.2,0.5)FG(0.2,0.5)A(0.005)'
# =============================================================================

def insert_str(string, str_to_insert, index):
    return string[:index] + str_to_insert + string[index:]
    
def findallin (cache,position,thestr,char1,char2):
#    initialize cache,postiton as []        
  
    if(thestr.find(char1) > -1 and thestr.find(char2) > -1):    
        substr = thestr[thestr.find(char1)+1:thestr.find(char2)]
        subsuperstr = thestr[thestr.find(char2)+1:len(thestr)]
        cache.append(substr)
        
        if not position: 
            position.append(thestr.find(char1)-1)
        else:
            idex = position[len(position)-1] + thestr.find(char1)
            position.append(idex)
        findallin(cache,position,subsuperstr,char1,char2)
        
    df = pd.DataFrame(cache)
    df.insert(1,1,position)
    return df
        
def puttogether(orgstr,bigstr,smallstr):
    df1 = findallin ([],[],bigstr,'(',')')
    df2 = findallin ([],[],smallstr,'(',')')

    merged = pd.merge(df1,df2,how = 'outer',on=[1])
    merged['0_y'].astype(str)
    for x in range(len(merged)):
        x = len(merged) - x - 1
        if float(merged['0_y'][x]) > 0:
            merged['0_x'][x] = str(merged['0_x'][x]) +','+ str(merged['0_y'][x])
        
        value = '('+merged['0_x'][x]+')'
        index = merged[1][x]
        orgstr = insert_str(orgstr,value,index+1)
                
    return orgstr

# =============================================================================
# Start a custom pivot of a MaxQuant .deppep file here
# =============================================================================
    
def CustomPivot():
    root = Tk()
    root.withdraw()
    a = filedialog.askopenfilename(initialdir = "D:\PhytonSCripts\MQOutput" \
           ,title = "Choose a MQ DP_peptides.deppep file to plot",\
           filetypes = (("deppep files","*.deppep"),("all files","*.*")))
    a = pd.read_table(a,low_memory=False)
    a = a.drop(0).reset_index(drop = True)
    
    b = filedialog.askopenfilename(initialdir = "D:\PhytonSCripts\MQOutput" \
           ,title = "Choose a MQ DP_peptides.deppep file to plot",\
           filetypes = (("deppep files","*.deppep"),("all files","*.*")))
    b = pd.read_table(b,low_memory=False)
    b = b.drop(0).reset_index(drop = True)

    a['Intensity'] = a['Intensity'].astype(float)
    b['Intensity'] = b['Intensity'].astype(float)

    raw_tablea = pd.pivot_table(a, values = 'Intensity',\
    index = ['DP Cluster Index','DP AA','DP Base Sequence'], columns = 'Raw file')
    raw_tableb = pd.pivot_table(b, values = 'Intensity',\
    index = ['DP Cluster Index','DP AA','DP Base Sequence'], columns = 'Raw file')
    
    raw_tablea = raw_tablea.reset_index()
    raw_tableb = raw_tableb.reset_index()
    
    # magic
    extractRows(raw_tablea)
    
    return raw_tablea,raw_tableb


