# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 16:07:07 2018

@author: heinzinger
"""
# =============================================================================
# Imports
# =============================================================================
import time as time 
from perseuspy import pd
from tkinter import Tk
from tkinter import filedialog
# =============================================================================
# Tools to extract the rows which to puttogether
# =============================================================================

def extractRowsandPivot(table):

    dupl = table['DP Cluster Index'].value_counts()>1
    dupl = dupl[dupl == True]
    for x in range(len(dupl)):
        table = extractdupls(dupl,x,table)   
 
    return table           
        
def extractdupls(dupl,x,table):            
        dupl = dupl[dupl == True]
        snip = table.loc[(table['DP Cluster Index'] == dupl.index[x])]
        duplsnip = snip['DP Base Sequence'].value_counts()>1
        duplsnip = duplsnip[duplsnip == True]
        for z in range(len(duplsnip)):
            seqs = snip.loc[(snip['DP Base Sequence'] == duplsnip.index[z])]
            # do work with sequences
            table = extractInformation(seqs,table)
        return table
    
def extractInformation(seqs,table):
    
    # new ordering
    newindex = seqs['DP Probabilities'].str.len().sort_values(ascending=False).index
    seqs = seqs.reindex(newindex)
    seqs = seqs.reset_index()
    base = seqs['DP Base Sequence'][0]
    big = seqs['DP Probabilities'][0]
    chars = set(seqs['DP AA'][0].replace(';',''))
            
    for x in range(len(seqs)-1):
        mods = seqs['DP AA'][x+1]
        small = seqs['DP Probabilities'][x+1]
        if any((c in chars) for c in mods):
            big = puttogether(base,big,small)
        else:
            seqs.drop(x)
            
    toinsert = seqs.loc[0]
    toinsert[4:] = seqs.mean()[1:]
    toinsert['DP Probabilities'] = big
    toinsert = toinsert.drop('index')
    seqs = seqs.set_index(seqs['index']) 
    seqs = seqs.drop(columns = 'index')
    
    table = table.drop(seqs.index)  
    table = table.append(toinsert)
    return table
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
    merged = merged.sort_values(by = [1])
    merged = merged.reset_index(drop = True)
    for x in range(len(merged)):
        x = len(merged) - x - 1
        pbig = merged['0_x']
        psmall = merged['0_y']       
        if float(psmall[x]) > 0 and not pbig.isnull()[x]:
            pbig[x] = str(pbig[x]) +','+ str(psmall[x])
        elif float(psmall[x]) > 0:
            pbig[x] = str(psmall[x]) 
        value = '('+pbig[x]+')'
        index = merged[1][x]
        orgstr = insert_str(orgstr,value,index+1)
                
    return orgstr

# =============================================================================
# Start a custom pivot of a MaxQuant .deppep file here
# =============================================================================
    
def CustomPivot():
        
    start_time = time.time();
    root = Tk()
    root.withdraw()
    a = filedialog.askopenfilename(initialdir = "D:\Thomas\PhytonSCripts\MQOutput" \
           ,title = "Choose a MQ DP_peptides.deppep file to plot",\
           filetypes = (("deppep files","*.deppep"),("all files","*.*")))
    a = pd.read_table(a,low_memory=False)
    a = a.drop(0).reset_index(drop = True)
    
    a['Intensity'] = a['Intensity'].astype(float)
    
    table = pd.pivot_table(a, values = 'Intensity',\
    index = ['DP Cluster Index','DP AA','DP Base Sequence','DP Probabilities'], columns = 'Raw file')
    
    table = table.reset_index()
    
    # magic
    table = extractRowsandPivot(table)
    table = table.sort_values(by = ['DP Cluster Index'])
    table = table.reset_index(drop=True)
    print("---Runtime =  %s minutes ---" % (time.time() - start_time))   

    return table

if __name__ == "__main__":    

    output = CustomPivot()
