# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 10:31:43 2018

@author: heinzinger
"""
from utilities import openfile
from perseuspy import pd

# plot different dependent Peptide against each other to see their abundancy 
# of found dependent peptides 
def pltcountperRaw():
    x = int(input('How many files do you want to plot together? '))
  
    a = openfile()
    withoutunmodified = a[a['DP Modification'] == 'Unmodified']
    unmod = pd.DataFrame(withoutunmodified['Raw file'].value_counts())
    values = pd.DataFrame(a['Raw file'].value_counts())  
    
# =============================================================================
#     #just for MSFragger files
#     b = openfile()
#     b['Raw file'] = b['Spectrum'].apply(lambda df: df.split('.')[0])
#     b = b[b['Observed Modifications']!='Unknown']
#     values['MSFragger identifications'] = b['Raw file'].value_counts().values
# =============================================================================
    for y in range(x-1):
        a = openfile()
        withoutunmodified = a[a['DP Modification'] == 'Unmodified']
        values['File '+ str(y)] = a['Raw file'].value_counts().values
        unmod['File '+ str(y)] = withoutunmodified['Raw file'].value_counts().values

    values = values.sort_index() 
    unmod = unmod.sort_index() 
    values = values.reset_index(drop = True)
    unmod = unmod.reset_index(drop = True)
    if x == 3:
        values = values.rename(index = str, columns = {'Raw file':'No matching',\
                                'File 0':'+/-1 fraction matching',\
                                'File 1':'Complete matching'})
        unmod = unmod.rename(index = str, columns = {'Raw file':'No matching',\
                                'File 0':'Matching restricted to one raw file',\
                                'File 1':'Matching between all raw files'})
    
    ax = values.plot.bar(rot=0,fontsize = 20)
    #unmod.plot.bar(ax=ax,color=['cornflowerblue', 'orange', 'limegreen'],fontsize = 15,legend = False,rot=0)
    ax.legend(fontsize = 20)
    ax.set_xticklabels(list(range(1,len(values)+1)))
    ax.set_xlabel('Raw file',fontsize = 20)
    ax.set_ylabel('Count',fontsize = 20)
    return unmod
   
if __name__ == '__main__':
    unmod = pltcountperRaw()