# -*- coding: utf-8 -*-
"""
Created on Wed May  2 14:09:24 2018

@author: heinzinger
"""

from janspivot import openfile
import matplotlib.pyplot as plt

def id_venn(file):
    idps = len(file[file['Modifications'].str.len() > 1])
    dps = len(file[file['DP Modification'].str.len() > 0])
    rest = len(file) - idps - dps
    labels = 'unidentified Spectra / unknown modification', 'identified spectra / unmodified',\
             'spectra identified by dep. Peptides / specific modification'
    sizes = [rest, idps, dps]
    explode = (0, 0, 0.1)
    
    # plot
    plt.pie(sizes, explode = explode, labels = labels, shadow = True, autopct = '%1.1f%%')
    plt.axis('equal')
    plt.show()
    
if __name__ == "__main__":  
    
    id_venn(openfile())