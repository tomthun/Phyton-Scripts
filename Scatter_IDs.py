# -*- coding: utf-8 -*-
"""
Created on Wed May 23 13:39:26 2018

@author: heinzinger
"""
from janspivot import openfile
import time as time 
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# =============================================================================
# %load_ext Cython
# %cython
# =============================================================================
def findrettime(raws,featureids,allpeptides):     
    # find matching retention time to raw and multiid
    rettime = allpeptides.loc[(int(raws),int(featureids)),'Retention time']
    return rettime
    
def preprocessing():
    # openfiles
    matchedFeatures = openfile()  
    allpeptides = openfile()    
    
    allpeptides = reshape(allpeptides)
    matchedFeaturesexpanded = pir2(matchedFeatures,'Raw files','Multiplet ids')
    matchedFeaturesexpanded['allPeptides retention time'] = \
        matchedFeaturesexpanded.apply(lambda df: findrettime(df['Raw files'],\
                                                     df['Multiplet ids'],\
                                                     allpeptides),axis = 1)                          
    return matchedFeaturesexpanded

def reshape(allpeptides):
    files = allpeptides['Raw file'].unique()
    for x in range(len(files)):
        allpeptides['Raw file'] = allpeptides['Raw file'].replace({files[x]:x})
    allpeptides = allpeptides.set_index(['Raw file','Feature id'])
    return allpeptides

def pir2(df, c1, c2):
    colc1 = df[c1].str.split(';')
    colc2 = df[c2].str.split(';')
    clst1 = colc1.values.tolist()
    clst2 = colc2.values.tolist()

    lens = [len(l) for l in clst1]
    j = df.columns.get_loc(c1)
    v = df.values
    n, m = v.shape
    r = np.arange(n).repeat(lens)
    result = pd.DataFrame(np.column_stack([v[r, 0:j], np.concatenate(clst1), v[r, j+1:]]),
                 columns=df.columns)
    result[c2] = np.concatenate(clst2)
    return result

def plotscatterbox(dfs):
    x = dfs['Calibrated retention time']
    y = dfs['allPeptides retention time']
    plt.figure()
    plt.scatter(x,y)
    plt.figure()
    plt.boxplot((x-y).astype(float))
    
if __name__ == "__main__":  
    
    start_time = time.time();
    data = preprocessing()
    plotscatterbox(data)
    print("---Runtime =  %s seconds ---" % (time.time() - start_time))   
