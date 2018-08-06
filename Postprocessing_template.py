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
import matplotlib.pyplot as plt
from janspivot import janspivot, openfile
# table = Parallel(n_jobs = -2)(extractdupls(dupl,x,table) for x in range(len(dupl)))

def recalculateAAProbabilities(prob_threshold,showplots,clusters):
    a = openfile()
    a['DP Positional Probability'] = a['DP Positional Probability'].astype(float)
    filtered = a[a['DP Positional Probability']>prob_threshold]
    x = filtered.groupby(['DP Cluster Index'])['DP Cluster Index'].count()
    x = x[x>1]
    regression = a
    regression['DP Probabilities after Regression'] = a['DP Probabilities']
    for cluster_at in range(len(x)):
        cluster = x.sort_values(ascending = False).index[cluster_at]
       
        #count the aminoacids per cluster where pos. Prob > 0.9 
        maxcluster = filtered.loc[lambda df: df['DP Cluster Index'] == cluster, :]
        maxcluster = maxcluster.reset_index()
        dpAAs = countAA(maxcluster)
        #adjoust Probabilites for the whole Cluster
        whole_cluster = a.loc[lambda df: df['DP Cluster Index'] == cluster, :]
        meanProb = maxcluster['DP Positional Probability'].mean()        
        whole_cluster['DP Probabilities after Regression'] = \
        whole_cluster['DP Probabilities'].apply(lambda df: doStatistics(df,meanProb,dpAAs))    
        regression = fillregressiontable(whole_cluster,regression)            
        # show counts of the  biggest clusters
        if showplots == True:
            if cluster_at < clusters:
                normalized = dpAAs
                normalized['Count'] = (normalized['Count']/sum(normalized['Count']))*100
                ax = normalized.plot.bar(x = 0,y = 'Count',rot=0,fontsize = 20)
                ax.legend().set_visible(False)
                ax.set_xlabel('Amino acid',fontsize = 20)
                ax.set_ylabel('Percentage share',fontsize = 20)
                plt.title('Cluster: ' + maxcluster['DP Modification'].values[0],fontsize = 20)
                         
    return regression

def countAA(maxcluster):
    dpAAs = pd.DataFrame(['nterm','cterm','G','A','V','L','I','P','F','Y','W','S','T',\
                           'C','M','N','Q','K','R','H','D','E'])
    dpAAs['Count'] = 0
    for x in range(len(maxcluster)-1):
        x = x + 1 
        i = dpAAs['Count'][dpAAs[0].isin(list(maxcluster['DP AA'][x]\
              .split(';')))].index
        dpAAs['Count'].loc[i] = dpAAs['Count'].loc[i] + 1
    return dpAAs
    
def doStatistics(whole_cluster_reg_row,meanProb,dpAAs):
    dfinfo = findallin([],[],[],whole_cluster_reg_row,'(',')')
    dfinfo = dfinfo.sort_values(by = 1)    
    
    beg = whole_cluster_reg_row.find('(') 
    end = whole_cluster_reg_row.find(')')  
    resArr = []
    # Formula atm: p  = p + (countAA * meanProb /max(countAA))
    #              P  = p/∑p's            
    for x in range(len(dfinfo)): 
        p = float(dfinfo[0][x]) # p = positional Prob of AA before Reg.
        countAA = float(dpAAs['Count'][dpAAs[0] == dfinfo[2][x]]) # count of specific AA 
        weight = 3
        if x == 0 and float(dpAAs['Count'][dpAAs[0] == 'nterm' ]) > 0:
            # do statistics here
            countAAnterm = countAA + float(dpAAs['Count'][dpAAs[0] == 'nterm'])
            res = weight * p + (countAAnterm * meanProb / dpAAs['Count'].max())
            res = round(res,3)
            resArr.append(res)
        elif x == len(dfinfo) and float(dpAAs['Count'][dpAAs[0] == 'cterm']) > 0:        
            # do statistics here
            countAActerm = countAA + float(dpAAs['Count'][dpAAs[0] == 'cterm'])
            res = weight * p + (countAActerm * meanProb / dpAAs['Count'].max())
            res = round(res,3)
            resArr.append(res)
        else:
            res = weight * p + (countAA * meanProb / dpAAs['Count'].max())
            res = round(res,3)
            resArr.append(res)         
    # caluculate sum for all Results: P  = p/∑p's 
    resSum = sum(resArr)
    for x in range(len(dfinfo)): 
        res = str(round(resArr[x]/resSum,3))
        # insert results here
        whole_cluster_reg_row = insert_str(whole_cluster_reg_row,res,beg+1,end)        
        beg = whole_cluster_reg_row.find('(',end) 
        end = whole_cluster_reg_row.find(')',beg) 
        
    return whole_cluster_reg_row  

def findallin (cache,position,AA,thestr,char1,char2):
#    initialize cache,postiton,AA as []        
  
    if(thestr.find(char1) > -1 and thestr.find(char2) > -1):    
        substr = thestr[thestr.find(char1)+1:thestr.find(char2)]
        subsuperstr = thestr[thestr.find(char2)+1:len(thestr)]
        cache.append(substr)
        
        if not position: 
            position.append(thestr.find(char1)-1)
            AA.append(thestr[thestr.find(char1)-1])
        else:
            idex = position[len(position)-1] + thestr.find(char1)
            position.append(idex)
            AA.append(thestr[thestr.find(char1)-1])

        findallin(cache,position,AA,subsuperstr,char1,char2)
        
    df = pd.DataFrame(cache)
    df.insert(1,1,position)
    df.insert(2,2,AA)
    return df

# insert a string into another one form x to y
def insert_str(string, str_to_insert, beg,end):
    return string[:beg] + str_to_insert + string[end:]      

# function to fill the regression table
def fillregressiontable(whole_cluster,regression):
    whole_cluster = whole_cluster.reset_index()
    regression['DP Probabilities after Regression'][whole_cluster['index']] = \
    whole_cluster['DP Probabilities after Regression']
    return regression

# =============================================================================
# Start a custom pivot of a MaxQuant .deppep file here
# =============================================================================    
if __name__ == "__main__":  
    
    start_time = time.time();
    regcalc = recalculateAAProbabilities(0.9,True,5)
    pivot = janspivot(regcalc,True)
    print("---Runtime =  %s seconds ---" % (time.time() - start_time))   
