# -*- coding: utf-8 -*-
"""
Created on Wed May 23 15:53:21 2018

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
    # to do: include nterm, cterm here:
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
    
    dfinfo[0] = dfinfo[0].astype(float)
    dfinfo['Ratio p/Aas'] = 1.0 # ratio is one except multiple p/Aa
    dpAAs = dpAAs.set_index(0)
    # calculate ratio of probabilities per Aminoacid
    multAAs = dfinfo[2].value_counts()
    dfinfo = dfinfo.set_index(2) # set index to Aas for better processing
    multAAs = multAAs[multAAs>1].index
    for aa in multAAs:
        AAs = dfinfo[0][aa]
        dfinfo['Ratio p/Aas'][aa] = AAs[aa]/sum(AAs[aa])
   
    # initialise counts
    dfinfo['posterior'] = 0
    s = whole_cluster_reg_row # parameter for insertion of prob in string
    end = len(dfinfo)-1 # parameter for insertion of prob in string
    
    probmean = pd.DataFrame(dfinfo.groupby(dfinfo.index)[0].sum())
    probmean['Data'] = 0
    for aa in probmean.index: 
        probmean['Data'][aa] = dpAAs['Count'][aa]

    # Dirichlet 1: without consideration of n/cterm
    # Formula atm: posterior = prior(1 + ∑Counts>0.9) + Likelihood(Counts>0.9 specifc Aa's)  
    # p = posterior_k / ∑posterior                   
    weight = sum(dpAAs['Count'].iloc[2:]) # weigth = sum of all counts (p>0.9)
    prior = 1 + weight*probmean[0] # Dirichlet prior distribution 
    posterior = prior + probmean['Data'] # post = alpha + multinomial
   
    for aa in dfinfo.index: 
        # caclculate real posterior out of ratio
        dfinfo['posterior'][aa] = dfinfo['Ratio p/Aas'][aa] * posterior[aa] 
        
    dfinfo['posterior p'] = dfinfo['posterior']/sum(dfinfo['posterior'])

    # Dirichlet 2 with n/cterm:
    dfinfo = dfinfo.reset_index()
    for idex, row in dfinfo.iterrows():
        if row[1] == 0: # if n-terminal 
            dfinfo['posterior'][idex] = dfinfo['posterior'][0] + dpAAs['Count']['nterm']
        elif s[len(s)-1] == ')' and row[1] == dfinfo[1].iloc[end]: # if c-terminal
            dfinfo['posterior'][idex] = dfinfo['posterior'][end] + dpAAs['Count']['cterm']
        else:   # if in the middle
            dfinfo['posterior'][idex] = dfinfo['posterior'][idex] \
            + weight -dpAAs['Count']['nterm']-dpAAs['Count']['cterm']
   
    dfinfo['posterior p'] = dfinfo['posterior']/sum(dfinfo['posterior'])        

    # here the probabilities are modified 
    beg = whole_cluster_reg_row.find('(') 
    end = whole_cluster_reg_row.find(')')  
    for idex in dfinfo.index: 
        res = str(dfinfo['posterior p'][idex])
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

# insert a string into another one from  to y
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
