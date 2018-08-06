# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import time as time 
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from perseuspy import pd
from tkinter import Tk
from tkinter import filedialog
from scipy.stats import kde
from pylab import polyfit, poly1d
import seaborn as sb
#    TEST: PARTITION Groups of Reports into Protein IDs  
#          --> Works, hard way tho A REMINDER FOR YOURSELF!!11!
# =============================================================================
# a = MS_report_open    
# arr = a.index
# results = []
# str = ''
# num = arr[0] 
# k = 0;
# for x in range(len(a)):
#     
#     str = a['Protein ID'][num] 
#     
#     if (k == (len(a)-1)):
#         results.append(str)            
#         break 
#     elif(arr[k]==arr[k+1]):
#         
#         num = arr[k]
#         str = a['Protein ID'][num] 
#         l = len(str)
#          
#         str2=str.values[0]    
#         d = 1
#         for z in range(l-1):
#             
#             str2 = str2 + ';' + str.values[d]      
#             d = d + 1
#             
#         results.append(str2)
#         num = arr[k+l]         
#         k = k + l   
#     else: 
#         k = k + 1
#         num = arr[k]
#         results.append(str)
# =============================================================================

# =============================================================================
#                                       FUNCTIONS:
# =============================================================================   
def compare(a,b):
    start_time = time.time()
    c = np.in1d(a,b)
    plt.figure()       
    venn2 (subsets = (len(a),len(b),c[c==True].size),set_labels=('MSFragger','MaxQuant'))
    print("---Runtime =  %s seconds ---" % (time.time() - start_time))         
 
def groupncompareprot(a,b):
    a = a.groupby(a.index)['SubGroup'].unique().str.join(';')
    b = b['Majority protein IDs']
    compare(a,b)
    
def groupncomparepep(a,b):
    
    a = a['Peptide']
    b = b['Sequence']
    compare(a.values,b.values)    
    
def extractMod(a):    
    assi_Mod = a['PSMs with Assigned Modifications'].values
    obs_Mod = a['PSMs with Observed Modifications'].values
    bin = a['Mass Bin'].values
    
    pos_assi =  np.where(assi_Mod>0)
    pos_obs = np.where(obs_Mod>0)
    num_assi = assi_Mod[pos_assi]
    num_obs = obs_Mod[pos_obs]
    bin_assi = bin[pos_assi]
    bin_obs = bin[pos_obs]    
  
 
    print(bin_obs,num_obs)
    print(bin_assi,num_assi)
    
# =============================================================================
# Test Area: PLOTs 
# =============================================================================
def plotsMS (ms_psmopen, range_int, zerovar, minbin, binsize, binrange, nomi ):   
    # Pie plot
    print('MSFragger Plot Parameters: min. modification: '+ str(range_int) + \
          '; removed '+str(zerovar)+ 'bins around zero'+'; min. binsize: '\
          +str(minbin)+'; binsize: '+str(binsize))
    plt.figure()
    dfMS = ms_psmopen.groupby('Observed Modifications')\
    .size().reset_index(name='No: of Modifications')
    
    range_intrest = range_int
    df2 = dfMS[dfMS['No: of Modifications']>range_intrest]
    df3 = dfMS[dfMS['No: of Modifications']<range_intrest]
    
    others = 0
    
    df2 = df2.reset_index(drop=True)
    df3 = df3.reset_index(drop=True)
    for x in range(len(df3)):
        others = others + df3['No: of Modifications'][x]
        
    df2.loc[len(df2)] = ['Others',others]
    df2.set_index('Observed Modifications', inplace=True)
    df2.drop('Unknown', inplace=True)
    df2 = df2['No: of Modifications'].sort_values()
    df2.plot.pie(y = 'No: of Modifications', legend=False).set_ylabel('')
    
    
    # Create new Bar Plot for Mass Bins ??not accurate == WHY??
    plt.figure()
    stepsize = binsize
    # set x bin range here
    bins = np.arange(-binrange,binrange,stepsize)
    var = round(len(bins)/2)
    
    if (nomi == True):
        sorted = round(ms_psmopen.sort_values('Adjusted Delta Mass'))
    else: 
        sorted = ms_psmopen.sort_values('Adjusted Delta Mass')

    hist,bin_edges = np.histogram(sorted['Adjusted Delta Mass'].values, bins = bins)
    hist[(var)] = 0 # remove unknown modifications around 0 
                                # adjusted Delta mass 
                  # --> prob. no modification
    #hist[hist>0] = np.log2(hist[hist>0])            
    plt.xlabel('Mass Bins', fontsize=15)  
    plt.ylabel('log-scale quantity', fontsize=13)                
             
    plt.bar(bin_edges[:-1], hist,log = True, width = 1)
#    plt.ylim([0,np.log(10000000000000)]) 
    plt.show()   
    
    
    # Pie PLot of Mass Bins
  
    # do only once per program execution 
    hist = np.append(hist,0)
    
    plt.figure()
    histMS = {'Histodata':hist, 'Bins':bins }
    dMS = pd.DataFrame(histMS)
    d2 = dMS[dMS['Histodata'] > minbin] 
    
    round_bin = np.around(d2['Bins'], decimals = 2)
    round_bin = round_bin.reset_index(drop = True)
    
    for x in range(len(round_bin)):
        round_bin[x] = str(round_bin[x])\
        + ' to ' + str(round((round_bin[x] + stepsize),2)) + ' bin'
        
    d2.set_index(round_bin,inplace = True)
    d2 = d2['Histodata'].sort_values()
    d2.plot.pie(y = 'Histodata', legend=False).set_ylabel('')
    return dMS
    
# =============================================================================
# MQ Plots
# =============================================================================
def plotsMQ (mq_dep, range_int, zerovar, minbin, binsize,binrange ):
    # Pie Plot of MQ found Modifications
    print('MaxQuant Plot Parameters: min. modification: '+ str(range_int) + \
          '; removed '+str(zerovar)+ 'bins around zero'+'; min. binsize : '\
          +str(minbin)+'; binsize: '+str(binsize))
    plt.figure()
    dfMQ = mq_dep.groupby('DP Modification')\
    .size().reset_index(name='No: of Modifications')
    
    #make range of intrest relative to datasize?
    range_intrest = range_int
    df2 = dfMQ[dfMQ['No: of Modifications']>range_intrest]
    df3 = dfMQ[dfMQ['No: of Modifications']<range_intrest]
    
    others = 0
    df2 = df2.reset_index(drop=True)
    df3 = df3.reset_index(drop=True)
    for x in range(len(df3)):
        others = others + df3['No: of Modifications'][x]
        
    df2.loc[len(df2)] = ['Others',others]
    df2.set_index('DP Modification', inplace=True)
    #df2.drop('Unknown', inplace=True) #dont kick them out
    df2 = df2['No: of Modifications'].sort_values()
    df2.plot.pie(y = 'No: of Modifications', legend=False).set_ylabel('')
    

    # Bar PLot for MQ Mass Bins 
    plt.figure()
    stepsize = binsize
    #adjust binrange here
    bins = np.arange(-binrange,binrange,stepsize)
    var = round(len(bins)/2)
    sorted = round(mq_dep.sort_values('DP Cluster Mass'))
    sorted['DP Cluster Mass'] = round(sorted['DP Cluster Mass'].astype(float))
    sorted = sorted['DP Cluster Mass'].values
    hist,bin_edges = np.histogram(sorted, bins = bins)
    hist[var] = 0 
                                # remove unknown modifications around 0 
                                # adjusted Delta mass 
                                # --> prob. no modification
    #hist[hist>0] = np.log2(hist[hist>0])                
    plt.xlabel('Mass Bins', fontsize=15)  
    plt.ylabel('log-scale quantity', fontsize=13)                   
    plt.bar(bin_edges[:-1], hist, width = 1,log = True)
#    plt.ylim([0,np.]) 

    plt.show()   
    
    
    # Pie Plot for Bins
    hist = np.append(hist,0)
    
    plt.figure()
    histMQ = {'Histodata':hist, 'Bins':bins }
    dMQ = pd.DataFrame(histMQ)
    d2 = dMQ[dMQ['Histodata'] > minbin] 
    
    round_bin = np.around(d2['Bins'], decimals = 2)
    round_bin = round_bin.reset_index(drop = True)
    
    for x in range(len(round_bin)):
        round_bin[x] = str(round_bin[x])\
        + ' to ' + str(round((round_bin[x] + stepsize),2)) + ' bin'
        
    d2.set_index(round_bin,inplace = True)
    d2 = d2['Histodata'].sort_values()
    d2.plot.pie(y = 'Histodata', legend=False).set_ylabel('')
    
    return dMQ

def idscatter(a,b,counts):

    df = a 
    df['Histodata MS'] = b['Histodata']
    sort = df[(df['Histodata'] > counts) & (df['Histodata MS'] > counts)]
    x = sort['Histodata'].reset_index(drop = True)
    y = sort['Histodata MS'].reset_index(drop = True)
    n = sort['Bins']
  
    plt.figure()
    ax = plt.subplot()
    plt.scatter(x,y)
    plt.xscale('log')
    plt.yscale('log')
    fit = np.polyfit(np.log(x), np.log(y), deg=1)
    plt.plot(x, fit[0] * x + fit[1], color='red')
    plt.xlabel("MaxQuant mass bins", fontsize = 20)
    plt.ylabel("MSFragger mass bins", fontsize = 20)
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    
    for i, txt in enumerate(n):
        ax.annotate(txt, (x[i],y[i]),fontsize = 20) 

       
def addRawid(MQ):
    # add the id of MQ Raw files (only for second peptides)
    # compare MSFragger to MQ       
# =============================================================================
#     MQ[MQ['mod. Raw File'] == 'b1928_293T_proteinID_08A_QE3_122212.35152']
#     ids = v['id'].values
#     MQ['mod. Raw File'].iloc[ids] = MQ['mod. Raw File'].iloc[ids] +'.'+ str(1234567890)
# =============================================================================
    
    mult = MQ['mod. Raw File'].value_counts()>1   
    mult = mult[mult == True]
    for x in range(len(mult)):
        rawidex = mult.index[x]
        dupl = MQ[MQ['mod. Raw File'] == rawidex]['Type'] == 'MULTI-SECPEP'
        duplidex = dupl[dupl == True].index
# =============================================================================
#         MQ['mod. Raw File'].iloc[duplidex.values] = \
#         MQ['mod. Raw File'].iloc[duplidex.values] + '.' + str(123456789)
# =============================================================================
        for y in range(len(duplidex)):
            idex = duplidex[y]
            MQ['mod. Raw File'][idex] = MQ['mod. Raw File'][idex] + '.' + \
            str(MQ['id'][idex])
            
def scoreplt(a,b):
    plt.figure()
    a = a.sort_values('ScanID')
    b = b.sort_values('mod. Raw File')
    
    x = a['Hyperscore']
    x.plot.density(legend=True,label ='MSFragger').set_xlim(0,(x.max()+1))
    
    c = np.in1d(a['ScanID'],b['mod. Raw File'])
    x = a[c]['Hyperscore']
    x.plot.density(secondary_y = True,legend=True,label='MSF and MQ identified',title = 'MSFragger Scores', style = 'y--').set_xlim(0,(x.max()+1))
    
    
    plt.figure()
    y = b['Score']
    y.plot.density(legend=True,label ='MaxQuant').set_xlim(0,(y.max()+1))
    c = np.in1d(b['mod. Raw File'],a['ScanID'])
    y = b[c]['Score']
    
    y.plot.density(secondary_y = True,legend=True,label='MSF and MQ identified',title = 'MaxQuant Scores', style = 'y--').set_xlim(0,(y.max()+1))
     
    plt.figure() 
    fit = np.polyfit(x,y,deg=1)
    plt.plot(x, fit[0] * x + fit[1], color='red')
    # Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
    nbins=150
    k = kde.gaussian_kde([x,y])
    xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))
     
    # Make the plot
    plt.pcolormesh(xi, yi, zi.reshape(xi.shape))
    
    plt.xlabel("MaxQuant")
    plt.ylabel("MSFragger")
    plt.colorbar()
    plt.show()   
     
    plt.figure()
    fit = np.polyfit(x,y,deg=1)
    plt.plot(x, fit[0] * x + fit[1], color='red')
    plt.scatter(x,y, alpha=0.1)
    plt.xlabel("MaxQuant")
    plt.ylabel("MSFragger")
    plt.show()   
    
def compareScanID(a,b):
    MS_Scans = []
    for x in range(len(a)):
        the,scanid,to,do = a['Spectrum'][x].split('.')      
        scanid = int(scanid)
        scanid = str(scanid)
        MS_Scans.append(the+'.'+scanid)    
    if 'ScanID' in a:
        scoreplt(a,b)
        compare(a['ScanID'],b['mod. Raw File'])
        
    else:          
        a.insert(loc = 0,column='ScanID',value=MS_Scans)
        arr = b['Raw file'] +'.'+ b['Scan number'].astype(str)
        b.insert(loc = 0,column='mod. Raw File',value=arr)
        start_time2 = time.time();
        addRawid(b) # dont really need it... yes you do meh
        print("---Runtime =  %s seconds ---" % (time.time() - start_time2))   
        compareScanID(a,b)
 
def testMQ():
    root = Tk()
    root.withdraw()
    a = filedialog.askopenfilename(initialdir = "D:\Thomas\PhytonSCripts\MQOutput" \
           ,title = "Choose a MQ nomatch.deppep file to plot",\
           filetypes = (("deppep files","*.deppep"),("all files","*.*")))
    a = pd.read_table(a,low_memory=False)
    a = a.drop(0).reset_index(drop = True)
    
    b = filedialog.askopenfilename(initialdir = "D:\Thomas\PhytonSCripts\MQOutput" \
           ,title = "Choose a MQ matching.deppep file to plot",\
           filetypes = (("deppep files","*.deppep"),("all files","*.*")))
    b = pd.read_table(b,low_memory=False)
    b = b.drop(0).reset_index(drop = True)
    
    arg=a[['Raw file','DP Base Raw File','DP Proteins','DP Base Sequence']] 
    argb=b[['Raw file','DP Base Raw File','DP Proteins','DP Base Sequence']]

    uniq = arg['DP Base Sequence'].unique()
    uniqb = argb['DP Base Sequence'].unique()
    out = uniq[np.in1d(uniq,uniqb)]
    df = pd.DataFrame({'Raw Files':arg['Raw file'].unique()})
    df['counta'] = 0
    df['countb'] = 0
    count = 0
    countb = 0
    
    for x in range(len(out)):
       coin = arg[arg['DP Base Sequence'] == out[x]]
       coinb = argb[argb['DP Base Sequence'] == out[x]]    
       coincount = np.in1d(coin['Raw file'].unique(),coinb['Raw file'].unique())
       coinbcount = np.in1d(coinb['Raw file'].unique(),coin['Raw file'].unique())
       count = count + coincount[coincount==False].size
       countb = countb + coinbcount[coinbcount==False].size
       
       for y in (coin['Raw file'].unique()[coincount==False]):
           df['counta'][df['Raw Files'] == y] = df['counta'][df['Raw Files'] == y] + 1
       for z in (coinb['Raw file'].unique()[coinbcount==False]):
           df['countb'][df['Raw Files'] == z] = df['countb'][df['Raw Files'] == z] + 1
    df.plot.bar()
           
    print('Amount of peptide sequences found in the Raw Files of file 1 but not'\
          ' found in the Raw Files of file 2 is ',count, '. For file 2: ', countb)  
    print('!!!Both peptide sequences need to be present in both files!!!')
    print('Amount of dp. Peptides in file 1: ',len(a), ' in file 2: ', len(b))
    
    values = pd.DataFrame(a['Raw file'].value_counts())
    values['Raw'] = b['Raw file'].value_counts().values
    values.columns = ['No match', 'Matching']
    values = values.sort_index()
    values.plot.bar()
    
    a['Intensity'] = a['Intensity'].astype(float)
    b['Intensity'] = b['Intensity'].astype(float)

    raw_tablea = pd.pivot_table(a, values = 'Intensity',\
    index = ['DP Cluster Index','DP AA','DP Base Sequence','DP Probabilities'], columns = 'Raw file')
    raw_tableb = pd.pivot_table(b, values = 'Intensity',\
    index = ['DP Cluster Index','DP AA','DP Base Sequence','DP Probabilities'], columns = 'Raw file')

    raw_tablea = raw_tablea.reset_index()
    raw_tableb = raw_tableb.reset_index()

    return raw_tablea,raw_tableb


# =============================================================================
# MAIN CONTROLLER
# =============================================================================

def mainpep():
    root = Tk()
    root.withdraw()
    a = filedialog.askopenfilename(initialdir = "D:\Thomas\PhytonSCripts\MSOutput" \
           ,title = "Choose a MSFragger peptide file",\
           filetypes = (("tsv files","*.tsv"),("all files","*.*")))
    b = filedialog.askopenfilename(initialdir = "D:\Thomas\PhytonSCripts\MQOutput" \
           ,title = "Choose a MQ peptide file",\
           filetypes = (("txt files","*.txt"),("all files","*.*")))
    a = pd.read_csv(a,sep = '\t', header = 0)
    b = pd.read_table(b)
    groupncomparepep(a,b)
    
def mainprot():
    root = Tk()
    root.withdraw()
    a = filedialog.askopenfilename(initialdir = "D:\Thomas\PhytonSCripts\MSOutput" \
           ,title = "Choose a MSFragger report file",\
           filetypes = (("tsv files","*.tsv"),("all files","*.*")))
    b = filedialog.askopenfilename(initialdir = "D:\Thomas\PhytonSCripts\MQOutput" \
           ,title = "Choose a MQ proteinGroups.txt file",\
           filetypes = (("txt files","*.txt"),("all files","*.*")))
    a = pd.read_csv(a,sep = '\t', header = 0)
    b = pd.read_table(b)
    groupncompareprot(a,b)
    
def mainmsms():
    root = Tk()
    root.withdraw()
    a = filedialog.askopenfilename(initialdir = "D:\Thomas\PhytonSCripts\MSOutput" \
           ,title = "Choose a MSFragger psm.tsv file",\
           filetypes = (("tsv files","*.tsv"),("all files","*.*")))
    b = filedialog.askopenfilename(initialdir = "D:\Thomas\PhytonSCripts\MQOutput" \
           ,title = "Choose a MQ msms.txt file",\
           filetypes = (("txt files","*.txt"),("all files","*.*")))
    a = pd.read_csv(a,sep = '\t', header = 0)
    b = pd.read_table(b)
    print('Comparing big psm/msms files can take a while.')
    compareScanID(a,b)

def maindeppep():
    root = Tk()
    root.withdraw()
    a = filedialog.askopenfilename(initialdir = "D:\Thomas\PhytonSCripts\MSOutput" \
           ,title = "Choose a MSFragger psm.tsv file to plot",\
           filetypes = (("tsv files","*.tsv"),("all files","*.*")))
    b = filedialog.askopenfilename(initialdir = "D:\Thomas\PhytonSCripts\MQOutput" \
           ,title = "Choose a MQ DP_peptides.deppep file to plot",\
           filetypes = (("deppep files","*.deppep"),("all files","*.*")))
    
    a = pd.read_csv(a,sep = '\t', header = 0)
    b = pd.read_table(b,low_memory=False)
    b = b.drop(0).reset_index(drop = True)
    
# =============================================================================
#     Tune Plots here: e.g: 
#     file, reduce to 1000 most common modificaions , remove 10 bins around 0,
#     reduce to 1000 most common bins, configure binsteps and size 
# =============================================================================

    histMS = plotsMS(a, 1500, 5, 1500, 1,500.5,False)
    histMQ = plotsMQ(b, 1500, 0, 1500, 1,500.5)
   
    idscatter(histMQ,histMS,110)
    return histMQ,histMS
    
def main():
    
    mainmsms()
    mainprot()
    mainpep()
    maindeppep()

    
    
    