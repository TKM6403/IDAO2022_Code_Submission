#!/usr/bin/env python
# coding: utf-8

# # Read in Public Data
# This portion of the code reads the public data files .

# In[ ]:


# Read in Public Data

import pandas as pd 
import numpy as np

#read in data
targets=pd.read_csv("C:/Users/aerob/Downloads/IDAO 2022 (1)/IDAO 2022/targets.csv")
targets._id=(targets._id).astype(str)

targets = targets.rename(columns={'_id': 'X_id'})

import os
os.chdir("C:/Users/aerob/Downloads/IDAO 2022 (1)/IDAO 2022/dichalcogenides_public/structures/")

import json
 
# Opening JSON file
f = open(targets.X_id[0]+".json")
 
# returns JSON object as
# a dictionary
data = json.load(f)


targets['a']=0
targets['b']=0
targets['c']=0
targets['alpha']=0
targets['beta']=0
targets['gamma']=0
targets['volume']=0


for id in range(0,len(targets['a'])):
    # Opening JSON file
    f = open(targets.X_id[id]+".json")

    # returns JSON object as
    # a dictionary
    samplist = json.load(f)
    targets['a'][id]=samplist['lattice']['a']
    targets['b'][id]=samplist['lattice']['b']
    targets['c'][id]=samplist['lattice']['c']
    targets['alpha'][id]=samplist['lattice']['alpha']
    targets['beta'][id]=samplist['lattice']['beta']
    targets['gamma'][id]=samplist['lattice']['gamma']
    targets['volume'][id]=samplist['lattice']['volume']


def listoflabels(id):
    listinit=[]
    f = open(targets.X_id[id]+".json")

    # returns JSON object as
    # a dictionary
    samplist = json.load(f)
    for i in range(0, len(samplist['sites'])):
        listinit.append(samplist['sites'][i]['label'])
    
    return pd.Series(listinit)


#substitutional impurities
targets['num_substitution']=0
for id in range(0,len(targets['num_substitution'])):
    f = open(targets.X_id[id]+".json")

    # returns JSON object as
    # a dictionary
    samplist = json.load(f)
    y=(listoflabels(id)).value_counts()
    if(len(y.index) > 2):
        targets['num_substitution'][id]=sum(y[2:])
    



#type of substitution
targets['sub1_location_x']=0
targets['sub1_location_y']=0
targets['sub1_location_z']=0
targets['sub1_element']=0

targets['sub2_location_x']=0
targets['sub2_location_y']=0
targets['sub2_location_z']=0
targets['sub2_element']= 0

targets['sub3_location_x']=0
targets['sub3_location_y']=0
targets['sub3_location_z']=0
targets['sub3_element']= 0
for id in range(0, len(targets['sub1_element'])) :
    f = open(targets.X_id[id]+".json")

    # returns JSON object as
    # a dictionary
    samplist = json.load(f)  
    x=listoflabels(id).value_counts()


    if targets['num_substitution'][id]==1 :
        targets.loc[id,'sub1_location_x']=samplist['sites'][np.where(listoflabels(id)==listoflabels(id).value_counts().index[2])[0][0]]['abc'][0]
        targets.loc[id,'sub1_location_y']=samplist['sites'][np.where(listoflabels(id)==listoflabels(id).value_counts().index[2])[0][0]]['abc'][1]
        targets.loc[id,'sub1_location_z']=samplist['sites'][np.where(listoflabels(id)==listoflabels(id).value_counts().index[2])[0][0]]['abc'][2]
        targets.loc[id, 'sub1_element']=listoflabels(id).value_counts().index[2]

    if targets['num_substitution'][id]==2:
        targets.loc[id,'sub1_location_x']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[0]]['abc'][0]
        targets.loc[id,'sub1_location_y']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[0]]['abc'][1]
        targets.loc[id,'sub1_location_z']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[0]]['abc'][2]
     
        targets.loc[id, 'sub1_element']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[0]]['label']

        targets.loc[id,'sub2_location_x']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[1]]['abc'][0]
        targets.loc[id,'sub2_location_y']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[1]]['abc'][1]
        targets.loc[id,'sub2_location_z']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[1]]['abc'][2]
        targets.loc[id, 'sub2_element']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[1]]['label']

  
    if targets['num_substitution'][id]==3:
        targets.loc[id,'sub1_location_x']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[0]]['abc'][0]
        targets.loc[id,'sub1_location_y']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[0]]['abc'][1]
        targets.loc[id,'sub1_location_z']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[0]]['abc'][2]
     
        targets.loc[id, 'sub1_element']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[0]]['label']

        targets.loc[id,'sub2_location_x']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[1]]['abc'][0]
        targets.loc[id,'sub2_location_y']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[1]]['abc'][1]
        targets.loc[id,'sub2_location_z']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[1]]['abc'][2]

        targets.loc[id, 'sub2_element']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[1]]['label']


        targets.loc[id,'sub3_location_x']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[2]]['abc'][0]
        targets.loc[id,'sub3_location_y']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[2]]['abc'][1]
        targets.loc[id,'sub3_location_z']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[2]]['abc'][2]

        targets.loc[id, 'sub3_element']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[2]]['label']







targets['num_vacancies']=0
total_points=192
for id in range(0, len(targets['num_vacancies'])):
    f = open(targets.X_id[id]+".json")

    # returns JSON object as
    # a dictionary
    samplist = json.load(f)    
    targets['num_vacancies'][id]=total_points-len(samplist['sites'])


def listofcoords(id):
    listinit=[]
    f = open(targets.X_id[id]+".json")

    # returns JSON object as
    # a dictionary
    samplist = json.load(f)
    for i in range(0, len(samplist['sites'])):
        listinit.append(samplist['sites'][i]['abc'])
    
    return pd.Series(listinit)

listofcoords(3)

targets['mo_vacancies']=0
targets['s_vacancies']=0
targets['s_bottom_vacancies']=0
targets['s_top_vacancies']=0

targets['vac1_element']=0
targets['vacancy1_location_x']=0
targets['vacancy1_location_y']=0
targets['vacancy1_location_z']=0

targets['vac2_element']=0
targets['vacancy2_location_x']=0
targets['vacancy2_location_y']=0
targets['vacancy2_location_z']=0

targets['vac3_element']=0
targets['vacancy3_location_x']=0
targets['vacancy3_location_y']=0
targets['vacancy3_location_z']=0

good_lattice=listofcoords(3)
good_df=pd.DataFrame(list(good_lattice), columns = ['x', 'y','z'])


for id in list(np.where(targets['num_vacancies'] > 0)[0]):
    f = open(targets.X_id[id]+".json")

    # returns JSON object as
    # a dictionary
    samplist = json.load(f)
    bad_lattice=listofcoords(id)
    bad_df=pd.DataFrame(list(bad_lattice), columns = ['x', 'y','z'])
    difference_df=pd.concat([good_df,bad_df]).drop_duplicates(keep=False)
    difference_df=pd.DataFrame(difference_df).reset_index()
    
#GET SET DIFFERENCES CORRECT
    for i in range(0,targets['num_vacancies'][id]):
        if difference_df['z'][i] < 0.2:
            targets['s_vacancies'][id]=targets['s_vacancies'][id]+1
            targets['s_bottom_vacancies'][id]=targets['s_bottom_vacancies'][id]+1
        
        if difference_df['z'][i] > 0.3:
            targets['s_vacancies'][id]=targets['s_vacancies'][id]+1
            targets['s_top_vacancies'][id]=targets['s_bottom_vacancies'][id]+1
        
        if difference_df['z'][i] > 0.2 and difference_df['z'][i] < 0.3:
            targets['mo_vacancies'][id]=targets['mo_vacancies'][id]+1
        
        
  
  
    if targets['num_vacancies'][id]==1:
        targets.loc[id,'vacancy1_location_x']=difference_df['x'][0]
        targets.loc[id,'vacancy1_location_y']=difference_df['y'][0]
        targets.loc[id,'vacancy1_location_z']=difference_df['z'][0]

      

    if targets['num_vacancies'][id]==2:
        targets.loc[id,'vacancy1_location_x']=difference_df['x'][0]
        targets.loc[id,'vacancy1_location_y']=difference_df['y'][0]
        targets.loc[id,'vacancy1_location_z']=difference_df['z'][0]
        targets.loc[id,'vacancy2_location_x']=difference_df['x'][1]
        targets.loc[id,'vacancy2_location_y']=difference_df['y'][1]
        targets.loc[id,'vacancy2_location_z']=difference_df['z'][1]
      
  
    if targets['num_vacancies'][id]==3:
        targets.loc[id,'vacancy1_location_x']=difference_df['x'][0]
        targets.loc[id,'vacancy1_location_y']=difference_df['y'][0]
        targets.loc[id,'vacancy1_location_z']=difference_df['z'][0]
        targets.loc[id,'vacancy2_location_x']=difference_df['x'][1]
        targets.loc[id,'vacancy2_location_y']=difference_df['y'][1]
        targets.loc[id,'vacancy2_location_z']=difference_df['z'][1]
        targets.loc[id,'vacancy3_location_x']=difference_df['x'][2]
        targets.loc[id,'vacancy3_location_y']=difference_df['y'][2]
        targets.loc[id,'vacancy3_location_z']=difference_df['z'][2]
        



for id in range(0,len(targets['vacancy1_location_z'])):
    if targets['vacancy1_location_z'][id] > 0.1 and targets['vacancy1_location_z'][id] < 0.2:
        targets.loc[id, 'vac1_element'] = "S"
    elif targets['vacancy1_location_z'][id] > 0.2 and targets['vacancy1_location_z'][id] < 0.3:
        targets.loc[id, 'vac1_element'] = "Mo"
    elif targets['vacancy1_location_z'][id] > 0.3 and targets['vacancy1_location_z'][id] < 0.4:
        targets.loc[id, 'vac1_element'] = "S"
    else:
        targets.loc[id, 'vac1_element'] = "None"

    if targets['vacancy2_location_z'][id] > 0.1 and targets['vacancy2_location_z'][id] < 0.2:
        targets.loc[id, 'vac2_element'] = "S"
    elif targets['vacancy2_location_z'][id] > 0.2 and targets['vacancy2_location_z'][id] < 0.3:
        targets.loc[id, 'vac2_element'] ="Mo"
    elif targets['vacancy2_location_z'][id] > 0.3 and targets['vacancy2_location_z'][id] < 0.4:
        targets.loc[id, 'vac2_element'] = "S"
    else:
        targets.loc[id, 'vac2_element'] = "None"
  
    if targets['vacancy3_location_z'][id] > 0.1 and targets['vacancy3_location_z'][id] < 0.2:
        targets.loc[id, 'vac3_element'] = "S"
    elif targets['vacancy3_location_z'][id] > 0.2 and targets['vacancy3_location_z'][id] < 0.3:
        targets.loc[id, 'vac3_element'] = "Mo"
    elif targets['vacancy3_location_z'][id] > 0.3 and targets['vacancy3_location_z'][id] < 0.4:
        targets.loc[id, 'vac3_element'] = "S"
    else:
        targets.loc[id, 'vac3_element'] = "None"

  
  



##new coordinate system
import math
targets['sub1_distance']=0
targets['sub2_distance']=0

targets['vac1_distance']=0
targets['vac2_distance']=0
targets['vac3_distance']=0

for i in range(0,len(targets['sub1_distance'])):
    targets.loc[i,'sub1_distance']= math.sqrt(targets['sub1_location_x'][i]**2+ targets['sub1_location_y'][i]**2+targets['sub1_location_z'][i]**2)

    targets.loc[i,'sub2_distance']= math.sqrt(targets['sub2_location_x'][i]**2+targets['sub2_location_y'][i]**2+targets['sub2_location_z'][i]**2)

    targets.loc[i,'sub3_distance']= math.sqrt(targets['sub3_location_x'][i]**2+targets['sub3_location_y'][i]**2+targets['sub3_location_z'][i]**2)

    targets.loc[i,'vac1_distance']= math.sqrt(targets['vacancy1_location_x'][i]**2+targets['vacancy1_location_y'][i]**2+targets['vacancy1_location_z'][i]**2)

    targets.loc[i,'vac2_distance']= math.sqrt(targets['vacancy2_location_x'][i]**2+targets['vacancy2_location_y'][i]**2+targets['vacancy2_location_z'][i]**2)
    targets.loc[i,'vac3_distance']= math.sqrt(targets['vacancy3_location_x'][i]**2+targets['vacancy3_location_y'][i]**2+targets['vacancy3_location_z'][i]**2)



##new coordinate system
import math
targets['sub1_distance2']=0
targets['sub2_distance2']=0

targets['vac1_distance2']=0
targets['vac2_distance2']=0
targets['vac3_distance2']=0

for i in range(0,len(targets['sub1_distance2'])):
    targets.loc[i,'sub1_distance2']= math.sqrt((targets['sub1_location_x'][i]-0.9583333)**2+ (targets['sub1_location_y'][i]-0.9583333)**2+(targets['sub1_location_z'][i]-0.355174)**2)

    targets.loc[i,'sub2_distance2']= math.sqrt((targets['sub2_location_x'][i]-0.9583333)**2+(targets['sub2_location_y'][i]-0.9583333)**2+(targets['sub2_location_z'][i]-0.355174)**2)

    targets.loc[i,'sub3_distance2']= math.sqrt((targets['sub3_location_x'][i]-0.9583333)**2+(targets['sub3_location_y'][i]-0.9583333)**2+(targets['sub3_location_z'][i]-0.355174)**2)

    targets.loc[i,'vac1_distance2']= math.sqrt((targets['vacancy1_location_x'][i]-0.9583333)**2+(targets['vacancy1_location_y'][i]-0.9583333)**2+(targets['vacancy1_location_z'][i]-0.355174)**2)

    targets.loc[i,'vac2_distance2']= math.sqrt((targets['vacancy2_location_x'][i]-0.9583333)**2+(targets['vacancy2_location_y'][i]-0.9583333)**2+(targets['vacancy2_location_z'][i]-0.355174)**2)

    targets.loc[i,'vac3_distance2']= math.sqrt((targets['vacancy3_location_x'][i]-0.9583333)**2+(targets['vacancy3_location_y'][i]-0.9583333)**2+(targets['vacancy3_location_z'][i]-0.355174)**2)


##new coordinate system
targets['sub1_distance3']=0
targets['sub2_distance3']=0

targets['vac1_distance3']=0
targets['vac2_distance3']=0
targets['vac3_distance3']=0

for i in range(0,len(targets['sub1_distance3'])):
    targets.loc[i,'sub1_distance3']= math.sqrt((targets['sub1_location_x'][i]-0.04166667)**2+ (targets['sub1_location_y'][i]-0.9583333)**2+(targets['sub1_location_z'][i]-0.355174)**2)

    targets.loc[i,'sub2_distance3']= math.sqrt((targets['sub2_location_x'][i]-0.04166667)**2+(targets['sub2_location_y'][i]-0.9583333)**2+(targets['sub2_location_z'][i]-0.355174)**2)

    targets.loc[i,'sub3_distance3']= math.sqrt((targets['sub3_location_x'][i]-0.04166667)**2+(targets['sub3_location_y'][i]-0.9583333)**2+(targets['sub3_location_z'][i]-0.355174)**2)

    targets.loc[i,'vac1_distance3']= math.sqrt((targets['vacancy1_location_x'][i]-0.04166667)**2+(targets['vacancy1_location_y'][i]-0.9583333)**2+(targets['vacancy1_location_z'][i]-0.355174)**2)

    targets.loc[i,'vac2_distance3']= math.sqrt((targets['vacancy2_location_x'][i]-0.04166667)**2+(targets['vacancy2_location_y'][i]-0.9583333)**2+(targets['vacancy2_location_z'][i]-0.355174)**2)

    targets.loc[i,'vac3_distance3']= math.sqrt((targets['vacancy3_location_x'][i]-0.04166667)**2+(targets['vacancy3_location_y'][i]-0.9583333)**2+(targets['vacancy3_location_z'][i]-0.355174)**2)


targets['sub1_element'][targets['sub1_element']==0]= "None"
targets['sub2_element'][targets['sub2_element']==0]= "None"
targets['sub3_element'][targets['sub3_element']==0]= "None"

# Read in Private Data

import pandas as pd 
import numpy as np

#read in data
results=pd.read_csv('C:/Users/aerob/Downloads/testingdata.csv')


results.id=(results.id).astype(str)

results = results.rename(columns={'id': 'X_id'})

import os
os.chdir("C:/Users/aerob/Downloads/IDAO 2022 (1)/IDAO 2022/dichalcogenides_private/structures/")

import json
 
# Opening JSON file
f = open(results.X_id[0]+".json")
 
# returns JSON object as
# a dictionary
data = json.load(f)


results['a']=0
results['b']=0
results['c']=0
results['alpha']=0
results['beta']=0
results['gamma']=0
results['volume']=0


for id in range(0,len(results['a'])):
    # Opening JSON file
    f = open(results.X_id[id]+".json")

    # returns JSON object as
    # a dictionary
    samplist = json.load(f)
    results['a'][id]=samplist['lattice']['a']
    results['b'][id]=samplist['lattice']['b']
    results['c'][id]=samplist['lattice']['c']
    results['alpha'][id]=samplist['lattice']['alpha']
    results['beta'][id]=samplist['lattice']['beta']
    results['gamma'][id]=samplist['lattice']['gamma']
    results['volume'][id]=samplist['lattice']['volume']


def listoflabels(id):
    listinit=[]
    f = open(results.X_id[id]+".json")

    # returns JSON object as
    # a dictionary
    samplist = json.load(f)
    for i in range(0, len(samplist['sites'])):
        listinit.append(samplist['sites'][i]['label'])
    
    return pd.Series(listinit)


#substitutional impurities
results['num_substitution']=0
for id in range(0,len(results['num_substitution'])):
    f = open(results.X_id[id]+".json")

    # returns JSON object as
    # a dictionary
    samplist = json.load(f)
    y=(listoflabels(id)).value_counts()
    if(len(y.index) > 2):
        results['num_substitution'][id]=sum(y[2:])
    



#type of substitution
results['sub1_location_x']=0
results['sub1_location_y']=0
results['sub1_location_z']=0
results['sub1_element']=0

results['sub2_location_x']=0
results['sub2_location_y']=0
results['sub2_location_z']=0
results['sub2_element']= 0

results['sub3_location_x']=0
results['sub3_location_y']=0
results['sub3_location_z']=0
results['sub3_element']= 0
for id in range(0, len(results['sub1_element'])) :
    f = open(results.X_id[id]+".json")

    # returns JSON object as
    # a dictionary
    samplist = json.load(f)  
    x=listoflabels(id).value_counts()


    if results['num_substitution'][id]==1 :
        results.loc[id,'sub1_location_x']=samplist['sites'][np.where(listoflabels(id)==listoflabels(id).value_counts().index[2])[0][0]]['abc'][0]
        results.loc[id,'sub1_location_y']=samplist['sites'][np.where(listoflabels(id)==listoflabels(id).value_counts().index[2])[0][0]]['abc'][1]
        results.loc[id,'sub1_location_z']=samplist['sites'][np.where(listoflabels(id)==listoflabels(id).value_counts().index[2])[0][0]]['abc'][2]
        results.loc[id, 'sub1_element']=listoflabels(id).value_counts().index[2]

    if results['num_substitution'][id]==2:
        results.loc[id,'sub1_location_x']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[0]]['abc'][0]
        results.loc[id,'sub1_location_y']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[0]]['abc'][1]
        results.loc[id,'sub1_location_z']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[0]]['abc'][2]
     
        results.loc[id, 'sub1_element']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[0]]['label']

        results.loc[id,'sub2_location_x']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[1]]['abc'][0]
        results.loc[id,'sub2_location_y']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[1]]['abc'][1]
        results.loc[id,'sub2_location_z']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[1]]['abc'][2]
        results.loc[id, 'sub2_element']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[1]]['label']

  
    if results['num_substitution'][id]==3:
        results.loc[id,'sub1_location_x']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[0]]['abc'][0]
        results.loc[id,'sub1_location_y']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[0]]['abc'][1]
        results.loc[id,'sub1_location_z']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[0]]['abc'][2]
     
        results.loc[id, 'sub1_element']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[0]]['label']

        results.loc[id,'sub2_location_x']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[1]]['abc'][0]
        results.loc[id,'sub2_location_y']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[1]]['abc'][1]
        results.loc[id,'sub2_location_z']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[1]]['abc'][2]

        results.loc[id, 'sub2_element']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[1]]['label']


        results.loc[id,'sub3_location_x']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[2]]['abc'][0]
        results.loc[id,'sub3_location_y']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[2]]['abc'][1]
        results.loc[id,'sub3_location_z']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[2]]['abc'][2]

        results.loc[id, 'sub3_element']=samplist['sites'][np.intersect1d(np.where(listoflabels(id) != "Mo")[0], (np.where(listoflabels(id) != "S"))[0])[2]]['label']







results['num_vacancies']=0
total_points=192
for id in range(0, len(results['num_vacancies'])):
    f = open(results.X_id[id]+".json")

    # returns JSON object as
    # a dictionary
    samplist = json.load(f)    
    results['num_vacancies'][id]=total_points-len(samplist['sites'])


def listofcoords(id):
    listinit=[]
    f = open(results.X_id[id]+".json")

    # returns JSON object as
    # a dictionary
    samplist = json.load(f)
    for i in range(0, len(samplist['sites'])):
        listinit.append(samplist['sites'][i]['abc'])
    
    return pd.Series(listinit)

len(listofcoords(2))

results['mo_vacancies']=0
results['s_vacancies']=0
results['s_bottom_vacancies']=0
results['s_top_vacancies']=0

results['vac1_element']=0
results['vacancy1_location_x']=0
results['vacancy1_location_y']=0
results['vacancy1_location_z']=0

results['vac2_element']=0
results['vacancy2_location_x']=0
results['vacancy2_location_y']=0
results['vacancy2_location_z']=0

results['vac3_element']=0
results['vacancy3_location_x']=0
results['vacancy3_location_y']=0
results['vacancy3_location_z']=0

good_lattice=listofcoords(2)
good_df=pd.DataFrame(list(good_lattice), columns = ['x', 'y','z'])


for id in list(np.where(results['num_vacancies'] > 0)[0]):
    f = open(results.X_id[id]+".json")

    # returns JSON object as
    # a dictionary
    samplist = json.load(f)
    bad_lattice=listofcoords(id)
    bad_df=pd.DataFrame(list(bad_lattice), columns = ['x', 'y','z'])
    difference_df=pd.concat([good_df,bad_df]).drop_duplicates(keep=False)
    difference_df=pd.DataFrame(difference_df).reset_index()
    
#GET SET DIFFERENCES CORRECT
    for i in range(0,results['num_vacancies'][id]):
        if difference_df['z'][i] < 0.2:
            results['s_vacancies'][id]=results['s_vacancies'][id]+1
            results['s_bottom_vacancies'][id]=results['s_bottom_vacancies'][id]+1
        
        if difference_df['z'][i] > 0.3:
            results['s_vacancies'][id]=results['s_vacancies'][id]+1
            results['s_top_vacancies'][id]=results['s_bottom_vacancies'][id]+1
        
        if difference_df['z'][i] > 0.2 and difference_df['z'][i] < 0.3:
            results['mo_vacancies'][id]=results['mo_vacancies'][id]+1
        
        
  
  
    if results['num_vacancies'][id]==1:
        results.loc[id,'vacancy1_location_x']=difference_df['x'][0]
        results.loc[id,'vacancy1_location_y']=difference_df['y'][0]
        results.loc[id,'vacancy1_location_z']=difference_df['z'][0]

      

    if results['num_vacancies'][id]==2:
        results.loc[id,'vacancy1_location_x']=difference_df['x'][0]
        results.loc[id,'vacancy1_location_y']=difference_df['y'][0]
        results.loc[id,'vacancy1_location_z']=difference_df['z'][0]
        results.loc[id,'vacancy2_location_x']=difference_df['x'][1]
        results.loc[id,'vacancy2_location_y']=difference_df['y'][1]
        results.loc[id,'vacancy2_location_z']=difference_df['z'][1]
      
  
    if results['num_vacancies'][id]==3:
        results.loc[id,'vacancy1_location_x']=difference_df['x'][0]
        results.loc[id,'vacancy1_location_y']=difference_df['y'][0]
        results.loc[id,'vacancy1_location_z']=difference_df['z'][0]
        results.loc[id,'vacancy2_location_x']=difference_df['x'][1]
        results.loc[id,'vacancy2_location_y']=difference_df['y'][1]
        results.loc[id,'vacancy2_location_z']=difference_df['z'][1]
        results.loc[id,'vacancy3_location_x']=difference_df['x'][2]
        results.loc[id,'vacancy3_location_y']=difference_df['y'][2]
        results.loc[id,'vacancy3_location_z']=difference_df['z'][2]
        


difference_df

for id in range(0,len(results['vacancy1_location_z'])):
    if results['vacancy1_location_z'][id] > 0.1 and results['vacancy1_location_z'][id] < 0.2:
        results.loc[id, 'vac1_element'] = "S"
    elif results['vacancy1_location_z'][id] > 0.2 and results['vacancy1_location_z'][id] < 0.3:
        results.loc[id, 'vac1_element'] = "Mo"
    elif results['vacancy1_location_z'][id] > 0.3 and results['vacancy1_location_z'][id] < 0.4:
        results.loc[id, 'vac1_element'] = "S"
    else:
        results.loc[id, 'vac1_element'] = "None"

    if results['vacancy2_location_z'][id] > 0.1 and results['vacancy2_location_z'][id] < 0.2:
        results.loc[id, 'vac2_element'] = "S"
    elif results['vacancy2_location_z'][id] > 0.2 and results['vacancy2_location_z'][id] < 0.3:
        results.loc[id, 'vac2_element'] ="Mo"
    elif results['vacancy2_location_z'][id] > 0.3 and results['vacancy2_location_z'][id] < 0.4:
        results.loc[id, 'vac2_element'] = "S"
    else:
        results.loc[id, 'vac2_element'] = "None"
  
    if results['vacancy3_location_z'][id] > 0.1 and results['vacancy3_location_z'][id] < 0.2:
        results.loc[id, 'vac3_element'] = "S"
    elif results['vacancy3_location_z'][id] > 0.2 and results['vacancy3_location_z'][id] < 0.3:
        results.loc[id, 'vac3_element'] = "Mo"
    elif results['vacancy3_location_z'][id] > 0.3 and results['vacancy3_location_z'][id] < 0.4:
        results.loc[id, 'vac3_element'] = "S"
    else:
        results.loc[id, 'vac3_element'] = "None"

  
  



##new coordinate system
import math
results['sub1_distance']=0
results['sub2_distance']=0

results['vac1_distance']=0
results['vac2_distance']=0
results['vac3_distance']=0

for i in range(0,len(results['sub1_distance'])):
    results.loc[i,'sub1_distance']= math.sqrt(results['sub1_location_x'][i]**2+ results['sub1_location_y'][i]**2+results['sub1_location_z'][i]**2)

    results.loc[i,'sub2_distance']= math.sqrt(results['sub2_location_x'][i]**2+results['sub2_location_y'][i]**2+results['sub2_location_z'][i]**2)

    results.loc[i,'sub3_distance']= math.sqrt(results['sub3_location_x'][i]**2+results['sub3_location_y'][i]**2+results['sub3_location_z'][i]**2)

    results.loc[i,'vac1_distance']= math.sqrt(results['vacancy1_location_x'][i]**2+results['vacancy1_location_y'][i]**2+results['vacancy1_location_z'][i]**2)

    results.loc[i,'vac2_distance']= math.sqrt(results['vacancy2_location_x'][i]**2+results['vacancy2_location_y'][i]**2+results['vacancy2_location_z'][i]**2)
    results.loc[i,'vac3_distance']= math.sqrt(results['vacancy3_location_x'][i]**2+results['vacancy3_location_y'][i]**2+results['vacancy3_location_z'][i]**2)



##new coordinate system
import math
results['sub1_distance2']=0
results['sub2_distance2']=0

results['vac1_distance2']=0
results['vac2_distance2']=0
results['vac3_distance2']=0

for i in range(0,len(results['sub1_distance2'])):
    results.loc[i,'sub1_distance2']= math.sqrt((results['sub1_location_x'][i]-0.9583333)**2+ (results['sub1_location_y'][i]-0.9583333)**2+(results['sub1_location_z'][i]-0.355174)**2)

    results.loc[i,'sub2_distance2']= math.sqrt((results['sub2_location_x'][i]-0.9583333)**2+(results['sub2_location_y'][i]-0.9583333)**2+(results['sub2_location_z'][i]-0.355174)**2)

    results.loc[i,'sub3_distance2']= math.sqrt((results['sub3_location_x'][i]-0.9583333)**2+(results['sub3_location_y'][i]-0.9583333)**2+(results['sub3_location_z'][i]-0.355174)**2)

    results.loc[i,'vac1_distance2']= math.sqrt((results['vacancy1_location_x'][i]-0.9583333)**2+(results['vacancy1_location_y'][i]-0.9583333)**2+(results['vacancy1_location_z'][i]-0.355174)**2)

    results.loc[i,'vac2_distance2']= math.sqrt((results['vacancy2_location_x'][i]-0.9583333)**2+(results['vacancy2_location_y'][i]-0.9583333)**2+(results['vacancy2_location_z'][i]-0.355174)**2)

    results.loc[i,'vac3_distance2']= math.sqrt((results['vacancy3_location_x'][i]-0.9583333)**2+(results['vacancy3_location_y'][i]-0.9583333)**2+(results['vacancy3_location_z'][i]-0.355174)**2)


##new coordinate system
results['sub1_distance3']=0
results['sub2_distance3']=0

results['vac1_distance3']=0
results['vac2_distance3']=0
results['vac3_distance3']=0

for i in range(0,len(results['sub1_distance3'])):
    results.loc[i,'sub1_distance3']= math.sqrt((results['sub1_location_x'][i]-0.04166667)**2+ (results['sub1_location_y'][i]-0.9583333)**2+(results['sub1_location_z'][i]-0.355174)**2)

    results.loc[i,'sub2_distance3']= math.sqrt((results['sub2_location_x'][i]-0.04166667)**2+(results['sub2_location_y'][i]-0.9583333)**2+(results['sub2_location_z'][i]-0.355174)**2)

    results.loc[i,'sub3_distance3']= math.sqrt((results['sub3_location_x'][i]-0.04166667)**2+(results['sub3_location_y'][i]-0.9583333)**2+(results['sub3_location_z'][i]-0.355174)**2)

    results.loc[i,'vac1_distance3']= math.sqrt((results['vacancy1_location_x'][i]-0.04166667)**2+(results['vacancy1_location_y'][i]-0.9583333)**2+(results['vacancy1_location_z'][i]-0.355174)**2)

    results.loc[i,'vac2_distance3']= math.sqrt((results['vacancy2_location_x'][i]-0.04166667)**2+(results['vacancy2_location_y'][i]-0.9583333)**2+(results['vacancy2_location_z'][i]-0.355174)**2)

    results.loc[i,'vac3_distance3']= math.sqrt((results['vacancy3_location_x'][i]-0.04166667)**2+(results['vacancy3_location_y'][i]-0.9583333)**2+(results['vacancy3_location_z'][i]-0.355174)**2)


results['sub1_element'][results['sub1_element']==0]= "None"
results['sub2_element'][results['sub2_element']==0]= "None"
results['sub3_element'][results['sub3_element']==0]= "None"

# One Vacancy Subset

One_vac=targets[targets['num_vacancies']==1]
One_vac=One_vac.reset_index()
One_vac.to_csv('One_vac_224.csv')

x=One_vac.drop(One_vac.columns[0:2],axis=1)
data=x.drop(['a', 'b','c','alpha','beta','gamma','volume'], axis=1)

One_vac_old=pd.read_csv('C:\\Users\\aerob\\Desktop\\One_vac_old.csv')

data_try=data[list(One_vac_old.columns)]

One_vac_old=One_vac_old.round(6)
data_try=data_try.round(6)

ne_stacked = (One_vac_old != data_try).stack()
changed = ne_stacked[ne_stacked]
changed.index.names = ['id', 'col']
changed

difference_locations = np.where(One_vac_old != data_try)
changed_from = One_vac_old.values[difference_locations]
changed_to = data_try.values[difference_locations]
y=pd.DataFrame({'from': changed_from, 'to': changed_to}, index=changed.index)


y

# Two Vacancy Subset

Two_vac=targets[targets['num_vacancies']==2]
Two_vac=Two_vac.reset_index()


Two_vac.columns



Two_vac['vac_distance_between']=0

for i in range(0, len(Two_vac['vac_distance_between'])):
    abs_x=(Two_vac['vacancy1_location_x'][i]-Two_vac['vacancy2_location_x'][i])**2
    abs_y=(Two_vac['vacancy1_location_y'][i]-Two_vac['vacancy2_location_y'][i])**2
    abs_z=(Two_vac['vacancy1_location_z'][i]-Two_vac['vacancy2_location_z'][i])**2
  
    Two_vac.loc[i,'vac_distance_between']=math.sqrt(abs_x+abs_y+abs_z)


Two_vac['vac_distance_between_xy']=0

for i in range(0, len(Two_vac['vac_distance_between'])):
    abs_x=(Two_vac['vacancy1_location_x'][i]-Two_vac['vacancy2_location_x'][i])**2
    abs_y=(Two_vac['vacancy1_location_y'][i]-Two_vac['vacancy2_location_y'][i])**2
  
    Two_vac.loc[i,'vac_distance_between_xy']=math.sqrt(abs_x+abs_y)


Two_vac['vac_distance_between_xz']=0

for i in range(0, len(Two_vac['vac_distance_between'])):
    abs_x=(Two_vac['vacancy1_location_x'][i]-Two_vac['vacancy2_location_x'][i])**2
    abs_y=(Two_vac['vacancy1_location_z'][i]-Two_vac['vacancy2_location_z'][i])**2
  
    Two_vac.loc[i,'vac_distance_between_xz']=math.sqrt(abs_x+abs_y)

Two_vac['vac_distance_between_yz']=0

for i in range(0, len(Two_vac['vac_distance_between'])):
    abs_x=(Two_vac['vacancy1_location_y'][i]-Two_vac['vacancy2_location_y'][i])**2
    abs_y=(Two_vac['vacancy1_location_z'][i]-Two_vac['vacancy2_location_z'][i])**2
  
    Two_vac.loc[i,'vac_distance_between_yz']=math.sqrt(abs_x+abs_y)

Two_vac["sub1_element"].replace({2.58: "S", 2.36: "W", 2.55:"Se", 2.16:"Mo", 0:"None"}, inplace=True)
Two_vac["sub2_element"].replace({2.58: "S", 2.36: "W", 2.55:"Se", 2.16:"Mo", 0:"None"}, inplace=True)
Two_vac["sub3_element"].replace({2.58: "S", 2.36: "W", 2.55:"Se", 2.16:"Mo", 0:"None"}, inplace=True)
Two_vac["vac1_element"].replace({2.58: "S", 2.36: "W", 2.55:"Se", 2.16:"Mo", 0:"None"}, inplace=True)
Two_vac["vac2_element"].replace({2.58: "S", 2.36: "W", 2.55:"Se", 2.16:"Mo", 0:"None"}, inplace=True)
Two_vac["vac3_element"].replace({2.58: "S", 2.36: "W", 2.55:"Se", 2.16:"Mo", 0:"None"}, inplace=True)

Two_vac['origin_xy_angle']=0
for i in range(0, len(Two_vac['vacancy2_location_x'])):
    a = np.array([Two_vac['vacancy2_location_x'][i],Two_vac['vacancy2_location_y'][i]])
    b = np.array([0,0])
    c = np.array([Two_vac['vacancy1_location_x'][i],Two_vac['vacancy1_location_y'][i]])

    ba = a - b
    bc = c - b

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)

    Two_vac.loc[i,'origin_xy_angle']=(np.degrees(angle))

for i in np.where(Two_vac['origin_xy_angle'].isnull()):
    Two_vac.loc[i, 'origin_xy_angle']=0

Two_vac['sub_vacs_xy_angle']=0
for i in range(0, len(Two_vac['vacancy2_location_x'])):
    a = np.array([Two_vac['vacancy2_location_x'][i],Two_vac['vacancy2_location_y'][i]])
    b = np.array([Two_vac['sub1_location_x'][i],Two_vac['sub1_location_y'][i]])
    c = np.array([Two_vac['vacancy3_location_x'][i],Two_vac['vacancy3_location_y'][i]])

    ba = a - b
    bc = c - b

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)

    Two_vac.loc[i,'sub_vacs_xy_angle']=(np.degrees(angle))

for i in np.where(Two_vac['sub_vacs_xy_angle'].isnull()):
    Two_vac.loc[i, 'sub_vacs_xy_angle']=0



Two_vac['xz_slope']=0
for i in range(0, len(Two_vac['vacancy2_location_x'])):
    a = abs(Two_vac['vacancy1_location_z'][i]-Two_vac['vacancy2_location_z'][i])
    b = abs(Two_vac['vacancy1_location_x'][i]-Two_vac['vacancy2_location_x'][i])
    if b==0:
        Two_vac.loc[i,'xz_slope']=1
    elif a==0:
        Two_vac.loc[i,'xz_slope']=0
    else:
        Two_vac.loc[i,'xz_slope']=a/b

Two_vac.to_csv('Two_vac_224.csv')

Two_vac['xz_angle']=0
for i in range(0, len(Two_vac['vacancy2_location_x'])):
    a = np.array([Two_vac['vacancy1_location_x'][i],Two_vac['vacancy2_location_z'][i]])
    b = np.array([Two_vac['vacancy1_location_x'][i],Two_vac['vacancy1_location_z'][i]])
    c = np.array([Two_vac['vacancy2_location_x'][i],Two_vac['vacancy2_location_z'][i]])

    ba = a - b
    bc = c - b

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)
    if Two_vac['vacancy2_location_z'][i]==Two_vac['vacancy1_location_z'][i]:
        Two_vac.loc[i,'xz_angle']=90.0
    else:
        Two_vac.loc[i,'xz_angle']=np.degrees(angle)
    
    


Two_vac['xy_angle']=0
for i in range(0, len(Two_vac['vacancy2_location_x'])):
    a = np.array([Two_vac['vacancy1_location_x'][i],Two_vac['vacancy2_location_y'][i]])
    b = np.array([Two_vac['vacancy1_location_x'][i],Two_vac['vacancy1_location_y'][i]])
    c = np.array([Two_vac['vacancy2_location_x'][i],Two_vac['vacancy2_location_y'][i]])

    ba = a - b
    bc = c - b

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)
    if Two_vac['vacancy2_location_y'][i]==Two_vac['vacancy1_location_y'][i]:
        Two_vac.loc[i,'xz_angle']=90.0
    else:
        Two_vac.loc[i,'xz_angle']=np.degrees(angle)
    
    

Two_vac['yz_angle']=0
for i in range(0, len(Two_vac['vacancy2_location_x'])):
    a = np.array([Two_vac['vacancy1_location_y'][i],Two_vac['vacancy2_location_z'][i]])
    b = np.array([Two_vac['vacancy1_location_y'][i],Two_vac['vacancy1_location_z'][i]])
    c = np.array([Two_vac['vacancy2_location_y'][i],Two_vac['vacancy2_location_z'][i]])

    ba = a - b
    bc = c - b

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)
    if Two_vac['vacancy2_location_z'][i]==Two_vac['vacancy1_location_z'][i]:
        Two_vac.loc[i,'xz_angle']=90.0
    else:
        Two_vac.loc[i,'xz_angle']=np.degrees(angle)
    
    

Two_vac['vac1_sub_distance_xy']=0
for i in range(0, len(Two_vac['vacancy2_location_x'])):
    Two_vac.loc[i, 'vac1_sub_distance_xy']= (Two_vac['vacancy1_location_x'][i]-Two_vac['sub1_location_x'][i])**2

Two_vac['vac2_sub_distance_xy']=0
for i in range(0, len(Two_vac['vacancy2_location_x'])):
    Two_vac.loc[i, 'vac2_sub_distance_xy']= (Two_vac['vacancy2_location_x'][i]-Two_vac['sub1_location_x'][i])**2+(Two_vac['vacancy2_location_y'][i]-Two_vac['sub1_location_y'][i])**2

#use this for two vacancies
x=Two_vac.drop(One_vac.columns[0:2],axis=1)
data=x.drop(['a', 'b','c','alpha','beta','gamma','volume'], axis=1)
data=data.round(decimals=6)
one_hot_encoded_data = pd.get_dummies(data, columns = ['sub1_element', 'vac1_element', 'vac2_element','vac3_element', 'sub2_element', 'sub3_element']) #use this setting for twoVacs
features = one_hot_encoded_data[['vac_distance_between', 'vac_distance_between_xy', 'vac1_element_Mo']]

Two_vac.to_csv('C:/Users/aerob/Downloads/Two_vac_227.csv')

# Three Vacancy Subset

Three_vac=targets[targets['num_vacancies']==3]
Three_vac=Three_vac.reset_index()

Three_vac['mo_s1_distance']=0
Three_vac['mo_s2_distance']=0
Three_vac['s1_s2_distance']=0

for i in range(0,len(Three_vac['mo_s1_distance'])):
    if Three_vac['vac1_element'][i]== "Mo":
        abs_mo_s1_x=(Three_vac['vacancy1_location_x'][i]-Three_vac['vacancy2_location_x'][i])**2
        abs_mo_s1_y=(Three_vac['vacancy1_location_y'][i]-Three_vac['vacancy2_location_y'][i])**2
        abs_mo_s1_z=(Three_vac['vacancy1_location_z'][i]-Three_vac['vacancy2_location_z'][i])**2

        Three_vac.loc[i,'mo_s1_distance']= (abs_mo_s1_x+abs_mo_s1_y+abs_mo_s1_z)**0.5

        abs_mo_s2_x=(Three_vac['vacancy1_location_x'][i]-Three_vac['vacancy3_location_x'][i])**2
        abs_mo_s2_y=(Three_vac['vacancy1_location_y'][i]-Three_vac['vacancy3_location_y'][i])**2
        abs_mo_s2_z=(Three_vac['vacancy1_location_z'][i]-Three_vac['vacancy3_location_z'][i])**2

        Three_vac.loc[i,'mo_s2_distance']= (abs_mo_s2_x+abs_mo_s2_y+abs_mo_s2_z)**0.5


        abs_s1_s2_x=abs(Three_vac['vacancy2_location_x'][i]-Three_vac['vacancy3_location_x'][i])**2
        abs_s1_s2_y=abs(Three_vac['vacancy2_location_y'][i]-Three_vac['vacancy3_location_y'][i])**2
        abs_s1_s2_z=abs(Three_vac['vacancy2_location_z'][i]-Three_vac['vacancy3_location_z'][i])**2

        Three_vac.loc[i,'s1_s2_distance']= (abs_s1_s2_x+abs_s1_s2_y+abs_s1_s2_z)**0.5
  
  
    if Three_vac['vac2_element'][i]== "Mo":
        abs_mo_s1_x=abs(Three_vac['vacancy1_location_x'][i]-Three_vac['vacancy2_location_x'][i])**2
        abs_mo_s1_y=abs(Three_vac['vacancy1_location_y'][i]-Three_vac['vacancy2_location_y'][i])**2
        abs_mo_s1_z=abs(Three_vac['vacancy1_location_z'][i]-Three_vac['vacancy2_location_z'][i])**2

        Three_vac.loc[i,'mo_s1_distance']= (abs_mo_s1_x+abs_mo_s1_y+abs_mo_s1_z)**0.5

        abs_mo_s2_x=abs(Three_vac['vacancy2_location_x'][i]-Three_vac['vacancy3_location_x'][i])**2
        abs_mo_s2_y=abs(Three_vac['vacancy2_location_y'][i]-Three_vac['vacancy3_location_y'][i])**2
        abs_mo_s2_z=abs(Three_vac['vacancy2_location_z'][i]-Three_vac['vacancy3_location_z'][i])**2

        Three_vac.loc[i,'mo_s2_distance']= (abs_mo_s2_x+abs_mo_s2_y+abs_mo_s2_z)**0.5

        abs_s1_s2_x=abs(Three_vac['vacancy1_location_x'][i]-Three_vac['vacancy3_location_x'][i])**2
        abs_s1_s2_y=abs(Three_vac['vacancy1_location_y'][i]-Three_vac['vacancy3_location_y'][i])**2
        abs_s1_s2_z=abs(Three_vac['vacancy1_location_z'][i]-Three_vac['vacancy3_location_z'][i])**2

        Three_vac.loc[i,'s1_s2_distance']= (abs_s1_s2_x+abs_s1_s2_y+abs_s1_s2_z)**0.5




    if Three_vac['vac3_element'][i]== "Mo":
        abs_mo_s1_x=abs(Three_vac['vacancy1_location_x'][i]-Three_vac['vacancy3_location_x'][i])**2
        abs_mo_s1_y=abs(Three_vac['vacancy1_location_y'][i]-Three_vac['vacancy3_location_y'][i])**2
        abs_mo_s1_z=abs(Three_vac['vacancy1_location_z'][i]-Three_vac['vacancy3_location_z'][i])**2

        Three_vac.loc[i,'mo_s1_distance']= (abs_mo_s1_x+abs_mo_s1_y+abs_mo_s1_z)**0.5

        abs_mo_s2_x=abs(Three_vac['vacancy2_location_x'][i]-Three_vac['vacancy3_location_x'][i])**2
        abs_mo_s2_y=abs(Three_vac['vacancy2_location_y'][i]-Three_vac['vacancy3_location_y'][i])**2
        abs_mo_s2_z=abs(Three_vac['vacancy2_location_z'][i]-Three_vac['vacancy3_location_z'][i])**2

        Three_vac.loc[i,'mo_s2_distance']= (abs_mo_s2_x+abs_mo_s2_y+abs_mo_s2_z)**0.5

        abs_s1_s2_x=abs(Three_vac['vacancy1_location_x'][i]-Three_vac['vacancy2_location_x'][i])**2
        abs_s1_s2_y=abs(Three_vac['vacancy1_location_y'][i]-Three_vac['vacancy2_location_y'][i])**2
        abs_s1_s2_z=abs(Three_vac['vacancy1_location_z'][i]-Three_vac['vacancy2_location_z'][i])**2

        Three_vac.loc[i,'s1_s2_distance']= (abs_s1_s2_x+abs_s1_s2_y+abs_s1_s2_z)**0.5

  
  
  



Three_vac['mo_s1_xydistance']=0
Three_vac['mo_s2_xydistance']=0
Three_vac['s1_s2_xydistance']=0

for i in range(0,len(Three_vac['mo_s1_xydistance'])):
    if Three_vac['vac1_element'][i]== "Mo":
        abs_mo_s1_x=(Three_vac['vacancy1_location_x'][i]-Three_vac['vacancy2_location_x'][i])**2
        abs_mo_s1_y=(Three_vac['vacancy1_location_y'][i]-Three_vac['vacancy2_location_y'][i])**2

        Three_vac.loc[i,'mo_s1_xydistance']= (abs_mo_s1_x+abs_mo_s1_y)**0.5

        abs_mo_s2_x=(Three_vac['vacancy1_location_x'][i]-Three_vac['vacancy3_location_x'][i])**2
        abs_mo_s2_y=(Three_vac['vacancy1_location_y'][i]-Three_vac['vacancy3_location_y'][i])**2

        Three_vac.loc[i,'mo_s2_xydistance']= (abs_mo_s2_x+abs_mo_s2_y)**0.5


        abs_s1_s2_x=abs(Three_vac['vacancy2_location_x'][i]-Three_vac['vacancy3_location_x'][i])**2
        abs_s1_s2_y=abs(Three_vac['vacancy2_location_y'][i]-Three_vac['vacancy3_location_y'][i])**2

        Three_vac.loc[i,'s1_s2_xydistance']= (abs_s1_s2_x+abs_s1_s2_y)**0.5
  
  
    if Three_vac['vac2_element'][i]== "Mo":
        abs_mo_s1_x=abs(Three_vac['vacancy1_location_x'][i]-Three_vac['vacancy2_location_x'][i])**2
        abs_mo_s1_y=abs(Three_vac['vacancy1_location_y'][i]-Three_vac['vacancy2_location_y'][i])**2

        Three_vac.loc[i,'mo_s1_xydistance']= (abs_mo_s1_x+abs_mo_s1_y)**0.5

        abs_mo_s2_x=abs(Three_vac['vacancy2_location_x'][i]-Three_vac['vacancy3_location_x'][i])**2
        abs_mo_s2_y=abs(Three_vac['vacancy2_location_y'][i]-Three_vac['vacancy3_location_y'][i])**2

        Three_vac.loc[i,'mo_s2_xydistance']= (abs_mo_s2_x+abs_mo_s2_y)**0.5

        abs_s1_s2_x=abs(Three_vac['vacancy1_location_x'][i]-Three_vac['vacancy3_location_x'][i])**2
        abs_s1_s2_y=abs(Three_vac['vacancy1_location_y'][i]-Three_vac['vacancy3_location_y'][i])**2

        Three_vac.loc[i,'s1_s2_xydistance']= (abs_s1_s2_x+abs_s1_s2_y)**0.5




    if Three_vac['vac3_element'][i]== "Mo":
        abs_mo_s1_x=abs(Three_vac['vacancy1_location_x'][i]-Three_vac['vacancy3_location_x'][i])**2
        abs_mo_s1_y=abs(Three_vac['vacancy1_location_y'][i]-Three_vac['vacancy3_location_y'][i])**2

        Three_vac.loc[i,'mo_s1_xydistance']= (abs_mo_s1_x+abs_mo_s1_y)**0.5

        abs_mo_s2_x=abs(Three_vac['vacancy2_location_x'][i]-Three_vac['vacancy3_location_x'][i])**2
        abs_mo_s2_y=abs(Three_vac['vacancy2_location_y'][i]-Three_vac['vacancy3_location_y'][i])**2

        Three_vac.loc[i,'mo_s2_xydistance']= (abs_mo_s2_x+abs_mo_s2_y)**0.5

        abs_s1_s2_x=abs(Three_vac['vacancy1_location_x'][i]-Three_vac['vacancy2_location_x'][i])**2
        abs_s1_s2_y=abs(Three_vac['vacancy1_location_y'][i]-Three_vac['vacancy2_location_y'][i])**2

        Three_vac.loc[i,'s1_s2_xydistance']= (abs_s1_s2_x+abs_s1_s2_y)**0.5

  
  
  

Three_vac['mo_s1_xzdistance']=0
Three_vac['mo_s2_xzdistance']=0
Three_vac['s1_s2_xzdistance']=0

for i in range(0,len(Three_vac['mo_s1_xzdistance'])):
    if Three_vac['vac1_element'][i]== "Mo":
        abs_mo_s1_x=(Three_vac['vacancy1_location_x'][i]-Three_vac['vacancy2_location_x'][i])**2
        abs_mo_s1_z=(Three_vac['vacancy1_location_z'][i]-Three_vac['vacancy2_location_z'][i])**2

        Three_vac.loc[i,'mo_s1_xzdistance']= (abs_mo_s1_x+abs_mo_s1_z)**0.5

        abs_mo_s2_x=(Three_vac['vacancy1_location_x'][i]-Three_vac['vacancy3_location_x'][i])**2
        abs_mo_s2_z=(Three_vac['vacancy1_location_z'][i]-Three_vac['vacancy3_location_z'][i])**2

        Three_vac.loc[i,'mo_s2_xzdistance']= (abs_mo_s2_x+abs_mo_s2_z)**0.5


        abs_s1_s2_x=abs(Three_vac['vacancy2_location_x'][i]-Three_vac['vacancy3_location_x'][i])**2
        abs_s1_s2_z=abs(Three_vac['vacancy2_location_z'][i]-Three_vac['vacancy3_location_z'][i])**2

        Three_vac.loc[i,'s1_s2_xzdistance']= (abs_s1_s2_x+abs_s1_s2_z)**0.5
  
  
    if Three_vac['vac2_element'][i]== "Mo":
        abs_mo_s1_x=abs(Three_vac['vacancy1_location_x'][i]-Three_vac['vacancy2_location_x'][i])**2
        abs_mo_s1_z=abs(Three_vac['vacancy1_location_z'][i]-Three_vac['vacancy2_location_z'][i])**2

        Three_vac.loc[i,'mo_s1_xzdistance']= (abs_mo_s1_x+abs_mo_s1_z)**0.5

        abs_mo_s2_x=abs(Three_vac['vacancy2_location_x'][i]-Three_vac['vacancy3_location_x'][i])**2
        abs_mo_s2_z=abs(Three_vac['vacancy2_location_z'][i]-Three_vac['vacancy3_location_z'][i])**2

        Three_vac.loc[i,'mo_s2_xzdistance']= (abs_mo_s2_x+abs_mo_s2_z)**0.5

        abs_s1_s2_x=abs(Three_vac['vacancy1_location_x'][i]-Three_vac['vacancy3_location_x'][i])**2
        abs_s1_s2_z=abs(Three_vac['vacancy1_location_z'][i]-Three_vac['vacancy3_location_z'][i])**2

        Three_vac.loc[i,'s1_s2_xzdistance']= (abs_s1_s2_x+abs_s1_s2_z)**0.5




    if Three_vac['vac3_element'][i]== "Mo":
        abs_mo_s1_x=abs(Three_vac['vacancy1_location_x'][i]-Three_vac['vacancy3_location_x'][i])**2
        abs_mo_s1_z=abs(Three_vac['vacancy1_location_z'][i]-Three_vac['vacancy3_location_z'][i])**2

        Three_vac.loc[i,'mo_s1_xzdistance']= (abs_mo_s1_x+abs_mo_s1_z)**0.5

        abs_mo_s2_x=abs(Three_vac['vacancy2_location_x'][i]-Three_vac['vacancy3_location_x'][i])**2
        abs_mo_s2_z=abs(Three_vac['vacancy2_location_z'][i]-Three_vac['vacancy3_location_z'][i])**2

        Three_vac.loc[i,'mo_s2_xzdistance']= (abs_mo_s2_x+abs_mo_s2_z)**0.5

        abs_s1_s2_x=abs(Three_vac['vacancy1_location_x'][i]-Three_vac['vacancy2_location_x'][i])**2
        abs_s1_s2_z=abs(Three_vac['vacancy1_location_z'][i]-Three_vac['vacancy2_location_z'][i])**2

        Three_vac.loc[i,'s1_s2_xzdistance']= (abs_s1_s2_x+abs_s1_s2_z)**0.5

  
  
  

Three_vac['mo_s1_yzdistance']=0
Three_vac['mo_s2_yzdistance']=0
Three_vac['s1_s2_yzdistance']=0

for i in range(0,len(Three_vac['mo_s1_yzdistance'])):
    if Three_vac['vac1_element'][i]== "Mo":
        abs_mo_s1_y=(Three_vac['vacancy1_location_y'][i]-Three_vac['vacancy2_location_y'][i])**2
        abs_mo_s1_z=(Three_vac['vacancy1_location_z'][i]-Three_vac['vacancy2_location_z'][i])**2

        Three_vac.loc[i,'mo_s1_yzdistance']= (abs_mo_s1_y+abs_mo_s1_z)**0.5

        abs_mo_s2_y=(Three_vac['vacancy1_location_y'][i]-Three_vac['vacancy3_location_y'][i])**2
        abs_mo_s2_z=(Three_vac['vacancy1_location_z'][i]-Three_vac['vacancy3_location_z'][i])**2

        Three_vac.loc[i,'mo_s2_yzdistance']= (abs_mo_s2_y+abs_mo_s2_z)**0.5


        abs_s1_s2_y=abs(Three_vac['vacancy2_location_y'][i]-Three_vac['vacancy3_location_y'][i])**2
        abs_s1_s2_z=abs(Three_vac['vacancy2_location_z'][i]-Three_vac['vacancy3_location_z'][i])**2

        Three_vac.loc[i,'s1_s2_yzdistance']= (abs_s1_s2_y+abs_s1_s2_z)**0.5
  
  
    if Three_vac['vac2_element'][i]== "Mo":
        abs_mo_s1_y=abs(Three_vac['vacancy1_location_y'][i]-Three_vac['vacancy2_location_y'][i])**2
        abs_mo_s1_z=abs(Three_vac['vacancy1_location_z'][i]-Three_vac['vacancy2_location_z'][i])**2

        Three_vac.loc[i,'mo_s1_yzdistance']= (abs_mo_s1_y+abs_mo_s1_z)**0.5

        abs_mo_s2_y=abs(Three_vac['vacancy2_location_y'][i]-Three_vac['vacancy3_location_y'][i])**2
        abs_mo_s2_z=abs(Three_vac['vacancy2_location_z'][i]-Three_vac['vacancy3_location_z'][i])**2

        Three_vac.loc[i,'mo_s2_yzdistance']= (abs_mo_s2_y+abs_mo_s2_z)**0.5

        abs_s1_s2_y=abs(Three_vac['vacancy1_location_y'][i]-Three_vac['vacancy3_location_y'][i])**2
        abs_s1_s2_z=abs(Three_vac['vacancy1_location_z'][i]-Three_vac['vacancy3_location_z'][i])**2

        Three_vac.loc[i,'s1_s2_yzdistance']= (abs_s1_s2_y+abs_s1_s2_z)**0.5




    if Three_vac['vac3_element'][i]== "Mo":
        abs_mo_s1_y=abs(Three_vac['vacancy1_location_y'][i]-Three_vac['vacancy3_location_y'][i])**2
        abs_mo_s1_z=abs(Three_vac['vacancy1_location_z'][i]-Three_vac['vacancy3_location_z'][i])**2

        Three_vac.loc[i,'mo_s1_yzdistance']= (abs_mo_s1_y+abs_mo_s1_z)**0.5

        abs_mo_s2_y=abs(Three_vac['vacancy2_location_y'][i]-Three_vac['vacancy3_location_y'][i])**2
        abs_mo_s2_z=abs(Three_vac['vacancy2_location_z'][i]-Three_vac['vacancy3_location_z'][i])**2

        Three_vac.loc[i,'mo_s2_yzdistance']= (abs_mo_s2_y+abs_mo_s2_z)**0.5

        abs_s1_s2_y=abs(Three_vac['vacancy1_location_y'][i]-Three_vac['vacancy2_location_y'][i])**2
        abs_s1_s2_z=abs(Three_vac['vacancy1_location_z'][i]-Three_vac['vacancy2_location_z'][i])**2

        Three_vac.loc[i,'s1_s2_yzdistance']= (abs_s1_s2_y+abs_s1_s2_z)**0.5

  
  
  

Three_vac['mo_angle']=0
def mo_angle(x1,y1,z1,x2,y2,z2,x3,y3,z3):
    num = (x2-x1)*(x3-x1)+(y2-y1)*(y3-y1)+(z2-z1)*(z3-z1)
    
    den_1 =((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)**0.5
    den_2=((x3-x1)**2+(y3-y1)**2+(z3-z1)**2)**0.5
    den=den_1*den_2
    
    angle = math.acos(num/den)
    
    return(angle)

for i in range(0, len(Three_vac['mo_angle'])):
    if Three_vac['vac1_element'][i]=="Mo":
        Three_vac.loc[i,'mo_angle']=mo_angle(Three_vac['vacancy1_location_x'][i], Three_vac['vacancy1_location_y'][i], Three_vac['vacancy1_location_z'][i],
                                       Three_vac['vacancy2_location_x'][i], Three_vac['vacancy2_location_y'][i],Three_vac['vacancy2_location_z'][i],
                                       Three_vac['vacancy3_location_x'][i], Three_vac['vacancy3_location_y'][i], Three_vac['vacancy3_location_z'][i])

  
    if Three_vac['vac2_element'][i]=="Mo":
        Three_vac[i,'mo_angle']=mo_angle(Three_vac['vacancy2_location_x'][i], Three_vac['vacancy2_location_y'][i], Three_vac['vacancy2_location_z'][i],
                                       Three_vac['vacancy1_location_x'][i], Three_vac['vacancy1_location_y'][i],Three_vac['vacancy1_location_z'][i],
                                       Three_vac['vacancy3_location_x'][i], Three_vac['vacancy3_location_y'][i], Three_vac['vacancy3_location_z'][i])

  

    if Three_vac['vac3_element'][i]=="Mo":
        Three_vac[i,'mo_angle']=mo_angle(Three_vac['vacancy3_location_x'][i], Three_vac['vacancy3_location_y'][i], Three_vac['vacancy3_location_z'][i],
                                       Three_vac['vacancy1_location_x'][i], Three_vac['vacancy1_location_y'][i],Three_vac['vacancy1_location_z'][i],
                                       Three_vac['vacancy2_location_x'][i], Three_vac['vacancy2_location_y'][i], Three_vac['vacancy2_location_z'][i])



Three_vac['mo_angle']=0
def mo_angle(x1,y1,z1,x2,y2,z2,x3,y3,z3):
    num = (x2-x1)*(x3-x1)+(y2-y1)*(y3-y1)+(z2-z1)*(z3-z1)
    
    den_1 =((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)**0.5
    den_2=((x3-x1)**2+(y3-y1)**2+(z3-z1)**2)**0.5
    den=den_1*den_2
    
    angle = math.acos(num/den)
    
    return(angle)

for i in range(0, len(Three_vac['mo_angle'])):
    if Three_vac['vac1_element'][i]=="Mo":
        Three_vac.loc[i,'mo_angle']=mo_angle(Three_vac['vacancy1_location_x'][i], Three_vac['vacancy1_location_y'][i], Three_vac['vacancy1_location_z'][i],
                                       Three_vac['vacancy2_location_x'][i], Three_vac['vacancy2_location_y'][i],Three_vac['vacancy2_location_z'][i],
                                       Three_vac['vacancy3_location_x'][i], Three_vac['vacancy3_location_y'][i], Three_vac['vacancy3_location_z'][i])

  
    if Three_vac['vac2_element'][i]=="Mo":
        Three_vac[i,'mo_angle']=mo_angle(Three_vac['vacancy2_location_x'][i], Three_vac['vacancy2_location_y'][i], Three_vac['vacancy2_location_z'][i],
                                       Three_vac['vacancy1_location_x'][i], Three_vac['vacancy1_location_y'][i],Three_vac['vacancy1_location_z'][i],
                                       Three_vac['vacancy3_location_x'][i], Three_vac['vacancy3_location_y'][i], Three_vac['vacancy3_location_z'][i])

  

    if Three_vac['vac3_element'][i]=="Mo":
        Three_vac[i,'mo_angle']=mo_angle(Three_vac['vacancy3_location_x'][i], Three_vac['vacancy3_location_y'][i], Three_vac['vacancy3_location_z'][i],
                                       Three_vac['vacancy1_location_x'][i], Three_vac['vacancy1_location_y'][i],Three_vac['vacancy1_location_z'][i],
                                       Three_vac['vacancy2_location_x'][i], Three_vac['vacancy2_location_y'][i], Three_vac['vacancy2_location_z'][i])



for i in range(0, len(Three_vac['vacancy2_location_x'])):
    a = np.array([Three_vac['vacancy2_location_x'][i],Three_vac['vacancy2_location_y'][i]])
    b = np.array([Three_vac['vacancy1_location_x'][i],Three_vac['vacancy1_location_y'][i]])
    c = np.array([Three_vac['vacancy3_location_x'][i],Three_vac['vacancy3_location_y'][i]])

    ba = a - b
    bc = c - b

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)

    Three_vac.loc[i,'mo_xy_angle']=(np.degrees(angle))

for i in range(0, len(Three_vac['vacancy2_location_x'])):
    b = np.array([Three_vac['vacancy2_location_x'][i],Three_vac['vacancy2_location_y'][i]])
    a = np.array([Three_vac['vacancy1_location_x'][i],Three_vac['vacancy1_location_y'][i]])
    c = np.array([Three_vac['vacancy3_location_x'][i],Three_vac['vacancy3_location_y'][i]])

    ba = a - b
    bc = c - b

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)

    Three_vac.loc[i,'s1_xy_angle']=(np.degrees(angle))

for i in range(0, len(Three_vac['vacancy2_location_x'])):
    centroid=((Three_vac['vacancy1_location_x'][i]+Three_vac['vacancy1_location_y'][i]+Three_vac['vacancy1_location_z'][i])/3,
             (Three_vac['vacancy2_location_x'][i]+Three_vac['vacancy2_location_y'][i]+Three_vac['vacancy2_location_z'][i])/3,
             (Three_vac['vacancy3_location_x'][i]+Three_vac['vacancy3_location_y'][i]+Three_vac['vacancy3_location_z'][i])/3)


    Three_vac.loc[i,'centroid_xy_distance']=centroid[0]**2+centroid[1]**2

Three_vac['mo_xy_angle'] = Three_vac['mo_xy_angle'].fillna(0)
Three_vac['s1_xy_angle'] = Three_vac['s1_xy_angle'].fillna(0)

Three_vac.to_csv("Three_vac_2-24.csv")

