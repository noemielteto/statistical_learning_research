# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import seaborn as sns
sns.set_context("poster")
from scipy import stats
from scipy.stats import exponnorm
from scipy.stats import mstats
import pandas as pd
import numpy as np
import collections
import os, errno
from pylab import *

def create_dummy_dataframe():
    '''
    This is for development purposes, to be able to try out new functions and features.
    '''
    df=pd.DataFrame(np.random.randn(6,4), columns=list('eventBCD'))
    return df

def input_formatting(datafile):

    df=pd.read_csv(datafile) # for future reference: pd.read_csv is recommended as opposed to pd.DataFrame.from_csv
    print df

    df=df.replace(r'\s+', np.nan, regex=True)
    df['Subject'] = pd.to_numeric(df['Subject'], errors='coerce')
    df['Block']=df['Block'].astype(basestring)
    df['epoch']=df['epoch'].astype(basestring)
    df['TripC']=df['TripC'].astype(basestring)
    df['firstRT'] = pd.to_numeric(df['firstRT'], errors='coerce')

    return df

def get_highs(pcode):
    table_=pd.read_csv('PCODE.csv', sep=',')
    high_triplets=[]
    pcode = pcode - 1 #index of the pcode, in order to refer to it in the pcode daframe (where the 0th row is the first pcode)
    high_triplets.append(table_['H1'][pcode])
    high_triplets.append(table_['H2'][pcode])
    high_triplets.append(table_['H3'][pcode])
    high_triplets.append(table_['H4'][pcode]) #get the four mappings for high freq triplets    

    return high_triplets    

def calc_overall_ACC(df):

    general_ACC = pd.pivot_table(df, values='firstACC', index='Subject', aggfunc=np.mean).rename(columns={'firstACC':'generalACC'})
    general_ACC['Subject'] = general_ACC.index
    df['Subject'] = pd.to_numeric(df['Subject'], errors='coerce')
    return general_ACC

def calc_overall_RT(df):

    general_RT = pd.pivot_table(df, values='firstRT', index='Subject', aggfunc=np.median).rename(columns={'firstRT':'generalRT'})
    print general_RT
    general_RT['Subject'] = general_RT.index
    df['Subject'] = pd.to_numeric(df['Subject'], errors='coerce')
    return general_RT

# ezt teszteld le!!!! adott szemely nyers es standard RT-jet plottold hisztogramra és nezd meg hogy egyezik-e a ketto
def standardize(df, method='z-score'):
    n=0

#TODO
#/home/noemi/Dropbox/Elteto_N/python_tools/learning_strategy_functions.py:57: SettingWithCopyWarning: 
#A value is trying to be set on a copy of a slice from a DataFrame.
#Try using .loc[row_indexer,col_indexer] = value instead

    df['standardized_RT'] = np.nan
    print('megvan az ures oszlop')
    subjects = df.Subject.unique()

    if method == 'z-score':

        for subject in subjects:
            n+=1

            x = np.asarray(df[df.Subject == subject].firstRT)
            
            # we check whether we have all datapoints
            if np.isnan(x)[np.isnan(x) == True].size > 0:
                # if we miss some, we leave them out, so that fitting does not crash
                print('WARNING! This subject had missing datapoints.')
                x = [n for n in x  if str(n) != 'nan']

            # get parameters of the subject
            m=np.mean(x)
            stdev=np.std(x)

            # compute standardized scores
            for index, row in df[df.Subject == subject].iterrows():
                raw_RT = row['firstRT']
                standardized_value = (raw_RT - m) / stdev
                df.loc[index, 'standardized_RT'] = standardized_value

            print "we are ready with " + str(n) + " subjects"        

    elif method == 'exgauss':

        for subject in subjects:
            n+=1

            x = np.asarray(df[df.Subject == subject].firstRT)
            
            # we check whether we have all datapoints
            if np.isnan(x)[np.isnan(x) == True].size > 0:
                # if we miss some, we leave them out, so that fitting does not crash
                print('WARNING! This subject had missing datapoints.')
                x = [n for n in x  if str(n) != 'nan']

            # we get the parameters of the subject
            shape_,loc_,scale_=exponnorm.fit(x)
            lambda_=(shape_*scale_)**(-1)
            tau_=1/lambda_

            # compute standardized scores
            for index, row in df[df.Subject == subject].iterrows():
                raw_RT = row['firstRT']
                standardized_value = (raw_RT - loc_) / scale_
                df.loc[index, 'standardized_RT'] = standardized_value

            print "we are ready with " + str(n) + " subjects"        

    return df

def remove_outlier_Q(df, col, threshold):
    subjects=df.Subject.unique().tolist()
    df_filt_all=pd.DataFrame()

    for s in subjects:
        df_g = df[ df['Subject'] == s ]
        df_g_filt = df_g[(df_g[col] < df_g[col].quantile((100-float(threshold))/100)) & (df_g[col] > df_g[col].quantile(float(threshold)/100))]
        df_filt_all = pd.concat([df_filt_all, df_g_filt], ignore_index=True)    

    return df_filt_all

def remove_outlier_absolute(df, col, lower, upper):
    df_filt_all = df[(df[col] > lower ) & (df[col] < upper)]    
    return df_filt_all

def select_correct_trials(df):
    # only accurate resps
    df = df[(df['firstACC']==1)]
    return df

def select_relevant_trials(df):
    '''
    Select only those trials which are not practice trials and not repetitions (AAA, ABB) or trills (ABA)
    '''
    # get the index of the event column so we can refer to this with iloc
    event_col = df.columns.get_loc("event") 

    # tag n-2 and n-1 repetitions
    repetition=[]
    for index, row in df.iterrows():
        if (df.iloc[index-1,event_col] == df.iloc[index,event_col]) or (df.iloc[index-2,event_col] == df.iloc[index,event_col]):
            repetition.append(1)
        else: repetition.append(0)

    df['repetition']=repetition

    # filter out repetitions and trills
    df = df[df['repetition']==0]

    return df

def invite_personal_data(df, path_to_personal_data, sep=","):
    df_personal = pd.read_csv(path_to_personal_data, sep=sep)

    # unify the datatype of the key
    df['Subject']=df['Subject'].astype(basestring)
    df_personal['Subject']=df_personal['Subject'].astype(basestring)

    df = df.merge(df_personal, on='Subject')

    return df

def calc_median_RTs(df, index=['Subject','epoch','TT'], values=['standardized_RT','firstRT']):
    # we merge the trialtype and triplet type data, so we will have PH, RH, RL trials    
    df['trialtype'] = df['TrialType'] + df['TT']

    medians = pd.pivot_table(df, values=values, index=index, aggfunc=np.median)
    medians = medians.rename(index=str, columns={'firstRT':'RT', "standardized_RT": "stand_RT"})

    # this is a great method to copy a level of multiindex to a column;
    medians['Subject'] = medians.index.get_level_values('Subject')
    medians['epoch'] = medians.index.get_level_values('epoch')
    medians[index[2]] = medians.index.get_level_values(index[2]) #trialtype or TT

    return medians

def calc_learning_scores(df_TT, df_trialtype, RT_value):
    df_learning_scores = pd.DataFrame()
    Subject=[]
    epoch=[]
    stat_learn=[]
    trip_learn=[]

    for index, subpart in df_TT.groupby(['Subject', 'epoch']):
        print index
        Subject.append(index[0])
        epoch.append(index[1])
        trip_learn.append(subpart[RT_value][1] - subpart[RT_value][0])
        # low - high

    for index, subpart in df_trialtype.groupby(['Subject', 'epoch']):
        stat_learn.append(subpart[RT_value][2] - subpart[RT_value][1])
        # randomlow - randomhigh

    df_learning_scores['Subject']=Subject
    df_learning_scores['Subject'] = pd.to_numeric(df_learning_scores['Subject'], errors='coerce')
    df_learning_scores['epoch']=epoch
    df_learning_scores['trip_learn']=trip_learn
    df_learning_scores['stat_learn']=stat_learn

    return df_learning_scores

def make_RT_dict(df, grouping_var):
# creates a {grouping variable:{subject:[list of RTs in that group]}} dict

    RTs={}
    for index, row in df.iterrows():
    
        #if first encounter with this agegroup
        if row[grouping_var] not in RTs:
            RTs[row[grouping_var]]=dict()

        #if first encounter with this subject in the agegroup
        if row['Subject'] not in RTs[row[grouping_var]]:
            RTs [row[grouping_var]] [row['Subject']] =list()

        RTs[row[grouping_var]][row['Subject']].append(row['firstRT'])

    return RTs

def make_params_df(RTs,all_high_low):

    params_df=pd.DataFrame()

    subjects=[]
    age_groups=[]
    locs=[]
    scales=[]
    ts=[]
    ps=[]
    for group,values in RTs.iteritems():

        for subject in values:

            x = np.asarray(RTs[group][subject])
            print 'Reaction time data of the participant - for fitting'                
            print x

            # we check whether we have all datapoints
            if np.isnan(x)[np.isnan(x) == True].size > 0:
                # if we miss some, we leave them out, so that fitting does not crash
                print('WARNING! This subject had missing datapoints.')
                x = [n for n in x  if str(n) != 'nan']

            # ExGaussian parameters for subject
            shape_,loc_,scale_=exponnorm.fit(x)
            lambda_=(shape_*scale_)**(-1)
            tau_=1/lambda_

            # test goodness of fit with Kolgomorov-Smirnov
            d,p = stats.kstest(x, 'exponnorm', args=(shape_,loc_,scale_))

            subjects.append(subject)
            age_groups.append(group)
            locs.append(loc_)
            scales.append(scale_)
            ts.append(tau_)
            ps.append(p)

    # add the results to the df

    params_df['Subject']            = np.asarray(subjects)
    params_df['age_group']            = np.asarray(age_groups)    
    params_df['mu_'+all_high_low]        = np.asarray(locs)
    params_df['sigma_'+all_high_low]    = np.asarray(scales)
    params_df['tau_'+all_high_low]        = np.asarray(ts)
    params_df['goodness_'+all_high_low]    = np.asarray(ps)

    return params_df

##############################################################

# calculate the variancess of the mu parameters within a schema (eg: 112, 122, 132, 142)
def sort_RTs_into_schemas(df,pcode):

    highs=get_highs(pcode)
    schema_RTs={highs[0]:{},highs[1]:{},highs[2]:{},highs[3]:{}}

    for index,row in df.iterrows():
        triplet=str(row['TripC'])
        schema = triplet[0] + '_' + triplet[2]

        # if we dont have this schema yet in the dictionary, we have to enter it
        if triplet not in schema_RTs[schema]:
            schema_RTs[schema][triplet]=[]

        schema_RTs[schema][triplet].append(row['firstRT'])

    return schema_RTs

def strategy_analysis(schema_RTs, pcode, method='divide_with_variance'):

    highs=get_highs(pcode)

    schema_medians={highs[0]:[],highs[1]:[],highs[2]:[],highs[3]:[]}
    all_locs=[]

    if method == 'zscore_triplets':

        for schema in schema_RTs.keys():
            for triplet in schema_RTs[schema].keys():

                l=np.asarray(schema_RTs[schema][triplet])
                med=np.median(l)
                all_locs.append(med) # here we simply gather the 12 triplet locs - for computing total variances
                
        total_std = np.std(all_locs)
        total_mean = np.mean(all_locs)
            
        for schema in schema_RTs.keys():
            for triplet in schema_RTs[schema].keys():

                l=np.asarray(schema_RTs[schema][triplet])
                med=np.median(l)
                stand = (med - total_mean)/total_std

                schema_medians[schema].append(stand)  # here we will have 3 loc values (eg. for 112, 132, 142) for every schema

    elif method == 'divide_with_variance':

        for schema in schema_RTs.keys():
            for triplet in schema_RTs[schema].keys():

                l=np.asarray(schema_RTs[schema][triplet])
                med=np.median(l)

                schema_medians[schema].append(med)  # here we will have 3 loc values (eg. for 112, 132, 142) for every schema
                all_locs.append(med)
        total_std = np.std(all_locs)
        total_mean = np.mean(all_locs)
    
    #####################################################
      
    within_schema_stds=[] # we calculate the 4 std values within the 4 schemas
    within_schema_means=[]
    for schema in schema_medians.keys(): # for every schema
        within_schema_stds.append(np.std(schema_medians[schema])) # we calculate the variances of the 4 loc values and add it
        within_schema_means.append(np.mean(schema_medians[schema]))

    ####################################################
    
    if method=='divide_with_variance':
        local_SOD_coherence = total_std/np.mean(within_schema_stds)
        global_SOD_coherence = total_std/np.var(within_schema_means)

    elif method=='zscore_triplets':
        local_SOD_coherence = 1/np.mean(within_schema_stds)
        global_SOD_coherence = 1/np.var(within_schema_means)

    return schema_medians, local_SOD_coherence, global_SOD_coherence

def display_schema_heatmap(schema_medians,dataset,subj,cmap):

    df = pd.DataFrame.from_dict(schema_medians)

    f, ax = plt.subplots(figsize=(3,3))

    sns.heatmap(df, vmin=-2, vmax=2, cmap=cmap, annot=False, yticklabels=False, square=True, ax=ax)
    ax.set_xlabel('schemas')
    ax.set_ylabel('triplets')
    ax.set_title('Subject: ' + str(subj))

    f.savefig('heatmap_' + str(dataset) + str(subj) + '.png')
    f.clf()

"""
########################################################################################################################
between_schema_stds=[] # we calculate the 4 variances values within the 3 rows grouped by 2nd element (thus the 4 triplets come from 4 different schemas)
between_schema_means=[]

# minden mindennel variaciot is visszahozni

for i in range(3):
    between_set = [ schema_medians[highs[0]][i],
                    schema_medians[highs[1]][i],
                    schema_medians[highs[2]][i],
                    schema_medians[highs[3]][i]]

    between_schema_stds.append(np.std(between_set))
    between_schema_means.append(np.mean(between_set))

"""
