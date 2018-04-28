# -*- coding: utf-8 -*-

import sys
sys.path.insert(0, '/home/noemi/Dropbox/Elteto_N/python_tools')
from asrt_analysis_toolbox import *
##################################

print('What dataset are you working with:')
dataset = raw_input()

if dataset == 'devsci':
    df_formatted=input_formatting('dev data/dev-sci_ASRT_all_final.csv')
if dataset == 'bart':
    df_formatted=input_formatting('dev data/bart-asrt-dev_ASRT_all_final.csv')

#df_prefilt = remove_outlier_absolute(df_formatted, 'firstRT', 80, 3000)

############ Classic Learning Scores

df_selected = select_relevant_trials(df_formatted)

general_ACC          = calc_overall_ACC(df_selected)
df_selected_corrects = select_correct_trials(df_selected)

general_RT           = calc_overall_RT(df_selected_corrects)
general_perf         = pd.merge(general_ACC,general_RT, on='Subject')

median_TT_RTs         = calc_median_RTs(df_selected_corrects, index=['Subject','epoch','TT'],values=['firstRT'])
median_trialtype_RTs  = calc_median_RTs(df_selected_corrects, index=['Subject','epoch','trialtype'],values=['firstRT'])

learning_scores = calc_learning_scores(median_TT_RTs, median_trialtype_RTs, 'RT')
learning_scores = learning_scores.merge(general_perf)

# aggregate across epoch, so we will have session means
learning_scores = learning_scores.pivot_table(index=['Subject'],values=['trip_learn','stat_learn','generalACC','generalRT'])
learning_scores['Subject'] = learning_scores.index


#if dataset == 'bart':
    # TODO invite personal data

############## Learning scores with exgauss

df_selected_H = df_selected_corrects[(df_selected_corrects['TT']=='H')]
df_selected_L = df_selected_corrects[(df_selected_corrects['TT']=='L')]

RTs_all = make_RT_dict(df_selected_corrects, 'age_group')
RTs_H = make_RT_dict(df_selected_H, 'age_group')
RTs_L = make_RT_dict(df_selected_L, 'age_group')

strategy_scores = make_params_df(RTs_all,'all')
strategy_scores_H = make_params_df(RTs_H,'H')
strategy_scores_L = make_params_df(RTs_L,'L')

strategy_scores = strategy_scores.merge(strategy_scores_H)
strategy_scores = strategy_scores.merge(strategy_scores_L)

strategy_scores['mu_HLdiff']         = strategy_scores['mu_L']         - strategy_scores['mu_H']
strategy_scores['sigma_HLdiff']         = strategy_scores['sigma_L']       - strategy_scores['sigma_H']
strategy_scores['tau_HLdiff']         = strategy_scores['tau_L']     - strategy_scores['tau_H']

strategy_scores['mu_HLdiff_standard']     = strategy_scores['mu_HLdiff']    / strategy_scores['mu_all']
strategy_scores['sigma_HLdiff_standard']     = strategy_scores['sigma_HLdiff'] / strategy_scores['sigma_all']
strategy_scores['tau_HLdiff_standard']     = strategy_scores['tau_HLdiff']   / strategy_scores['tau_all']

df_selected_H = df_selected_corrects[(df_selected_corrects['TT']=='H')]
for subj in df_selected_H.Subject.unique():
    df = df_selected_H[(df_selected_H['Subject'] == subj)]
    pcode = int(df.PCode.iloc[0])

    print ('similarity analysis for participant: ' + str(subj))

    schema_RTs = sort_RTs_into_schemas(df,pcode)

    # WITHIN SCHEMA SIMILARITY

    average_within_schema_similarities = []
    for schema in schema_RTs.keys():
        within_schema_similarity=[]
        
        x=1 # triplet index to correlate with; we set it to 1 so one won't correlate with self
        for triplet in schema_RTs[schema].keys():

            for i in range(x,len(schema_RTs[schema].keys())):

                triplet_to_compare = schema_RTs[schema].keys()[i]

                # Kolmogorov Smirnov for the comparison of the two distributions
                k_value = ks_2samp(schema_RTs[schema][triplet] , schema_RTs[schema][triplet_to_compare]) [0] # we take only the k value, and disregard the p, which is just the ratio of k and sample size 
                
                within_schema_similarity.append(1-k_value)

            x+=1 # we increment it so we won't compare to those which we have already been compared (duplicates)
        
        average_within_schema_similarities.append(np.mean (within_schema_similarity) )

    within_NAD_similarity = np.mean(average_within_schema_similarities)

    # TOTAL SIMILARITY

    all_similarities = []
    triplet_RTs = sort_RTs_into_triplets(df,pcode)
    for triplet in triplet_RTs.keys():
            for triplet_to_compare in triplet_RTs.keys():
                k_value = ks_2samp(triplet_RTs[triplet] , triplet_RTs[triplet_to_compare]) [0] # we take only the k value, and disregard the p, which is just the ratio of k and sample size

                all_similarities.append(1-k_value)    

    # take out self-comparisons, which are values of 1:
    all_similarities=[x for x in all_similarities if x != 1]

    overall_similarity = np.mean(all_similarities)

    strategy_scores.loc[strategy_scores['Subject'] == subj, 'within_NAD_similarity'] = within_NAD_similarity
    strategy_scores.loc[strategy_scores['Subject'] == subj, 'overall_similarity'] = overall_similarity

    # TODO heatmap
    #triplet_RTs = sort_RTs_into_triplets(df,pcode)
    #triplets_df=pd.DataFrame.from_dict(triplet_RTs)

strategy_scores['NAD_learning'] = strategy_scores['within_NAD_similarity'] - strategy_scores['overall_similarity']

"""
############ Learning strategies

cmap=sns.diverging_palette(220, 20, n=30)

for subj in df_selected_H.Subject.unique():
    df = df_selected_H[(df_selected_H['Subject'] == subj)]
    pcode = int(df.PCode.iloc[0])

    print ('learning strategy for participant' + str(subj) + ':')

    schema_RTs = sort_RTs_into_schemas(df,pcode)

    schema_medians, local_SOD_coherence, global_SOD_coherence = strategy_analysis(schema_RTs, pcode, method='zscore_triplets')

    #display_schema_heatmap(schema_medians,dataset=dataset,subj=subj,cmap=cmap)
    
    strategy_scores.ix[strategy_scores.Subject==subj, 'local_SOD_coherence'] = local_SOD_coherence
    strategy_scores.ix[strategy_scores.Subject==subj, 'global_SOD_coherence'] = global_SOD_coherence

dev_learning_strategy_analysis.py:73: DeprecationWarning: 
.ix is deprecated. Please use
.loc for label based indexing or
.iloc for positional indexing

See the documentation here:
http://pandas.pydata.org/pandas-docs/stable/indexing.html#ix-indexer-is-deprecated
  strategy_scores.ix[strategy_scores.Subject==subj, 'total_std'] = total_std
"""

# Data merge

learning_and_strategy_scores = strategy_scores.merge(learning_scores, on='Subject')

if dataset == 'devsci':
    learning_and_strategy_scores.to_excel('learning_and_strategy_scores_devsci.xlsx')
    #learning_and_strategy_scores=pd.read_excel('learning_and_strategy_scores_devsci.xlsx')
    learning_and_strategy_scores = learning_and_strategy_scores[(learning_and_strategy_scores['generalACC']>0.8)]

if dataset == 'bart':
    #learning_and_strategy_scores.to_excel('learning_and_strategy_scores_bart.xlsx')
    #personal=pd.read_excel('dev data/BART-ASRT-DEV.xlsx')
    learning_and_strategy_scores=pd.read_excel('learning_and_strategy_scores_bart.xlsx')
    #learning_and_strategy_scores = pd.merge(learning_and_strategy_scores, personal, on='Subject')
    learning_and_strategy_scores['WM_composit'] = (learning_and_strategy_scores['ZDSPAN'] + learning_and_strategy_scores['ZCORSI'])/2

    # Filtering

    learning_and_strategy_scores = learning_and_strategy_scores[(learning_and_strategy_scores['generalACC']>0.8)]

    learning_and_strategy_scores = learning_and_strategy_scores[learning_and_strategy_scores["Subject"].isin([
     '267', # mérési hiba, adat nagy része elveszett
     '269', # erõsen figyelmetlen, osztályfõnök szerint figyelemzastd  
     '271', # megvágta az ujját vizsgálat elõtt, ezért nehezen és lassan nyomta a billenytûket 
     '220', # szülõ adatfelvétel közben elvitte, adat nagy része nem lett felvéve 
     '233', # > 1s RT,figyelmetlen
     '47', # volt ASRT-n 3 éve 
     '33', # depresszió 
     '72', # 3 óra alvás éjjel, szubjektív minõsége 1, megébredések száma 5 
    ]) == False]

################################### FIGURES

#learning_and_strategy_scores = learning_and_strategy_scores[(learning_and_strategy_scores['local_SOD_coherence']<4)]

#### the lineplots

#cmap = sns.dark_palette(color = "#61CC8E", n_colors = 11, reverse = True)
#cmap = sns.diverging_palette(120, 12, l=60, n=11, center="dark")
cmap = sns.color_palette("Set1", n_colors=11, desat=.5)

if dataset == 'devsci':
    age_order=['04-06', '07-08', '09-10', '11-12', '14-17', '18-29', '30-44', '45-59', '60-85']
    cmap = cmap[:5] + cmap[7:]
    
if dataset == 'bart':
    age_order=['4-6', '7-8', '9-10', '11-12', '13-14', '15-16', '17-18', '19-29', '30-44', '45-59', '60-85']
    cmap = cmap[1:8]
    cmap.insert(0,np.asarray([0,0,0]))
    cmap.append(np.asarray([0,0,0]))
    cmap.append(np.asarray([0,0,0]))
    cmap.append(np.asarray([0,0,0]))

sns.factorplot(x = 'age_group', y = 'trip_learn', data=learning_and_strategy_scores, order=age_order, color='black', size=9, aspect=2, fontsize=30)
sns.swarmplot(x = 'age_group', y = 'trip_learn', data=learning_and_strategy_scores, order=age_order, size=6, alpha=0.7, palette=cmap)
plt.xlabel(u'age group')
plt.ylabel(u'RT (ms)')
plt.ylim((-25,50))
plt.title(u'Overall statistical learning across development')
plt.savefig('trip_learn_' + dataset + '_.png')
plt.clf()

sns.lmplot(x = 'NAD_learning', y = 'trip_learn', data=learning_and_strategy_scores, order=age_order, hue='age_group', palette = cmap, fit_reg=False, size = 9)
#sns.regplot(x = 'NAD_learning', y = 'trip_learn', data=learning_and_strategy_scores, scatter_kws={'alpha':0}, line_kws={'color':'black'})
plt.xlabel(u'NAD learning')
plt.ylabel(u'overall statistical learning')
plt.ylim((-25,50))
plt.title(u'The relationship between\noverall statistical learning and NAD learning')
plt.savefig('trip_learn_NAD_learn_corr' + dataset + '_.png')
plt.clf()

"""
#### scatter plot

sns.lmplot(x="stat_learn", y="local_SOD_coherence", hue='age_group', palette=cmap, fit_reg=False, hue_order=age_order, data=learning_and_strategy_scores, size=7, legend_out=True)
plt.xlabel(u'predictable vs. less predictable difference (RT)')
plt.ylabel(u'second-order dependence learning\n(standardized schema coherence)')
plt.title(u'The relationship between\npredictability difference score and SOD learning')
plt.savefig('scatter_local_statlearn' + '.png')
plt.clf()

sns.lmplot(x="DSPAN", y="local_SOD_coherence", hue='age_group', palette=cmap, fit_reg=True, ci=None, hue_order=age_order, data=learning_and_strategy_scores, size=9, aspect=0.6, legend=False)
plt.xlabel(u'verbal working memory\n(digit span)')
plt.ylabel(u'second-order dependence learning \n (standardized schema coherence)')
plt.ylim((1,3.5))
plt.title(u'The relationship between\nverbal WM and SOD learning')
plt.savefig('scatter_local_vWM_new' + '.png')
plt.clf()

sns.lmplot(x="CORSI", y="local_SOD_coherence", hue='age_group', palette=cmap, fit_reg=True, ci=None, hue_order=age_order, data=learning_and_strategy_scores, size=9, aspect=0.6, legend=False)
plt.xlabel(u'spatial working memory\n(Corsi span)')
plt.ylabel(u' \n ')
plt.ylim((1,3.5))
plt.title(u'The relationship between\nspatial WM and SOD learning')
plt.savefig('scatter_local_sWM_new' + '.png')
plt.clf()

sns.lmplot(x="comp_letter", y="local_SOD_coherence", hue='age_group', palette=cmap, fit_reg=True, ci=None, hue_order=age_order, data=learning_and_strategy_scores, size=9, aspect=0.6, legend=False)
plt.xlabel(u'executive function\n(composit score of\ncounting span and letter fluency)')
plt.ylabel(u' \n  ')
plt.ylim((1,3.5))
plt.title(u'The relationship between\nEF and SOD learning')
plt.savefig('scatter_local_EF_new' + '.png')
plt.clf()

sns.lmplot(x="comp_fluency", y="local_SOD_coherence", hue='age_group', palette=cmap, fit_reg=True, ci=None, hue_order=age_order, data=learning_and_strategy_scores, size=7, aspect = 1.5, legend_out=True)
plt.xlabel(u'executive function (composit score)')
plt.ylabel(u' \n  ')
plt.title(u'The relationship between \n executive functions and SOD learning')
plt.savefig('legend' + '.png')
plt.clf()
"""

sns.lmplot(x = 'DSPAN', y = 'NAD_learning', data=learning_and_strategy_scores, scatter_kws={'color':'black'}, line_kws={'color':'black'}, size=5, aspect=1)
plt.xlabel(u'digit span')
plt.ylabel(u'NAD learning')
plt.title(u'Verbal short-term memory')
plt.savefig('NAD_learn_DSPAN_corr_' + dataset + '_.png')
plt.clf()
sns.lmplot(x = 'CORSI', y = 'NAD_learning', data=learning_and_strategy_scores, scatter_kws={'color':'black'}, line_kws={'color':'black'}, size=5, aspect=1)
plt.xlabel(u'Corsi span')
plt.ylabel(u'NAD learning')
plt.title(u'Spatial short-term memory')
plt.savefig('NAD_learn_CORSI_corr_' + dataset + '_.png')
plt.clf()
sns.lmplot(x = 'CSPAN', y = 'NAD_learning', data=learning_and_strategy_scores, scatter_kws={'color':'black'}, line_kws={'color':'black'}, size=5, aspect=1)
plt.xlabel(u'counting span')
plt.ylabel(u'NAD learning')
plt.title(u'Working memory')
plt.savefig('NAD_learn_EF_corr_' + dataset + '_.png')
plt.clf()



#######################################################################################################x
df=learning_and_strategy_scores
#age_order=age_order[1:8]

f, ax = plt.subplots(3, 7, figsize=(20,8), sharex=True, sharey='row')
#plt.xlim((-50,50))
#plt.ylim((1,3.5))
ax[0,0].scatter(df['NAD_learning'][df['age_group'] == age_order[0]], df['DSPAN'][df['age_group'] == age_order[0]], c=cmap[1])
ax[0,0].set_title(age_order[0])
ax[0,1].scatter(df['NAD_learning'][df['age_group'] == age_order[1]], df['DSPAN'][df['age_group'] == age_order[1]], c=cmap[2])
ax[0,1].set_title(age_order[1])
ax[0,2].scatter(df['NAD_learning'][df['age_group'] == age_order[2]], df['DSPAN'][df['age_group'] == age_order[2]], c=cmap[3])
ax[0,2].set_title(age_order[2])
ax[0,3].scatter(df['NAD_learning'][df['age_group'] == age_order[3]], df['DSPAN'][df['age_group'] == age_order[3]], c=cmap[4])
ax[0,3].set_title(age_order[3])
ax[0,4].scatter(df['NAD_learning'][df['age_group'] == age_order[4]], df['DSPAN'][df['age_group'] == age_order[4]], c=cmap[5])
ax[0,4].set_title(age_order[4])
ax[0,5].scatter(df['NAD_learning'][df['age_group'] == age_order[5]], df['DSPAN'][df['age_group'] == age_order[5]], c=cmap[6])
ax[0,5].set_title(age_order[5])
ax[0,6].scatter(df['NAD_learning'][df['age_group'] == age_order[6]], df['DSPAN'][df['age_group'] == age_order[6]], c=cmap[7])
ax[0,6].set_title(age_order[6])

ax[1,0].scatter(df['NAD_learning'][df['age_group'] == age_order[0]], df['CORSI'][df['age_group'] == age_order[0]], c=cmap[1])
ax[1,1].scatter(df['NAD_learning'][df['age_group'] == age_order[1]], df['CORSI'][df['age_group'] == age_order[1]], c=cmap[2])
ax[1,2].scatter(df['NAD_learning'][df['age_group'] == age_order[2]], df['CORSI'][df['age_group'] == age_order[2]], c=cmap[3])
ax[1,3].scatter(df['NAD_learning'][df['age_group'] == age_order[3]], df['CORSI'][df['age_group'] == age_order[3]], c=cmap[4])
ax[1,4].scatter(df['NAD_learning'][df['age_group'] == age_order[4]], df['CORSI'][df['age_group'] == age_order[4]], c=cmap[5])
ax[1,5].scatter(df['NAD_learning'][df['age_group'] == age_order[5]], df['CORSI'][df['age_group'] == age_order[5]], c=cmap[6])
ax[1,6].scatter(df['NAD_learning'][df['age_group'] == age_order[6]], df['CORSI'][df['age_group'] == age_order[6]], c=cmap[7])

ax[2,0].scatter(df['NAD_learning'][df['age_group'] == age_order[0]], df['CSPAN'][df['age_group'] == age_order[0]], c=cmap[1])
ax[2,1].scatter(df['NAD_learning'][df['age_group'] == age_order[1]], df['CSPAN'][df['age_group'] == age_order[1]], c=cmap[2])
ax[2,2].scatter(df['NAD_learning'][df['age_group'] == age_order[2]], df['CSPAN'][df['age_group'] == age_order[2]], c=cmap[3])
ax[2,3].scatter(df['NAD_learning'][df['age_group'] == age_order[3]], df['CSPAN'][df['age_group'] == age_order[3]], c=cmap[4])
ax[2,4].scatter(df['NAD_learning'][df['age_group'] == age_order[4]], df['CSPAN'][df['age_group'] == age_order[4]], c=cmap[5])
ax[2,5].scatter(df['NAD_learning'][df['age_group'] == age_order[5]], df['CSPAN'][df['age_group'] == age_order[5]], c=cmap[6])
ax[2,6].scatter(df['NAD_learning'][df['age_group'] == age_order[6]], df['CSPAN'][df['age_group'] == age_order[6]], c=cmap[7])

ax[2,3].set_xlabel(u'NAD learning')
ax[0,0].set_ylabel(u'digit span', rotation=90)
ax[1,0].set_ylabel(u'Corsi span', rotation=90)
ax[2,0].set_ylabel(u'counting span', rotation=90)
f.suptitle(u'The relationship between NAD learning and memory\nacross the development', y=1)
f.savefig('scatter_NAD_memory' + '.png')
plt.clf()


"""
# szimulacio
# TODO egy szemely exgaussijahoz viszonyítva csinalni, nem pedig whatever normal parameterekkel

a1=np.random.normal(400,100,size=62)
b1=np.random.normal(400,100,size=62)
c1=np.random.normal(400,100,size=62)

a2=np.random.normal(400,100,size=62)
b2=np.random.normal(400,100,size=62)
c2=np.random.normal(400,100,size=62)

a3=np.random.normal(400,100,size=62)
b3=np.random.normal(400,100,size=62)
c3=np.random.normal(400,100,size=62)

a4=np.random.normal(400,100,size=62)
b4=np.random.normal(400,100,size=62)
c4=np.random.normal(400,100,size=62)

triplet_RTs=[a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4]

triplet_medians=[]
for t in triplet_RTs:
    triplet_medians.append(np.median(t))

mean_=np.mean(triplet_medians)
std_=np.mean(triplet_medians)

#z-score
triplet_zscores=[]
for t in triplet_medians:
    triplet_zscores.append((t-mean_)/std_)

subset_coherences=[]
subset=[]
for z in triplet_zscores[:3]:
    subset.append(z)
subset_coherences.append(np.std(subset))

subset=[]
for z in triplet_zscores[3:6]:
    subset.append(z)
subset_coherences.append(np.std(subset))

subset=[]
for z in triplet_zscores[6:9]:
    subset.append(z)
subset_coherences.append(np.std(subset))

subset=[]
for z in triplet_zscores[9:12]:
    subset.append(z)
subset_coherences.append(np.std(subset))

"""

