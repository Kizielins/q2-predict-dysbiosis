# -----------------------------------------------------------------------------
# Copyright (c) 2023, Bioinformatics at Małopolska Centre of Biotechnology
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

import numpy as np
import pandas as pd
import pickle


    
def list_intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

def calculate_index(table, pathways_stratified, pathways_unstratified):
    
    core = list(pd.read_csv("data/core_functions.txt", sep="\t").iloc[:,0])
    positive_pairs = list(pd.read_csv("data/positive_pairs.txt", sep="\t").iloc[:,0])
    contributions = list(pd.read_csv("data/Top_strat_pathway_contributions_in_healthy.txt", sep="\t").iloc[:,0])
    ml_model = pickle.load(open('data/rf_model_v1.sav', 'rb'))

    # Keep rows as samples, columns as taxonomical species names
    table_df = table
    pathways_stratified_df = pathways_stratified
    pathways_unstratified_df = pathways_unstratified


    sample_list = list_intersection(list(table_df.columns),list(pathways_unstratified_df.columns))
    sample_list = list_intersection(sample_list,list(pathways_stratified_df.columns))
    params_df = pd.DataFrame(columns=["Sample","Gupta_good","Gupta_bad","Frac_of_core_functions_among_all","Frac_of_core_functions_found","Species_found_together","Contributions_per_species"])
    sample_col = []
    gupta_good_col = []
    gupta_bad_col = []
    frac_core_fun_col = []
    frac_core_found_col = []
    spec_found_tog_col = []
    contributions_col = []
    contr_per_spec_col = []
    
    for sample in sample_list:
        sample_row = []
        try:
            sample_taxonomy = table_df[sample]
            sample_paths_strat = pathways_stratified_df[sample]
            sample_paths_unstrat = pathways_unstratified_df[sample]

            sample_paths_unstrat_cumsum = pd.DataFrame(data = list(sample_paths_unstrat), columns=['sample'], index=sample_paths_unstrat.index).sort_values(by='sample', ascending=True)
            sample_paths_unstrat_cumsum['cumsum'] = sample_paths_unstrat_cumsum['sample'].cumsum()
            sample_paths_unstrat_cumsum = sample_paths_unstrat_cumsum.loc[sample_paths_unstrat_cumsum['cumsum'] > 0]
            sample_paths_unstrat = sample_paths_unstrat_cumsum.loc[sample_paths_unstrat_cumsum['cumsum'] > 0.00001].drop('cumsum', axis=1)
            
            sample_paths_strat_cumsum = pd.DataFrame(data = list(sample_paths_strat), columns=['sample'], index=sample_paths_strat.index).sort_values(by='sample', ascending=True)
            sample_paths_strat_cumsum['cumsum'] = sample_paths_strat_cumsum['sample'].cumsum()
            sample_paths_strat_cumsum = sample_paths_strat_cumsum.loc[sample_paths_strat_cumsum['cumsum'] > 0]
            sample_paths_strat = sample_paths_strat_cumsum.loc[sample_paths_strat_cumsum['cumsum'] > 0.001].drop('cumsum', axis=1)
            
            
            sample_taxonomy_cumsum = pd.DataFrame(data = table_df[sample])
            sample_taxonomy_cumsum.columns=['sample']
            sample_taxonomy_cumsum = sample_taxonomy_cumsum.sort_values(by='sample', ascending=True)
            sample_taxonomy_cumsum['cumsum'] = sample_taxonomy_cumsum['sample'].cumsum()
            sample_taxonomy_cumsum = sample_taxonomy_cumsum.loc[sample_taxonomy_cumsum['cumsum'] != 0]
            sample_taxonomy_cumsum = sample_taxonomy_cumsum.loc[sample_taxonomy_cumsum['cumsum'] > 0.1].drop('cumsum', axis=1)
            sample_taxonomy_cumsum.columns = [sample]
            sample_taxonomy = sample_taxonomy_cumsum[sample]
            
            
            # Calculate index
            
            # Gupta species found
            gupta_good_species = ['Alistipes_senegalensis','Bacteroidales_bacterium_ph8','Bifidobacterium_adolescentis','Bifidobacterium_angulatum','Bifidobacterium_catenulatum','Lachnospiraceae_bacterium_8_1_57FAA','Sutterella_wadsworthensis']
            gupta_bad_species = ['Anaerotruncus_colihominis','Atopobium_parvulum','Bifidobacterium_dentium','Blautia_producta','candidate_division_TM7_single_cell_isolate_TM7c','Clostridiales_bacterium_1_7_47FAA','Clostridium_asparagiforme','Clostridium_bolteae','Clostridium_citroniae','Clostridium_clostridioforme','Clostridium_hathewayi','Clostridium_nexile','Clostridium_ramosum','Clostridium_symbiosum','Eggerthella_lenta','Erysipelotrichaceae_bacterium_2_2_44A','Flavonifractor_plautii','Fusobacterium_nucleatum','Gemella_morbillorum','Gemella_sanguinis','Granulicatella_adiacens','Holdemania_filiformis','Klebsiella_pneumoniae','Lachnospiraceae_bacterium_1_4_56FAA','Lachnospiraceae_bacterium_2_1_58FAA','Lachnospiraceae_bacterium_3_1_57FAA_CT1','Lachnospiraceae_bacterium_5_1_57FAA','Lachnospiraceae_bacterium_9_1_43BFAA','Lactobacillus_salivarius','Peptostreptococcus_stomatis','Ruminococcaceae_bacterium_D16','Ruminococcus_gnavus','Solobacterium_moorei','Streptococcus_anginosus','Streptococcus_australis','Streptococcus_gordonii','Streptococcus_infantis','Streptococcus_mitis_oralis_pneumoniae','Streptococcus_sanguinis','Streptococcus_vestibularis','Subdoligranulum_sp_4_3_54A2FAA','Subdoligranulum_variabile','Veillonella_atypica'] 
            
            gupta_good = 0
            gupta_bad = 0

            for c in gupta_good_species:
                if c in list(sample_taxonomy[sample_taxonomy > 0].index):
                    gupta_good += 1
                    
            for c in gupta_bad_species:
                if c in list(sample_taxonomy[sample_taxonomy > 0].index):
                    gupta_bad += 1
                    
            sample_row.append(sample)
            sample_row.append(gupta_good)
            sample_row.append(gupta_bad)
            
            # fraction of core functions found, fraction of core functions among all functions
            core_found = 0
            for c in core:
                if c in list(sample_paths_unstrat[sample_paths_unstrat > 0].index):
                    core_found += 1
            core_functions_fraction_in_all = core_found/len(sample_paths_unstrat[sample_paths_unstrat > 0].index)
            core_functions_found = core_found/len(core)
            ##
            if core_functions_fraction_in_all > 0:
                sample_row.append(core_functions_fraction_in_all)
            else:
                sample_row.append(0)
            if core_functions_found > 0:
                sample_row.append(core_functions_found)
            else:
                sample_row.append(0)
            
            # Common occurrence of species positively correlated in health in at least 2/3 studies (fraction)
            pairs_count = 0
            for pair in positive_pairs:
                if pair.split("+")[0] in list(sample_taxonomy.index) and pair.split("+")[1] in list(sample_taxonomy.index):
                    if sample_taxonomy.loc[pair.split("+")[0]] > 0 and sample_taxonomy.loc[pair.split("+")[1]] > 0:
                        pairs_count += 1
            pairs_frac = pairs_count/len(positive_pairs)
            sample_row.append(pairs_frac)
            
            # Found in healthy contributions to pathways
            ## NOT USED IN THE ACTUAL MODEL ATM
            
            sample_paths_strat_nozero = sample_paths_strat.loc[sample_paths_strat['sample'] != 0]
            contributions_count = 0
            
            for a in list(sample_paths_strat_nozero.index):
                if a in contributions:
                    contributions_count += 1
            #if contributions_count == 0:
            #    contributions_count = 0.0000001
            
            contributions_final = contributions_count/len(contributions)
            if contributions_final > 0:
                sample_row.append(contributions_final)
            else:
                sample_row.append(0)
            # Average number of contributions to all functions per species
            
            all_species = list(sample_taxonomy.index)
            contributions_count = []
            for a in all_species:
                tmp = sample_paths_strat[sample_paths_strat.index.str.contains(a)]
                tmp = tmp[tmp > 0]
                contributions_count.append(tmp.shape[0])
            pathway_contributions_per_species = np.mean(contributions_count)
            contributions_count_nonneg = []
            for b in contributions_count:
                if b > 0:
                    contributions_count_nonneg.append(b)
            pathway_contributions_per_species_nonneg = np.mean(contributions_count_nonneg)
            if pathway_contributions_per_species_nonneg > 0:
                sample_row.append(pathway_contributions_per_species_nonneg)
            else:
                sample_row.append(1)
            
            ## Merging all

            if len(sample_row) == 8:
                
                sample_col.append(sample_row[0])
                #print(sample_col)
                gupta_good_col.append(sample_row[1])
                #print(gupta_good_col)
                gupta_bad_col.append(sample_row[2])
                #print(gupta_bad_col)
                frac_core_fun_col.append(sample_row[3])
                #print(frac_core_fun_col)
                frac_core_found_col.append(sample_row[4])
                #print(frac_core_found_col)
                spec_found_tog_col.append(sample_row[5])
                #print(spec_found_tog_col)
                contributions_col.append(sample_row[6])
                #print(contr_per_spec_col)
                contr_per_spec_col.append(sample_row[7])
                #print(contr_per_spec_col)
            
        except:
            pass

    params_df['Sample'] = sample_col
    params_df["Gupta_good"] = gupta_good_col
    params_df["Gupta_bad"] = gupta_bad_col
    params_df["Frac_of_core_functions_among_all"] = frac_core_fun_col
    params_df["Frac_of_core_functions_found"] = frac_core_found_col
    params_df["Species_found_together"] = spec_found_tog_col
    params_df["Contributions_to_pathways"] = contributions_col
    params_df["Contributions_per_species"] = contr_per_spec_col

    preds = ml_model.predict_proba(params_df.iloc[:,1:].values)
    scores_pred = []
    for a in list(preds):
        scores_pred.append(a[1])
    
    
    scores_pred_df = pd.DataFrame()
    scores_pred_df['SampleID'] = list(params_df['Sample'])
    scores_pred_df['Score'] = scores_pred
    scores_pred_df = scores_pred_df.set_index('SampleID')
    
    score_df = scores_pred_df['Score']
    score_df.name = "Dysbiosis_score"


    return score_df