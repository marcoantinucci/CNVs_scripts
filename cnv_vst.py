#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import itertools

data_del = pd.read_csv("./samples_vst_input_DEL.csv", header = 0)
data_dup = pd.read_csv("./samples_vst_input_DUP.csv", header = 0)

groups_SV_frequencies = pd.read_csv("./pops_SVs_frequencies.txt", header = 0, sep = "\t")

groups_SV_locations = pd.read_csv("./pops_SVs_locations.txt", header = 0, sep = "\t")

def vst(pop1, pop2, n_pop1, n_pop2):
    pop1_n = n_pop1
    pop2_n = n_pop2
    n_tot = pop1_n + pop2_n
    pop1_var = np.var(pop1)
    pop2_var = np.var(pop2)
    var_tot = np.var(np.r_[pop1, pop2])
    if var_tot == 0:
        vst = 0
    else:
        vst = 1 - (((pop1_var*pop1_n) + (pop2_var*pop2_n))/(n_tot*var_tot))
    return vst


groups = ["Roma", "Europe", "Middle_East", "South_Asia"]
group_pairs = [i for i in itertools.combinations(groups, 2)]

del_sv_tags = data_del["SV_tag"].unique().tolist()
dup_sv_tags = data_dup["SV_tag"].unique().tolist()

roma_n = len(data_del['ID'][(data_del['Group'] == 'Roma')].unique())
europe_n = len(data_del['ID'][(data_del['Group'] == 'Europe')].unique())
middleeast_n = len(data_del['ID'][(data_del['Group'] == 'Middle_East')].unique())
southasia_n = len(data_del['ID'][(data_del['Group'] == 'South_Asia')].unique())

grouppairs_vstXsv_dels = pd.DataFrame(columns = ['Group_pair', 'SV_id', 'Vst', 'Allele_f_DEL_group1', 'Allele_f_DEL_group2', 'Genomic_location'])
grouppairs_vstXsv_dups = pd.DataFrame(columns = ['Group_pair', 'SV_id', 'Vst', 'Allele_f_DUP_group1', 'Allele_f_DUP_group2', 'Genomic_location'])

del_index_iterator = 0
dup_index_iterator = 0

for i in group_pairs:
    pop1 = i[0]
    pop2 = i[1]
    pops = "-".join(i)
    if pop1 == "Roma":
        n_pop1 = roma_n
    if pop1 == "Europe":
        n_pop1 = europe_n
    if pop1 == "Middle_East":
        n_pop1 = middleeast_n
    if pop1 == "South_Asia":
        n_pop1 = southasia_n
    
    if pop2 == "Roma":
        n_pop2 = roma_n
    if pop2 == "Europe":
        n_pop2 = europe_n
    if pop2 == "Middle_East":
        n_pop2 = middleeast_n
    if pop2 == "South_Asia":
        n_pop2 = southasia_n
    
    for sv in del_sv_tags:
        pop1_go = True
        pop2_go = True
        pop1_df_ncopies_DEL = data_del['N_copies'][(data_del["Group"] == pop1) & (data_del["SV_tag"] == sv)].tolist()
        pop2_df_ncopies_DEL = data_del['N_copies'][(data_del["Group"] == pop2) & (data_del["SV_tag"] == sv)].tolist()
        
        uniq_cn_pop1 = set(pop1_df_ncopies_DEL)
        uniq_cn_pop2 = set(pop2_df_ncopies_DEL)
        l_uniq_cn_pop1 = len(uniq_cn_pop1)
        l_uniq_cn_pop2 = len(uniq_cn_pop2)
        if l_uniq_cn_pop1 == 1 and 2 in uniq_cn_pop1:
            pop1_go = False
        if l_uniq_cn_pop2 == 1 and 2 in uniq_cn_pop2:
            pop2_go = False
        
        if pop1_go == True and pop2_go == True:
            groups_current_sv = data_del['Group'][data_del['SV_tag'] == sv].unique()

            vst_currentSV = round(vst(pop1_df_ncopies_DEL, pop2_df_ncopies_DEL, n_pop1, n_pop2),6)
            grouppairs_vstXsv_dels.loc[del_index_iterator, 'Vst'] = vst_currentSV

            freq_sv_pop1 = float(groups_SV_frequencies.loc[(groups_SV_frequencies['Group'] == pop1) & 
                                                        (groups_SV_frequencies['SV_id'] == sv), 'f_a'].values[0])
            freq_sv_pop2 = float(groups_SV_frequencies.loc[(groups_SV_frequencies['Group'] == pop2) & 
                                                        (groups_SV_frequencies['SV_id'] == sv), 'f_a'].values[0])

            if pop1 in groups_current_sv and len(groups_current_sv) == 1:
                grouppairs_vstXsv_dels.loc[del_index_iterator, 'Allele_f_DEL_group1'] = "private in " + pop1
            else:
                grouppairs_vstXsv_dels.loc[del_index_iterator, 'Allele_f_DEL_group1'] = freq_sv_pop1

            if pop2 in groups_current_sv and len(groups_current_sv) == 1:
                grouppairs_vstXsv_dels.loc[del_index_iterator, 'Allele_f_DEL_group2'] = "private in " + pop2
            else:
                grouppairs_vstXsv_dels.loc[del_index_iterator, 'Allele_f_DEL_group2'] = freq_sv_pop2

            current_sv_location_pop1 = groups_SV_locations.loc[(groups_SV_locations['Group'] == pop1) & 
                                                               (groups_SV_locations['SV_id'] == sv), 'Genomic_Location'].values[0:]
            
            current_sv_location_pop2 = groups_SV_locations.loc[(groups_SV_locations['Group'] == pop2) & 
                                                               (groups_SV_locations['SV_id'] == sv), 'Genomic_Location'].values[0:]
            
            current_sv_location_pop_1_2 = list(set(current_sv_location_pop1) & set(current_sv_location_pop2))
            len_sv_loc_pop1_2 = len(current_sv_location_pop_1_2)
            if len_sv_loc_pop1_2 >= 2:
                for index, genoloc in enumerate(current_sv_location_pop_1_2):
                    if index == 0:
                        grouppairs_vstXsv_dels.loc[del_index_iterator, 'Genomic_location'] = genoloc + "_"
                    elif index == (len_sv_loc_pop1_2 - 1):
                        grouppairs_vstXsv_dels.loc[del_index_iterator, 'Genomic_location'] += genoloc
                    else:
                        grouppairs_vstXsv_dels.loc[del_index_iterator, 'Genomic_location'] += genoloc + "_"
            else:
                for genoloc in current_sv_location_pop1:
                    grouppairs_vstXsv_dels.loc[del_index_iterator, 'Genomic_location'] = genoloc
            
            grouppairs_vstXsv_dels.loc[del_index_iterator, 'Group_pair'] = pops
            grouppairs_vstXsv_dels.loc[del_index_iterator, 'SV_id'] = sv

            del_index_iterator += 1
    
    for sv in dup_sv_tags:
        pop1_go = True
        pop2_go = True

        pop1_df_ncopies_DUP = data_dup['N_copies'][(data_dup["Group"] == pop1) & (data_dup["SV_tag"] == sv)].tolist()
        pop2_df_ncopies_DUP = data_dup['N_copies'][(data_dup["Group"] == pop2) & (data_dup["SV_tag"] == sv)].tolist()
        
        uniq_cn_pop1 = set(pop1_df_ncopies_DUP)
        uniq_cn_pop2 = set(pop2_df_ncopies_DUP)
        l_uniq_cn_pop1 = len(uniq_cn_pop1)
        l_uniq_cn_pop2 = len(uniq_cn_pop2)
        if l_uniq_cn_pop1 == 1 and 2 in uniq_cn_pop1:
            pop1_go = False
        if l_uniq_cn_pop2 == 1 and 2 in uniq_cn_pop2:
            pop2_go = False

        if pop1_go == True and pop2_go == True:
            groups_current_sv = data_dup['Group'][data_dup['SV_tag'] == sv].unique()

            vst_currentSV = round(vst(pop1_df_ncopies_DUP, pop2_df_ncopies_DUP, n_pop1, n_pop2),6)
            grouppairs_vstXsv_dups.loc[dup_index_iterator, 'Vst'] = vst_currentSV

            freq_sv_pop1 = float(groups_SV_frequencies.loc[(groups_SV_frequencies['Group'] == pop1) & 
                                                        (groups_SV_frequencies['SV_id'] == sv), 'f_a'].values[0])
            freq_sv_pop2 = float(groups_SV_frequencies.loc[(groups_SV_frequencies['Group'] == pop2) & 
                                                        (groups_SV_frequencies['SV_id'] == sv), 'f_a'].values[0])

            if pop1 in groups_current_sv and len(groups_current_sv) == 1:
                grouppairs_vstXsv_dups.loc[dup_index_iterator, 'Allele_f_DUP_group1'] = "private in " + pop1
            else:
                grouppairs_vstXsv_dups.loc[dup_index_iterator, 'Allele_f_DUP_group1'] = freq_sv_pop1

            if pop2 in groups_current_sv and len(groups_current_sv) == 1:
                grouppairs_vstXsv_dups.loc[dup_index_iterator, 'Allele_f_DUP_group2'] = "private in " + pop2
            else:
                grouppairs_vstXsv_dups.loc[dup_index_iterator, 'Allele_f_DUP_group2'] = freq_sv_pop2

            current_sv_location_pop1 = groups_SV_locations.loc[(groups_SV_locations['Group'] == pop1) & 
                                                               (groups_SV_locations['SV_id'] == sv), 'Genomic_Location'].values[0:]
            
            current_sv_location_pop2 = groups_SV_locations.loc[(groups_SV_locations['Group'] == pop2) & 
                                                               (groups_SV_locations['SV_id'] == sv), 'Genomic_Location'].values[0:]
            
            current_sv_location_pop_1_2 = list(set(current_sv_location_pop1) & set(current_sv_location_pop2))
            len_sv_loc_pop1_2 = len(current_sv_location_pop_1_2)
            if len_sv_loc_pop1_2 >= 2:
                for index, genoloc in enumerate(current_sv_location_pop_1_2):
                    if index == 0:
                        grouppairs_vstXsv_dups.loc[dup_index_iterator, 'Genomic_location'] = genoloc + "_"
                    elif index == (len_sv_loc_pop1_2 - 1):
                        grouppairs_vstXsv_dups.loc[dup_index_iterator, 'Genomic_location'] += genoloc
                    else:
                        grouppairs_vstXsv_dups.loc[dup_index_iterator, 'Genomic_location'] += genoloc + "_"
            else:
                for genoloc in current_sv_location_pop1:
                    grouppairs_vstXsv_dups.loc[dup_index_iterator, 'Genomic_location'] = genoloc
            
            grouppairs_vstXsv_dups.loc[dup_index_iterator, 'Group_pair'] = pops
            grouppairs_vstXsv_dups.loc[dup_index_iterator, 'SV_id'] = sv

            dup_index_iterator += 1

grouppairs_vstXsv_dels.to_csv("./group_pairs_VstXSV_DELs.txt", index = False, sep = "\t")
grouppairs_vstXsv_dups.to_csv("./group_pairs_VstXSV_DUPs.txt", index = False, sep = "\t")