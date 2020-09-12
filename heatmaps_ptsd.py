#%% Housekeeping
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import numpy as np
import os
import csv
#%% Function def

def substr_after(s, delim): # partition finds the matchand splits the string into a tuple of [0/1/2] = before/match/after, s = string, delim = 
    return s.partition(delim)[2] # return everything afterwards

#%% Import t-scores and p-values 
root_dir = os.path.join(os.getcwd(),'analyses')
# dictionary for converting string iterable to integer {ENTER DEP VARS OF YOUR CHOICE HERE}
vars_dict = {"dmri_dti.full.fa.gm_cort.desikan" : 0, "dmri_dti.full.fa.wm_cort.desikan" : 1, "dmri_dti.full.fa.gwc_cort.desikan" : 2,
             "dmri_dti.full.md.gm_cort.desikan" : 3, "dmri_dti.full.md.wm_cort.desikan" : 4, "dmri_dti.full.md.gwc_cort.desikan" : 5, 
             "dmri_dti.full.ld.gm_cort.desikan" : 6, "dmri_dti.full.ld.wm_cort.desikan" : 7, "dmri_dti.full.ld.gwc_cort.desikan" : 8, 
             "dmri_rsi.nd.gm_cort.desikan" : 9, "dmri_rsi.nd.wm_cort.desikan" : 10, "dmri_rsi.nd.gwc_cort.desikan" : 11,  
             "dmri_rsi.n0.gm_cort.desikan" : 12, "dmri_rsi.n0.wm_cort.desikan" : 13, "dmri_rsi.n0.gwc_cort.desikan" : 14,             
             "smri_thick_cort.desikan" : 15, "smri_area_cort.desikan" : 16, "smri_vol_cort.desikan" : 17}

# initialize dataframes {MAKE SURE THESE MATCH THE # OF VARS ABOVE}
tvalues_one = np.zeros((71,18), dtype=float)
tvalues_two = np.zeros((71,18), dtype=float)
pvalues_one = np.zeros((71,18), dtype=float)
pvalues_two = np.zeros((71,18), dtype=float)

# iterate through folders and build dataframe
for file in os.listdir(os.path.join(root_dir, "ptsd_categorized", "exposures", "tables", "table_coefs")): # {"ptsd_categorized is what I named my custom output folder"}
    filename = os.fsdecode(file)
    for x in vars_dict:
        if filename.startswith(x):
            df = pd.read_csv(os.path.join(root_dir, "ptsd_categorized", "exposures", "tables", "table_coefs", filename), usecols=["parameter_comp","t value","Pr(>|t|)"])
            df_ones = df[df['parameter_comp'].str.match('exposuresone')] # slice rows matching strings
            df_two = df[df['parameter_comp'].str.match('exposurestwo+')]
            
            temp = df_ones.loc[:,"t value"] # need to slice x:x+1 rather than x
            tvalues_one[:,vars_dict[x]:(vars_dict[x] + 1)] = temp[:, np.newaxis] # need to create new axis to correctly broadcast 71 to 71,1
            temp = df_two.loc[:,"t value"]
            tvalues_two[:,vars_dict[x]:(vars_dict[x] + 1)] = temp[:, np.newaxis] 
            temp = df_ones.loc[:,"Pr(>|t|)"]
            pvalues_one[:,vars_dict[x]:(vars_dict[x] + 1)] = temp[:, np.newaxis] 
            temp = df_two.loc[:,"Pr(>|t|)"]
            pvalues_two[:,vars_dict[x]:(vars_dict[x] + 1)] = temp[:, np.newaxis] 
            print(filename)
            
# read whatever the last .csv was to grab ordered Desikan atlas parcellations
desikan_parc = pd.read_csv(os.path.join(root_dir, "ptsd_categorized", "exposures", "tables", "table_coefs", filename), usecols=["parameter_comp", "dep_var"]) 
desikan_parc = desikan_parc.loc[desikan_parc['parameter_comp'] == 'exposuresone'].dropna() # extract parcellations from just one variable to avoid dups
desikan_parc = desikan_parc["dep_var"].values # convert to np array
# partition ordered desikan parcellation names
i = 0 
for value in desikan_parc:
    desikan_parc[i] = substr_after(desikan_parc[i], "desikan_")
    i = i + 1

roi_one = np.zeros((71,15), dtype=float)
roi_two = np.zeros((71,15), dtype=float)

# Generate heat maps
cmap = sns.diverging_palette(220, 20, sep=20, as_cmap=True)

xlabels = ["FA (GM)", "FA (WM)", "FA (GWC)", # {COLUMN NAMES FOR HEATMAP, MATCH DEP VARS}
           "MD (GM)", "MD (WM)", "MD (GWC)", 
           "LD (GM)", "LD (WM)", "LD (GWC)", 
           "ND (GM)", "ND (WM)", "ND (GWC)", 
           "N0 (GM)", "N0 (WM)", "N0 (GWC)", 
           "Thick", "Area", "Volume"]

fig, [ax1, ax2] = plt.subplots(1, 2, sharey=True, figsize=(15,10))
fig.tight_layout()
ax3 = fig.add_axes([1, 0.0275, 0.05, 0.95]) # add new axes for cbar [horizontal, vertical, width, height]
ax1.tick_params(labelsize=14, size=0, rotation=45)
ax2.tick_params(labelsize=14, size=0, rotation=45)
ax3.tick_params(labelsize=16, size=0, rotation=270)
spec = gridspec.GridSpec(ncols=2, nrows=1, figure=fig)
sns.heatmap(tvalues_one, cmap=cmap, vmin=-4, vmax=4, cbar=False, ax=ax1, xticklabels=xlabels, yticklabels=False)
sns.heatmap(tvalues_two, cmap=cmap, vmin=-4, vmax=4, ax=ax2, xticklabels=xlabels, yticklabels=False, cbar=True, cbar_ax=ax3)

save_dir = os.path.join(os.getcwd(),'plots', 'ptsd', 'ptsd_heatmaps.pdf') # {CAN ENTER CUSTOM SAVE DIR HERE}
fig.savefig(save_dir, bbox_inches='tight')

# Generate masked heat maps
fig, [ax1, ax2] = plt.subplots(1, 2, sharey=True, figsize=(15,10))
fig.tight_layout()
ax3 = fig.add_axes([1, 0.0275, 0.05, 0.95]) # add new axes for cbar [horizontal, vertical, width, height]
ax1.tick_params(labelsize=14, size=0, rotation=45)
ax2.tick_params(labelsize=14, size=0, rotation=45)
ax3.tick_params(labelsize=16, size=0, rotation=270)
sns.heatmap(tvalues_one, cmap=cmap, vmin=-4, vmax=4, cbar=False, ax=ax1, xticklabels=xlabels, yticklabels=False, mask=(pvalues_one>0.05))
sns.heatmap(tvalues_two, cmap=cmap, vmin=-4, vmax=4, ax=ax2, xticklabels=xlabels, yticklabels=False, 
            mask=(pvalues_two>0.05), cbar=True, cbar_ax=ax3)

save_dir = os.path.join(os.getcwd(),'plots', 'ptsd', 'ptsd_heatmaps_masked.pdf')
fig.savefig(save_dir, bbox_inches='tight')

#%% 2D heatmaps for DTI atlas-based measures
root_dir = os.path.join(os.getcwd(),'analyses')

# dictionary for converting string iterable to integer
vars_dict = {"dmri_rsi.vol_fiber" : 0, "dmri_dti.fa_fiber.at" : 1, "dmri_dti.md_fiber.at" : 2, 
             "dmri_dti.ld_fiber.at" : 3, "dmri_dti.td_fiber.at" : 4, "dmri_rsi.n0_fiber.at" : 5, "dmri_rsi.nd_fiber.at": 6}

# initialize dataframes
# df_total = pd.DataFrame(data=None, columns=dep_vars)
# df_severity = pd.DataFrame(data=None, columns=dep_vars)
tvalues_one = np.zeros((42,7), dtype=float)
tvalues_two = np.zeros((42,7), dtype=float)
pvalues_one = np.zeros((42,7), dtype=float)
pvalues_two = np.zeros((42,7), dtype=float)

# iterate through folders and build dataframe
for file in os.listdir(os.path.join(root_dir, "ptsd_categorized", "exposures", "tables", "table_coefs")): 
    filename = os.fsdecode(file)
    for x in vars_dict:
        if filename.startswith(x):
            df = pd.read_csv(os.path.join(root_dir, "ptsd_categorized", "exposures", "tables", "table_coefs", filename), usecols=["parameter_comp","t value","Pr(>|t|)"])
            df_ones = df[df['parameter_comp'].str.match('exposuresone')] # slice rows matching strings
            df_two = df[df['parameter_comp'].str.match('exposurestwo+')]
            
            temp = df_ones.loc[:,"t value"] # need to slice x:x+1 rather than x
            tvalues_one[:,vars_dict[x]:(vars_dict[x] + 1)] = temp[:, np.newaxis] # need to create new axis to correctly broadcast 71 to 71,1
            temp = df_two.loc[:,"t value"]
            tvalues_two[:,vars_dict[x]:(vars_dict[x] + 1)] = temp[:, np.newaxis] 
            temp = df_ones.loc[:,"Pr(>|t|)"]
            pvalues_one[:,vars_dict[x]:(vars_dict[x] + 1)] = temp[:, np.newaxis] 
            temp = df_two.loc[:,"Pr(>|t|)"]
            pvalues_two[:,vars_dict[x]:(vars_dict[x] + 1)] = temp[:, np.newaxis] 
            print(filename)
        
    
# Generate heat maps
cmap = sns.diverging_palette(220, 20, sep=20, as_cmap=True)

xlabels = ["Vol", "FA", "MD", "LD", "TD", "N0", "ND"]
ylabels = list(desikan_parc)

fig, [ax1, ax2] = plt.subplots(1, 2, sharey=True, figsize=(15,10))
fig.tight_layout()
ax3 = fig.add_axes([1, 0.0275, 0.05, 0.95]) # add new axes for cbar [horizontal, vertical, width, height]
ax1.tick_params(labelsize=14, size=0, rotation=45)
ax2.tick_params(labelsize=14, size=0, rotation=45)
ax3.tick_params(labelsize=16, size=0, rotation=270)
spec = gridspec.GridSpec(ncols=2, nrows=1, figure=fig)
sns.heatmap(tvalues_one, cmap=cmap, vmin=-4, vmax=4, cbar=False, ax=ax1, xticklabels=xlabels, yticklabels=False)
sns.heatmap(tvalues_two, cmap=cmap, vmin=-4, vmax=4, ax=ax2, xticklabels=xlabels, yticklabels=False, cbar=True, cbar_ax=ax3)

save_dir = os.path.join(os.getcwd(),'plots', 'ptsd', 'ptsd_heatmaps_at.pdf')
fig.savefig(save_dir, bbox_inches='tight')

# Generate masked heat maps
fig, [ax1, ax2] = plt.subplots(1, 2, sharey=True, figsize=(15,10))
fig.tight_layout()
ax3 = fig.add_axes([1, 0.0275, 0.05, 0.95]) # add new axes for cbar [horizontal, vertical, width, height]
ax1.tick_params(labelsize=14, size=0, rotation=45)
ax2.tick_params(labelsize=14, size=0, rotation=45)
ax3.tick_params(labelsize=16, size=0, rotation=270)
sns.heatmap(tvalues_one, cmap=cmap, vmin=-4, vmax=4, cbar=False, ax=ax1, xticklabels=xlabels, yticklabels=False, mask=(pvalues_one>0.05))
sns.heatmap(tvalues_two, cmap=cmap, vmin=-4, vmax=4, ax=ax2, xticklabels=xlabels, yticklabels=False, 
            mask=(pvalues_two>0.05), cbar=True, cbar_ax=ax3)

save_dir = os.path.join(os.getcwd(),'plots', 'ptsd', 'ptsd_heatmaps_at_masked.pdf')
fig.savefig(save_dir, bbox_inches='tight')


#%% 2D heatmaps for ASEG

root_dir = os.path.join(os.getcwd(),'analyses')

# dictionary for converting string iterable to integer
vars_dict = {"smri_vol_subcort.aseg" : 0, "dmri_dti.fa_subcort.aseg" : 1, "dmri_dti.md_subcort.aseg" : 2, 
             "dmri_dti.ld_subcort.aseg" : 3, "dmri_dti.td_subcort.aseg" : 4, "dmri_rsi.n0_subcort.aseg" : 5, "dmri_rsi.nd_subcort.aseg": 6}

# initialize dataframes
# df_total = pd.DataFrame(data=None, columns=dep_vars)
# df_severity = pd.DataFrame(data=None, columns=dep_vars)
tvalues_one = np.zeros((30,7), dtype=float)
tvalues_two = np.zeros((30,7), dtype=float)
pvalues_one = np.zeros((30,7), dtype=float)
pvalues_two = np.zeros((30,7), dtype=float)

# iterate through folders and build dataframe
for file in os.listdir(os.path.join(root_dir, "ptsd_categorized", "exposures", "tables", "table_coefs")): 
    filename = os.fsdecode(file)
    for x in vars_dict:
        if filename.startswith(x):
            df = pd.read_csv(os.path.join(root_dir, "ptsd_categorized", "exposures", "tables", "table_coefs", filename), usecols=["parameter_comp","t value","Pr(>|t|)"])
            df_ones = df[df['parameter_comp'].str.match('exposuresone')] # slice rows matching strings
            df_two = df[df['parameter_comp'].str.match('exposurestwo+')]
            
            temp = df_ones.loc[:,"t value"] # need to slice x:x+1 rather than x
            tvalues_one[:,vars_dict[x]:(vars_dict[x] + 1)] = temp[:, np.newaxis] # need to create new axis to correctly broadcast 71 to 71,1
            temp = df_two.loc[:,"t value"]
            tvalues_two[:,vars_dict[x]:(vars_dict[x] + 1)] = temp[:, np.newaxis] 
            temp = df_ones.loc[:,"Pr(>|t|)"]
            pvalues_one[:,vars_dict[x]:(vars_dict[x] + 1)] = temp[:, np.newaxis] 
            temp = df_two.loc[:,"Pr(>|t|)"]
            pvalues_two[:,vars_dict[x]:(vars_dict[x] + 1)] = temp[:, np.newaxis] 
            print(filename)
        
    
# Generate heat maps
cmap = sns.diverging_palette(220, 20, sep=20, as_cmap=True)

xlabels = ["Vol", "FA", "MD", "LD", "TD", "N0", "ND"]
ylabels = list(desikan_parc)

fig, [ax1, ax2] = plt.subplots(1, 2, sharey=True, figsize=(15,10))
fig.tight_layout()
ax3 = fig.add_axes([1, 0.0275, 0.05, 0.95]) # add new axes for cbar [horizontal, vertical, width, height]
ax1.tick_params(labelsize=14, size=0, rotation=45)
ax2.tick_params(labelsize=14, size=0, rotation=45)
ax3.tick_params(labelsize=16, size=0, rotation=270)
spec = gridspec.GridSpec(ncols=2, nrows=1, figure=fig)
sns.heatmap(tvalues_one, cmap=cmap, vmin=-4, vmax=4, cbar=False, ax=ax1, xticklabels=xlabels, yticklabels=False)
sns.heatmap(tvalues_two, cmap=cmap, vmin=-4, vmax=4, ax=ax2, xticklabels=xlabels, yticklabels=False, cbar=True, cbar_ax=ax3)

save_dir = os.path.join(os.getcwd(),'plots', 'ptsd', 'ptsd_heatmaps_ASEG.pdf')
fig.savefig(save_dir, bbox_inches='tight')

# Generate masked heat maps
fig, [ax1, ax2] = plt.subplots(1, 2, sharey=True, figsize=(15,10))
fig.tight_layout()
ax3 = fig.add_axes([1, 0.0275, 0.05, 0.95]) # add new axes for cbar [horizontal, vertical, width, height]
ax1.tick_params(labelsize=14, size=0, rotation=45)
ax2.tick_params(labelsize=14, size=0, rotation=45)
ax3.tick_params(labelsize=16, size=0, rotation=270)
sns.heatmap(tvalues_one, cmap=cmap, vmin=-4, vmax=4, cbar=False, ax=ax1, xticklabels=xlabels, yticklabels=False, mask=(pvalues_one>0.05))
sns.heatmap(tvalues_two, cmap=cmap, vmin=-4, vmax=4, ax=ax2, xticklabels=xlabels, yticklabels=False, 
            mask=(pvalues_two>0.05), cbar=True, cbar_ax=ax3)

save_dir = os.path.join(os.getcwd(),'plots', 'ptsd', 'ptsd_heatmaps_ASEG_masked.pdf')
fig.savefig(save_dir, bbox_inches='tight')


