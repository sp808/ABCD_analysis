#%% Housekeeping
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import numpy as np
import os
#%% Import t-scores and p-values (DESIKAN)
root_dir = os.path.join(os.getcwd(),'analyses')

ind_vars = {"psych_total" : "prodrom_psych_ss_number", "psych_severity" : "prodrom_psych_ss_severity_score"} # pair folders with variables

vars_dict = {"dmri_dti.full.fa.gm_cort.desikan" : 0, "dmri_dti.full.fa.wm_cort.desikan" : 1, "dmri_dti.full.fa.gwc_cort.desikan" : 2,
             "dmri_dti.full.md.gm_cort.desikan" : 3, "dmri_dti.full.md.wm_cort.desikan" : 4, "dmri_dti.full.md.gwc_cort.desikan" : 5, 
             "dmri_dti.full.ld.gm_cort.desikan" : 6, "dmri_dti.full.ld.wm_cort.desikan" : 7, "dmri_dti.full.ld.gwc_cort.desikan" : 8, 
             "dmri_rsi.nd.gm_cort.desikan" : 9, "dmri_rsi.nd.wm_cort.desikan" : 10, "dmri_rsi.nd.gwc_cort.desikan" : 11,  
             "dmri_rsi.n0.gm_cort.desikan" : 12, "dmri_rsi.n0.wm_cort.desikan" : 13, "dmri_rsi.n0.gwc_cort.desikan" : 14,             
             "smri_thick_cort.desikan" : 15, "smri_area_cort.desikan" : 16, "smri_vol_cort.desikan" : 17}

df_total = np.zeros((71,18), dtype=float)
df_severity = np.zeros((71,18), dtype=float)
total_pvalues = np.zeros((71,18), dtype=float)
severity_pvalues = np.zeros((71,18), dtype=float)

# iterate through folders and build dataframe
for names in ind_vars:
    for file in os.listdir(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs")): # iterate through a file list in dir
        filename = os.fsdecode(file)
        for x in vars_dict:
            if filename.startswith(x):
                if names == "psych_total":
                    df = pd.read_csv(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs", filename), usecols=["t value"])
                    df_total[:,vars_dict[x]:(vars_dict[x] + 1)] = df.values # need to slice x:x+1 rather than x
                    df = pd.read_csv(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs", filename), usecols=["Pr(>|t|)"])
                    total_pvalues[:,vars_dict[x]:(vars_dict[x] + 1)] = df.values 
                    print(filename)
                if names == "psych_severity":
                    df = pd.read_csv(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs", filename), usecols=["t value"])
                    df_severity[:,vars_dict[x]:(vars_dict[x] + 1)] = df.values
                    df = pd.read_csv(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs", filename), usecols=["Pr(>|t|)"])
                    severity_pvalues[:,vars_dict[x]:(vars_dict[x] + 1)] = df.values 
                    print(filename)

# Generate heat maps
cmap = sns.diverging_palette(220, 20, sep=20, as_cmap=True)

xlabels = ["FA (GM)", "FA (WM)", "FA (GWC)", 
           "MD (GM)", "MD (WM)", "MD (GWC)", 
           "LD (GM)", "LD (WM)", "LD (GWC)", 
           "ND (GM)", "ND (WM)", "ND (GWC)", 
           "N0 (GM)", "N0 (WM)", "N0 (GWC)", 
           "Thick", "Area", "Volume"]

fig, [ax1, ax2] = plt.subplots(1, 2, sharey=True, figsize=(20,15))
fig.tight_layout()
ax3 = fig.add_axes([1, 0.0275, 0.05, 0.95]) # add new axes for cbar [horizontal, vertical, width, height]
ax1.tick_params(labelsize=14, size=0, rotation=45)
ax2.tick_params(labelsize=14, size=0, rotation=45)
ax3.tick_params(labelsize=16, size=0, rotation=270)
spec = gridspec.GridSpec(ncols=2, nrows=1, figure=fig)
sns.heatmap(df_total, cmap=cmap, vmin=-4, vmax=4, cbar=False, ax=ax1, xticklabels=xlabels, yticklabels=False)
sns.heatmap(df_severity, cmap=cmap, vmin=-4, vmax=4, ax=ax2, xticklabels=xlabels, yticklabels=False, cbar=True, cbar_ax=ax3)

save_dir = os.path.join(os.getcwd(),'plots', 'psychosis', 'psychosis_heatmaps.pdf')
fig.savefig(save_dir, bbox_inches='tight')

# Generate masked heat maps
fig, [ax1, ax2] = plt.subplots(1, 2, sharey=True, figsize=(15,10))
fig.tight_layout()
ax3 = fig.add_axes([1, 0.0275, 0.05, 0.95]) # add new axes for cbar [horizontal, vertical, width, height]
ax1.tick_params(labelsize=14, size=0, rotation=45)
ax2.tick_params(labelsize=14, size=0, rotation=45)
ax3.tick_params(labelsize=16, size=0, rotation=270)
sns.heatmap(df_total, cmap=cmap, vmin=-4, vmax=4, cbar=False, ax=ax1, xticklabels=xlabels, yticklabels=False, mask=(total_pvalues>0.05))
sns.heatmap(df_severity, cmap=cmap, vmin=-4, vmax=4, ax=ax2, xticklabels=xlabels, yticklabels=False, 
            mask=(severity_pvalues>0.05), cbar=True, cbar_ax=ax3)

save_dir = os.path.join(os.getcwd(),'plots', 'psychosis', 'psychosis_heatmaps_masked.pdf')
fig.savefig(save_dir, bbox_inches='tight')

#%% Import t-scores and p-values (DTI ATLAS)
root_dir = os.path.join(os.getcwd(),'analyses')

ind_vars = {"psych_total" : "prodrom_psych_ss_number", "psych_severity" : "prodrom_psych_ss_severity_score"} # pair folders with variables

vars_dict = {"dmri_rsi.vol_fiber" : 0, "dmri_dti.fa_fiber.at" : 1, "dmri_dti.md_fiber.at" : 2, 
             "dmri_dti.ld_fiber.at" : 3, "dmri_dti.td_fiber.at" : 4, "dmri_rsi.n0_fiber.at" : 5, "dmri_rsi.nd_fiber.at": 6}

df_total = np.zeros((42,7), dtype=float)
df_severity = np.zeros((42,7), dtype=float)
total_pvalues = np.zeros((42,7), dtype=float)
severity_pvalues = np.zeros((42,7), dtype=float)

# iterate through folders and build dataframe
for names in ind_vars:
    for file in os.listdir(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs")): # iterate through a file list in dir
        filename = os.fsdecode(file)
        for x in vars_dict:
            if filename.startswith(x):
                if names == "psych_total":
                    df = pd.read_csv(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs", filename), usecols=["t value"])
                    df_total[:,vars_dict[x]:(vars_dict[x] + 1)] = df.values # need to slice x:x+1 rather than x
                    df = pd.read_csv(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs", filename), usecols=["Pr(>|t|)"])
                    total_pvalues[:,vars_dict[x]:(vars_dict[x] + 1)] = df.values 
                    print(filename)
                if names == "psych_severity":
                    df = pd.read_csv(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs", filename), usecols=["t value"])
                    df_severity[:,vars_dict[x]:(vars_dict[x] + 1)] = df.values
                    df = pd.read_csv(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs", filename), usecols=["Pr(>|t|)"])
                    severity_pvalues[:,vars_dict[x]:(vars_dict[x] + 1)] = df.values 
                    print(filename)

#Generate heat maps
cmap = sns.diverging_palette(220, 20, sep=20, as_cmap=True)

xlabels = ["Vol", "FA", "MD", "LD", "TD", "N0", "ND"]

fig, [ax1, ax2] = plt.subplots(1, 2, sharey=True, figsize=(20,15))
fig.tight_layout()
ax3 = fig.add_axes([1, 0.0275, 0.05, 0.95]) # add new axes for cbar [horizontal, vertical, width, height]
ax1.tick_params(labelsize=14, size=0, rotation=45)
ax2.tick_params(labelsize=14, size=0, rotation=45)
ax3.tick_params(labelsize=16, size=0, rotation=270)
spec = gridspec.GridSpec(ncols=2, nrows=1, figure=fig)
sns.heatmap(df_total, cmap=cmap, vmin=-4, vmax=4, cbar=False, ax=ax1, xticklabels=xlabels, yticklabels=False)
sns.heatmap(df_severity, cmap=cmap, vmin=-4, vmax=4, ax=ax2, xticklabels=xlabels, yticklabels=False, cbar=True, cbar_ax=ax3)

save_dir = os.path.join(os.getcwd(),'plots', 'psychosis', 'psychosis_heatmaps_AT.pdf')
fig.savefig(save_dir, bbox_inches='tight')

#%Generate masked heat maps
fig, [ax1, ax2] = plt.subplots(1, 2, sharey=True, figsize=(15,10))
fig.tight_layout()
ax3 = fig.add_axes([1, 0.0275, 0.05, 0.95]) # add new axes for cbar [horizontal, vertical, width, height]
ax1.tick_params(labelsize=14, size=0, rotation=45)
ax2.tick_params(labelsize=14, size=0, rotation=45)
ax3.tick_params(labelsize=16, size=0, rotation=270)
sns.heatmap(df_total, cmap=cmap, vmin=-4, vmax=4, cbar=False, ax=ax1, xticklabels=xlabels, yticklabels=False, mask=(total_pvalues>0.05))
sns.heatmap(df_severity, cmap=cmap, vmin=-4, vmax=4, ax=ax2, xticklabels=xlabels, yticklabels=False, 
            mask=(severity_pvalues>0.05), cbar=True, cbar_ax=ax3)

save_dir = os.path.join(os.getcwd(),'plots', 'psychosis', 'psychosis_heatmaps_masked_AT.pdf')
fig.savefig(save_dir, bbox_inches='tight')

#%% Import t-scores and p-values (ASEG)
root_dir = os.path.join(os.getcwd(),'analyses')

ind_vars = {"psych_total" : "prodrom_psych_ss_number", "psych_severity" : "prodrom_psych_ss_severity_score"} # pair folders with variables

vars_dict = {"smri_vol_subcort.aseg" : 0, "dmri_dti.fa_subcort.aseg" : 1, "dmri_dti.md_subcort.aseg" : 2, 
             "dmri_dti.ld_subcort.aseg" : 3, "dmri_dti.td_subcort.aseg" : 4, "dmri_rsi.n0_subcort.aseg" : 5, "dmri_rsi.nd_subcort.aseg": 6}

df_total = np.zeros((30,7), dtype=float)
df_severity = np.zeros((30,7), dtype=float)
total_pvalues = np.zeros((30,7), dtype=float)
severity_pvalues = np.zeros((30,7), dtype=float)

# iterate through folders and build dataframe
for names in ind_vars:
    for file in os.listdir(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs")): # iterate through a file list in dir
        filename = os.fsdecode(file)
        for x in vars_dict:
            if filename.startswith(x):
                if names == "psych_total":
                    df = pd.read_csv(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs", filename), usecols=["t value"])
                    df_total[:,vars_dict[x]:(vars_dict[x] + 1)] = df.values # need to slice x:x+1 rather than x
                    df = pd.read_csv(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs", filename), usecols=["Pr(>|t|)"])
                    total_pvalues[:,vars_dict[x]:(vars_dict[x] + 1)] = df.values 
                    print(filename)
                if names == "psych_severity":
                    df = pd.read_csv(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs", filename), usecols=["t value"])
                    df_severity[:,vars_dict[x]:(vars_dict[x] + 1)] = df.values
                    df = pd.read_csv(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs", filename), usecols=["Pr(>|t|)"])
                    severity_pvalues[:,vars_dict[x]:(vars_dict[x] + 1)] = df.values 
                    print(filename)

# Generate heat maps
cmap = sns.diverging_palette(220, 20, sep=20, as_cmap=True)

xlabels = ["Vol", "FA", "MD", "LD", "TD", "N0", "ND"]

fig, [ax1, ax2] = plt.subplots(1, 2, sharey=True, figsize=(20,15))
fig.tight_layout()
ax3 = fig.add_axes([1, 0.0275, 0.05, 0.95]) # add new axes for cbar [horizontal, vertical, width, height]
ax1.tick_params(labelsize=14, size=0, rotation=45)
ax2.tick_params(labelsize=14, size=0, rotation=45)
ax3.tick_params(labelsize=16, size=0, rotation=270)
spec = gridspec.GridSpec(ncols=2, nrows=1, figure=fig)
sns.heatmap(df_total, cmap=cmap, vmin=-4, vmax=4, cbar=False, ax=ax1, xticklabels=xlabels, yticklabels=False)
sns.heatmap(df_severity, cmap=cmap, vmin=-4, vmax=4, ax=ax2, xticklabels=xlabels, yticklabels=False, cbar=True, cbar_ax=ax3)

save_dir = os.path.join(os.getcwd(),'plots', 'psychosis', 'psychosis_heatmaps_ASEG.pdf')
fig.savefig(save_dir, bbox_inches='tight')

# Generate masked heat maps
fig, [ax1, ax2] = plt.subplots(1, 2, sharey=True, figsize=(15,10))
fig.tight_layout()
ax3 = fig.add_axes([1, 0.0275, 0.05, 0.95]) # add new axes for cbar [horizontal, vertical, width, height]
ax1.tick_params(labelsize=14, size=0, rotation=45)
ax2.tick_params(labelsize=14, size=0, rotation=45)
ax3.tick_params(labelsize=16, size=0, rotation=270)
sns.heatmap(df_total, cmap=cmap, vmin=-4, vmax=4, cbar=False, ax=ax1, xticklabels=xlabels, yticklabels=False, mask=(total_pvalues>0.05))
sns.heatmap(df_severity, cmap=cmap, vmin=-4, vmax=4, ax=ax2, xticklabels=xlabels, yticklabels=False, 
            mask=(severity_pvalues>0.05), cbar=True, cbar_ax=ax3)

save_dir = os.path.join(os.getcwd(),'plots', 'psychosis', 'psychosis_heatmaps_masked_ASEG.pdf')
fig.savefig(save_dir, bbox_inches='tight')
