#%% Housekeeping
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import numpy as np
import os
#%% Import t-scores and p-values
root_dir = os.path.join(os.getcwd(),'analyses')

folders = ["psych_total", "psych_severity"]

ind_vars = {"psych_total" : "prodrom_psych_ss_number", "psych_severity" : "prodrom_psych_ss_severity_score"}

dep_vars = ["dmri_dti.full.fa.gm_cort.desikan", "dmri_dti.full.fa.wm_cort.desikan", "dmri_rsi.nds2.gm_cort.desikan", 
             "dmri_rsi.nds2.wm_cort.desikan", "smri_t2w.gray02_cort.desikan", "smri_t2w.white02_cort.desikan", 
             "smri_vol_cort.desikan"]

vars_dict = {"dmri_dti.full.fa.gm_cort.desikan" : 0, "dmri_dti.full.fa.wm_cort.desikan" : 1, "dmri_rsi.nds2.gm_cort.desikan" : 2, 
             "dmri_rsi.nds2.wm_cort.desikan" : 3, "smri_t2w.gray02_cort.desikan" : 4, "smri_t2w.white02_cort.desikan" : 5, 
             "smri_vol_cort.desikan" : 6}

# initialize dataframes
# df_total = pd.DataFrame(data=None, columns=dep_vars)
# df_severity = pd.DataFrame(data=None, columns=dep_vars)
df_total = np.zeros((71,7), dtype=float)
df_severity = np.zeros((71,7), dtype=float)
total_pvalues = np.zeros((71,7), dtype=float)
severity_pvalues = np.zeros((71,7), dtype=float)

# iterate through folders and build dataframe
for names in folders:
    for file in os.listdir(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs")): # iterate through a file list in dir
        filename = os.fsdecode(file)
        for x in dep_vars:
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

#%% Generate heat maps

cmap = sns.diverging_palette(220, 20, sep=20, as_cmap=True)

xlabels = ["FA (GM)", "FA (WM)", "NODDI (GM)", "NODDI (WM)", "T2 (GM)", "T2 (WM)", "Volume"]

fig, [ax1, ax2] = plt.subplots(1, 2, sharey=True, figsize=(15,10))
fig.tight_layout()
ax3 = fig.add_axes([1, 0.0275, 0.05, 0.95]) # add new axes for cbar [horizontal, vertical, width, height]
ax1.tick_params(labelsize=14, size=0, rotation=45)
ax2.tick_params(labelsize=14, size=0, rotation=45)
ax3.tick_params(labelsize=16, size=0, rotation=270)
spec = gridspec.GridSpec(ncols=2, nrows=1, figure=fig)
sns.heatmap(df_total, cmap=cmap, vmin=-4, vmax=4, cbar=False, ax=ax1, xticklabels=xlabels, yticklabels=False)
sns.heatmap(df_severity, cmap=cmap, vmin=-4, vmax=4, ax=ax2, xticklabels=xlabels, yticklabels=False, cbar=True, cbar_ax=ax3)

save_dir = os.path.join(os.getcwd(),'plots', 'psychosis_heatmaps.pdf')
fig.savefig(save_dir, bbox_inches='tight')

#%% Generate masked heat maps

fig, [ax1, ax2] = plt.subplots(1, 2, sharey=True, figsize=(15,10))
fig.tight_layout()
ax3 = fig.add_axes([1, 0.0275, 0.05, 0.95]) # add new axes for cbar [horizontal, vertical, width, height]
ax1.tick_params(labelsize=14, size=0, rotation=45)
ax2.tick_params(labelsize=14, size=0, rotation=45)
ax3.tick_params(labelsize=16, size=0, rotation=270)
sns.heatmap(df_total, cmap=cmap, vmin=-4, vmax=4, cbar=False, ax=ax1, xticklabels=xlabels, yticklabels=False, mask=(total_pvalues>0.05))
sns.heatmap(df_severity, cmap=cmap, vmin=-4, vmax=4, ax=ax2, xticklabels=xlabels, yticklabels=False, 
            mask=(severity_pvalues>0.05), cbar=True, cbar_ax=ax3)

save_dir = os.path.join(os.getcwd(),'plots', 'psychosis_heatmaps_masked.pdf')
fig.savefig(save_dir, bbox_inches='tight')

#%% 2D histograms for FA WM vs. NODDI GM

sns.scatterplot(x=df_total[:,2], y=df_total[:,3])
            
