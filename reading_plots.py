#%% README

# Exploratory data analysis & visualization tool to output 2D heatmaps and tabulated lists of significantly affected ROIs across 
# a range of user-defined sMRI/dMRI DTI + RSI structural parameters, parcellated by Desikan/DTI fiber/ASEG atlases

# Each codeblock corresponds to each of those 3 atlases and requires specifying two independent continuous variables 
# and a list of structural parameters.

# Assumptions:
#     - Both independent variables are continuous
#     - .csv outputs from DEAPext.R are in a folder in root directory called 'analyses' and output to the same folder for each independent variable
#     - Heatmap files are saved in a folder in root directory called 'plots'/custom folder name
#     - Tabulated outputs are saved in a folder in root directory called 'plots'/custom folder name/'lists'

# Outputs (for each of the 3 atlases):
    

# Specify each custom input whereever there are curly brackets with capslock text

#%% Housekeeping
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import numpy as np
import os
import csv
from itertools import zip_longest

#%% Function def

def substr_after(s, delim): # partition finds the matchand splits the string into a tuple of [0/1/2] = before/match/after, s = string, delim = 
    return s.partition(delim)[2] # return everything afterwards

def sig_rois(tvalues, pvalues, columns, rois, filepath, filename): # inputs: columns = list of structural parameter names, rois = list of parcellations
    output = [[] for i in range(len(columns))] # initialize list of lists
    #
    for (i,row) in enumerate(tvalues):
        for (j,value) in enumerate(row):
            if pvalues[i,j] < 0.05:
                string = rois[i] + ', ' + str(round(tvalues[i,j],5)) + ', ' + str(round(pvalues[i,j],5))
                output[j].append(string)
        
    export_data = zip_longest(*output, fillvalue = '') # transpose to vertical columns
    with open(os.path.join(filepath, filename), 'w', encoding="ISO-8859-1", newline='') as csvfile:
      wr = csv.writer(csvfile)
      wr.writerow((columns))
      wr.writerows(export_data)
    
#%% Import t-scores and p-values (DESIKAN)
root_dir = os.path.join(os.getcwd(),'analyses') # {ENTER DIR OF ANALYSIS FOLDERS}

ind_vars = {"reading_hours" : "sports_activity_ss_read_hours_p", "reading_years" : "sports_activity_ss_read_years_p"} # {ENTER FOLDER NAME : VAR NAME}

vars_dict = {"dmri_dti.full.fa.gm_cort.desikan" : 0, "dmri_dti.full.fa.wm_cort.desikan" : 1, "dmri_dti.full.fa.gwc_cort.desikan" : 2, # {ENTER ANY NUMBER OF DEP VARS BELOW WITH INCREMENTAL INTEGERS}
             "dmri_dti.full.md.gm_cort.desikan" : 3, "dmri_dti.full.md.wm_cort.desikan" : 4, "dmri_dti.full.md.gwc_cort.desikan" : 5, 
             "dmri_dti.full.ld.gm_cort.desikan" : 6, "dmri_dti.full.ld.wm_cort.desikan" : 7, "dmri_dti.full.ld.gwc_cort.desikan" : 8, 
             "dmri_rsi.nd.gm_cort.desikan" : 9, "dmri_rsi.nd.wm_cort.desikan" : 10, "dmri_rsi.nd.gwc_cort.desikan" : 11,  
             "dmri_rsi.n0.gm_cort.desikan" : 12, "dmri_rsi.n0.wm_cort.desikan" : 13, "dmri_rsi.n0.gwc_cort.desikan" : 14,             
             "smri_thick_cort.desikan" : 15, "smri_area_cort.desikan" : 16, "smri_vol_cort.desikan" : 17}

xlabels = ["FA (GM)", "FA (WM)", "FA (GWC)", # {CHANGE THESE TO MATCH DEP VAR NAMES}
           "MD (GM)", "MD (WM)", "MD (GWC)", 
           "LD (GM)", "LD (WM)", "LD (GWC)", 
           "ND (GM)", "ND (WM)", "ND (GWC)", 
           "N0 (GM)", "N0 (WM)", "N0 (GWC)", 
           "Thick", "Area", "Volume"]

df_var1 = np.zeros((71,18), dtype=float) # {MAKE SURE # OF COLUMNS MATCH # OF DEP VARS BELOW}
df_var2 = np.zeros((71,18), dtype=float)
pvalues_var1 = np.zeros((71,18), dtype=float)
pvalues_var2 = np.zeros((71,18), dtype=float)
run_once = 0
# iterate through folders and build dataframe
for names in ind_vars:
    for file in os.listdir(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs")): # iterate through a file list in dir
        filename = os.fsdecode(file)
        for x in vars_dict:
            if "desikan" in filename and run_once == 0:
                roi_parc = pd.read_csv(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs", filename), usecols=["dep_var"])
                roi_parc = roi_parc["dep_var"].values # convert to np array
                run_once = 1 # only grab ROIs once
            if filename.startswith(x):
                if names == "reading_hours": # {CHANGE THIS TO FOLDER FOR FIRST IND VAR}
                    df = pd.read_csv(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs", filename), usecols=["t value"])
                    df_var1[:,vars_dict[x]:(vars_dict[x] + 1)] = df.values # need to slice x:x+1 rather than x
                    df = pd.read_csv(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs", filename), usecols=["Pr(>|t|)"])
                    pvalues_var1[:,vars_dict[x]:(vars_dict[x] + 1)] = df.values 
                    print(filename)
                if names == "reading_years": # {CHANGE THIS TO FOLDER FOR SECOND IND VAR}
                    df = pd.read_csv(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs", filename), usecols=["t value"])
                    df_var2[:,vars_dict[x]:(vars_dict[x] + 1)] = df.values
                    df = pd.read_csv(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs", filename), usecols=["Pr(>|t|)"])
                    pvalues_var2[:,vars_dict[x]:(vars_dict[x] + 1)] = df.values 
                    print(filename)
                
# Generate heat maps
cmap = sns.diverging_palette(220, 20, sep=20, as_cmap=True)

fig, [ax1, ax2] = plt.subplots(1, 2, sharey=True, figsize=(20,15))
fig.tight_layout()
ax3 = fig.add_axes([1, 0.0275, 0.05, 0.95]) # add new axes for cbar [horizontal, vertical, width, height]
ax1.tick_params(labelsize=14, size=0, rotation=45)
ax2.tick_params(labelsize=14, size=0, rotation=45)
ax3.tick_params(labelsize=16, size=0, rotation=270)
spec = gridspec.GridSpec(ncols=2, nrows=1, figure=fig)
sns.heatmap(df_var1, cmap=cmap, vmin=-4, vmax=4, cbar=False, ax=ax1, xticklabels=xlabels, yticklabels=False)
sns.heatmap(df_var2, cmap=cmap, vmin=-4, vmax=4, ax=ax2, xticklabels=xlabels, yticklabels=False, cbar=True, cbar_ax=ax3)

save_dir = os.path.join(os.getcwd(),'plots', 'reading', 'reading_heatmaps.pdf') # {ENTER SAVE DIR HERE}
fig.savefig(save_dir, bbox_inches='tight')

# Generate masked heat maps
fig, [ax1, ax2] = plt.subplots(1, 2, sharey=True, figsize=(15,10))
fig.tight_layout()
ax3 = fig.add_axes([1, 0.0275, 0.05, 0.95]) # add new axes for cbar [horizontal, vertical, width, height]
ax1.tick_params(labelsize=14, size=0, rotation=45)
ax2.tick_params(labelsize=14, size=0, rotation=45)
ax3.tick_params(labelsize=16, size=0, rotation=270)
sns.heatmap(df_var1, cmap=cmap, vmin=-4, vmax=4, cbar=False, ax=ax1, xticklabels=xlabels, yticklabels=False, mask=(pvalues_var1>0.05))
sns.heatmap(df_var2, cmap=cmap, vmin=-4, vmax=4, ax=ax2, xticklabels=xlabels, yticklabels=False, 
            mask=(pvalues_var2>0.05), cbar=True, cbar_ax=ax3)

save_dir = os.path.join(os.getcwd(),'plots', 'reading', 'reading_heatmaps_masked.pdf')  # {ENTER SAVE DIR HERE}
fig.savefig(save_dir, bbox_inches='tight')

# partition ordered parcellation names
i = 0 
for value in roi_parc:
    roi_parc[i] = substr_after(roi_parc[i], "desikan_")
    i = i + 1
# output .csv files with significant ROIs
save_dir = os.path.join(os.getcwd(),'plots', 'reading', 'lists')  # {ENTER SAVE DIR HERE}
output = sig_rois(df_var1, pvalues_var1, xlabels, roi_parc, save_dir, 'reading_hours.csv')
output = sig_rois(df_var2, pvalues_var2, xlabels, roi_parc, save_dir, 'reading_years.csv')

#%% Import t-scores and p-values (DTI ATLAS)
root_dir = os.path.join(os.getcwd(),'analyses') # {ENTER DIR OF ANALYSIS FOLDERS}

ind_vars = {"reading_hours" : "sports_activity_ss_read_hours_p", "reading_years" : "sports_activity_ss_read_years_p"} # {ENTER FOLDER NAME : VAR NAME}
# {ENTER ANY NUMBER OF DEP VARS BELOW WITH INCREMENTAL INTEGERS}
vars_dict = {"dmri_rsi.vol_fiber" : 0, "dmri_dti.fa_fiber.at" : 1, "dmri_dti.md_fiber.at" : 2, 
             "dmri_dti.ld_fiber.at" : 3, "dmri_dti.td_fiber.at" : 4, "dmri_rsi.n0_fiber.at" : 5, "dmri_rsi.nd_fiber.at": 6}
# {MAKE SURE # OF COLUMNS MATCH # OF DEP VARS BELOW}
df_var1 = np.zeros((42,7), dtype=float)
df_var2 = np.zeros((42,7), dtype=float)
pvalues_var1 = np.zeros((42,7), dtype=float)
pvalues_var2 = np.zeros((42,7), dtype=float)
run_once = 0
# iterate through folders and build dataframe
for names in ind_vars:
    for file in os.listdir(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs")): # iterate through a file list in dir
        filename = os.fsdecode(file)
        for x in vars_dict:
            if "fiber.at" in filename and run_once == 0:
                roi_parc = pd.read_csv(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs", filename), usecols=["dep_var"])
                roi_parc = roi_parc["dep_var"].values # convert to np array
                run_once = 1 # only grab ROIs once
            if filename.startswith(x):
                if names == "reading_hours": # {CHANGE THIS TO FOLDER FOR FIRST IND VAR}
                    df = pd.read_csv(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs", filename), usecols=["t value"])
                    df_var1[:,vars_dict[x]:(vars_dict[x] + 1)] = df.values # need to slice x:x+1 rather than x
                    df = pd.read_csv(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs", filename), usecols=["Pr(>|t|)"])
                    pvalues_var1[:,vars_dict[x]:(vars_dict[x] + 1)] = df.values 
                    print(filename)
                if names == "reading_years": # {CHANGE THIS TO FOLDER FOR SECOND IND VAR}
                    df = pd.read_csv(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs", filename), usecols=["t value"])
                    df_var2[:,vars_dict[x]:(vars_dict[x] + 1)] = df.values
                    df = pd.read_csv(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs", filename), usecols=["Pr(>|t|)"])
                    pvalues_var2[:,vars_dict[x]:(vars_dict[x] + 1)] = df.values 
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
sns.heatmap(df_var1, cmap=cmap, vmin=-4, vmax=4, cbar=False, ax=ax1, xticklabels=xlabels, yticklabels=False)
sns.heatmap(df_var2, cmap=cmap, vmin=-4, vmax=4, ax=ax2, xticklabels=xlabels, yticklabels=False, cbar=True, cbar_ax=ax3)

save_dir = os.path.join(os.getcwd(),'plots', 'reading', 'reading_heatmaps_AT.pdf') # {ENTER SAVE DIR HERE}
fig.savefig(save_dir, bbox_inches='tight')

#%Generate masked heat maps
fig, [ax1, ax2] = plt.subplots(1, 2, sharey=True, figsize=(15,10))
fig.tight_layout()
ax3 = fig.add_axes([1, 0.0275, 0.05, 0.95]) # add new axes for cbar [horizontal, vertical, width, height]
ax1.tick_params(labelsize=14, size=0, rotation=45)
ax2.tick_params(labelsize=14, size=0, rotation=45)
ax3.tick_params(labelsize=16, size=0, rotation=270)
sns.heatmap(df_var1, cmap=cmap, vmin=-4, vmax=4, cbar=False, ax=ax1, xticklabels=xlabels, yticklabels=False, mask=(pvalues_var1>0.05))
sns.heatmap(df_var2, cmap=cmap, vmin=-4, vmax=4, ax=ax2, xticklabels=xlabels, yticklabels=False, 
            mask=(pvalues_var2>0.05), cbar=True, cbar_ax=ax3)

save_dir = os.path.join(os.getcwd(),'plots', 'reading', 'reading_heatmaps_masked_AT.pdf') # {ENTER SAVE DIR HERE}
fig.savefig(save_dir, bbox_inches='tight')

# partition ordered parcellation names
i = 0 
for value in roi_parc:
    roi_parc[i] = substr_after(roi_parc[i], "fiber.at_")
    i = i + 1
# output .csv files with significant ROIs
save_dir = os.path.join(os.getcwd(),'plots', 'reading', 'lists')  # {ENTER SAVE DIR HERE}
output = sig_rois(df_var1, pvalues_var1, xlabels, roi_parc, save_dir, 'reading_hours_AT.csv')
output = sig_rois(df_var2, pvalues_var2, xlabels, roi_parc, save_dir, 'reading_years_AT.csv')

#%% Import t-scores and p-values (ASEG)
root_dir = os.path.join(os.getcwd(),'analyses')
# {ENTER ANY NUMBER OF DEP VARS BELOW WITH INCREMENTAL INTEGERS}
ind_vars = {"reading_hours" : "sports_activity_ss_read_hours_p", "reading_years" : "sports_activity_ss_read_years_p"} # {ENTER FOLDER NAME : VAR NAME}

vars_dict = {"smri_vol_subcort.aseg" : 0, "dmri_dti.fa_subcort.aseg" : 1, "dmri_dti.md_subcort.aseg" : 2, 
             "dmri_dti.ld_subcort.aseg" : 3, "dmri_dti.td_subcort.aseg" : 4, "dmri_rsi.n0_subcort.aseg" : 5, "dmri_rsi.nd_subcort.aseg": 6}
# {MAKE SURE # OF COLUMNS MATCH # OF DEP VARS BELOW}
df_var1 = np.zeros((30,7), dtype=float)
df_var2 = np.zeros((30,7), dtype=float)
pvalues_var1 = np.zeros((30,7), dtype=float)
pvalues_var2 = np.zeros((30,7), dtype=float)
run_once = 0
# iterate through folders and build dataframe
for names in ind_vars:
    for file in os.listdir(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs")): # iterate through a file list in dir
        filename = os.fsdecode(file)
        for x in vars_dict:
            if "aseg" in filename and run_once == 0:
                roi_parc = pd.read_csv(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs", filename), usecols=["dep_var"])
                roi_parc = roi_parc["dep_var"].values # convert to np array
                run_once = 1 # only grab ROIs once
            if filename.startswith(x):
                if names == "reading_hours": # {CHANGE THIS TO FOLDER FOR FIRST IND VAR}
                    df = pd.read_csv(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs", filename), usecols=["t value"])
                    df_var1[:,vars_dict[x]:(vars_dict[x] + 1)] = df.values # need to slice x:x+1 rather than x
                    df = pd.read_csv(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs", filename), usecols=["Pr(>|t|)"])
                    pvalues_var1[:,vars_dict[x]:(vars_dict[x] + 1)] = df.values 
                    print(filename)
                if names == "reading_years": # {CHANGE THIS TO FOLDER FOR SECOND IND VAR}
                    df = pd.read_csv(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs", filename), usecols=["t value"])
                    df_var2[:,vars_dict[x]:(vars_dict[x] + 1)] = df.values
                    df = pd.read_csv(os.path.join(root_dir, names, ind_vars[names], "tables", "table_coefs", filename), usecols=["Pr(>|t|)"])
                    pvalues_var2[:,vars_dict[x]:(vars_dict[x] + 1)] = df.values 
                    print(filename)

# Generate heat maps
cmap = sns.diverging_palette(220, 20, sep=20, as_cmap=True)

xlabels = ["Vol", "FA", "MD", "LD", "TD", "N0", "ND"]

fig, [ax1, ax2] = plt.subplots(1, 2, sharey=True, figsize=(20,15))
fig.tight_layout()
ax3 = fig.add_axes([1, 0.0275, 0.05, 0.95]) # add new axes for cbar [horizontal, vertical, width, height]
ax1.tick_params(labelsize=14, size=4, rotation=45)
ax2.tick_params(labelsize=14, size=4, rotation=45)
ax3.tick_params(labelsize=16, size=0, rotation=270)
spec = gridspec.GridSpec(ncols=2, nrows=1, figure=fig)
sns.heatmap(df_var1, cmap=cmap, vmin=-4, vmax=4, cbar=False, ax=ax1, xticklabels=xlabels, yticklabels=False)
sns.heatmap(df_var2, cmap=cmap, vmin=-4, vmax=4, ax=ax2, xticklabels=xlabels, yticklabels=False, cbar=True, cbar_ax=ax3)

save_dir = os.path.join(os.getcwd(),'plots', 'reading', 'reading_heatmaps_ASEG.pdf') # {ENTER SAVE DIR HERE}
fig.savefig(save_dir, bbox_inches='tight')

# Generate masked heat maps
fig, [ax1, ax2] = plt.subplots(1, 2, sharey=True, figsize=(15,10))
fig.tight_layout()
ax3 = fig.add_axes([1, 0.0275, 0.05, 0.95]) # add new axes for cbar [horizontal, vertical, width, height]
ax1.tick_params(labelsize=14, size=0, rotation=45)
ax2.tick_params(labelsize=14, size=0, rotation=45)
ax3.tick_params(labelsize=16, size=0, rotation=270)
sns.heatmap(df_var1, cmap=cmap, vmin=-4, vmax=4, cbar=False, ax=ax1, xticklabels=xlabels, yticklabels=False, mask=(pvalues_var1>0.05))
sns.heatmap(df_var2, cmap=cmap, vmin=-4, vmax=4, ax=ax2, xticklabels=xlabels, yticklabels=False, 
            mask=(pvalues_var2>0.05), cbar=True, cbar_ax=ax3)

save_dir = os.path.join(os.getcwd(),'plots', 'reading', 'reading_heatmaps_masked_ASEG.pdf') # {ENTER SAVE DIR HERE}
fig.savefig(save_dir, bbox_inches='tight')

# partition ordered parcellation names
i = 0 
for value in roi_parc:
    roi_parc[i] = substr_after(roi_parc[i], "aseg_")
    i = i + 1
# output .csv files with significant ROIs
save_dir = os.path.join(os.getcwd(),'plots', 'reading', 'lists')  # {ENTER SAVE DIR HERE}
output = sig_rois(df_var1, pvalues_var1, xlabels, roi_parc, save_dir, 'reading_hours_ASEG.csv')
output = sig_rois(df_var2, pvalues_var2, xlabels, roi_parc, save_dir, 'reading_years_ASEG.csv')