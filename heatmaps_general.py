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
import pandas as pd
import numpy as np
import os
import csv
from itertools import zip_longest

#%% Function def

def substr_after(s, delim): # partition finds the matchand splits the string into a tuple of [0/1/2] = before/match/after, s = string, delim = 
    return s.partition(delim)[2] # return everything afterwards

def sig_rois(tvalues, pvalues, columns, rois, filepath, filename): # inputs: columns = list of structural parameter ind_var, rois = list of parcellations
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
      
def heatmap(t1, p1, save_dir, xlabels):
    cmap = sns.diverging_palette(220, 20, sep=20, as_cmap=True)
    fig, [ax1, ax2] = plt.subplots(1, 2, sharey=True, figsize=(15,10))
    fig.tight_layout()
    ax3 = fig.add_axes([1, 0.0275, 0.05, 0.95]) # add new axes for cbar [horizontal, vertical, width, height]
    ax1.tick_params(labelsize=14, size=0, rotation=45)
    ax2.tick_params(labelsize=14, size=0, rotation=45)
    ax3.tick_params(labelsize=16, size=0, rotation=270)
    sns.heatmap(t1, cmap=cmap, vmin=-4, vmax=4, cbar=False, ax=ax1, xticklabels=xlabels, yticklabels=False)
    sns.heatmap(t1, cmap=cmap, vmin=-4, vmax=4, ax=ax2, xticklabels=xlabels, yticklabels=False, mask=(p1>0.05), cbar=True, cbar_ax=ax3)

    fig.savefig(save_dir, bbox_inches='tight')
    
#%% Import t-scores and p-values (DESIKAN)
root_dir = os.path.join(os.getcwd(), "analyses") # {ENTER DIR OF ANALYSIS FOLDERS}

ind_var = 'affected' # {VAR NAME}
category = 'affectedYes' # {NAME OF THE ACTUAL VARIABLE IN .CSV (pick the factor level you want to compare for a categorical, otherwise same as ind_var)}

vars_dict = {"dmri_dti.full.fa.gm_cort.desikan" : 0, "dmri_dti.full.fa.wm_cort.desikan" : 1, "dmri_dti.full.fa.gwc_cort.desikan" : 2, # {ENTER ANY NUMBER OF DEP VARS BELOW WITH INCREMENTAL INTEGERS}
             "dmri_dti.full.md.gm_cort.desikan" : 3, "dmri_dti.full.md.wm_cort.desikan" : 4, "dmri_dti.full.md.gwc_cort.desikan" : 5, 
             "dmri_dti.full.ld.gm_cort.desikan" : 6, "dmri_dti.full.ld.wm_cort.desikan" : 7, "dmri_dti.full.ld.gwc_cort.desikan" : 8, 
             "dmri_rsi.nd.gm_cort.desikan" : 9, "dmri_rsi.nd.wm_cort.desikan" : 10, "dmri_rsi.nd.gwc_cort.desikan" : 11,  
             "dmri_rsi.n0.gm_cort.desikan" : 12, "dmri_rsi.n0.wm_cort.desikan" : 13, "dmri_rsi.n0.gwc_cort.desikan" : 14,             
             "smri_thick_cort.desikan" : 15, "smri_area_cort.desikan" : 16, "smri_vol_cort.desikan" : 17}

xlabels = ["FA (GM)", "FA (WM)", "FA (GWC)", # {CHANGE THESE TO MATCH DEP VAR ind_var}
           "MD (GM)", "MD (WM)", "MD (GWC)", 
           "LD (GM)", "LD (WM)", "LD (GWC)", 
           "ND (GM)", "ND (WM)", "ND (GWC)", 
           "N0 (GM)", "N0 (WM)", "N0 (GWC)", 
           "Thick", "Area", "Volume"]

tvalues_one = np.zeros((71,18), dtype=float) # {MAKE SURE # OF COLUMNS MATCH # OF DEP VARS BELOW}
pvalues_one = np.zeros((71,18), dtype=float)
run_once = 0
roi_parc = []
# iterate through folders and build dataframe
for file in os.listdir(os.path.join(root_dir, ind_var, "tables", "table_coefs")): # iterate through a file list in dir
    filename = os.fsdecode(file)
    for x in vars_dict:
        if "desikan" in filename and run_once == 0:
            roi_parc = pd.read_csv(os.path.join(root_dir, ind_var, "tables", "table_coefs", filename), usecols=["dep_var"])
            roi_parc = roi_parc["dep_var"].values # convert to np array
            run_once = 1 # only grab ROIs once
        if filename.startswith(x):
            df = pd.read_csv(os.path.join(root_dir, ind_var, "tables", "table_coefs", filename), usecols=["parameter_comp", "t value"])
            df = df.loc[df['parameter_comp'] == category]
            temp = df.loc[:,"t value"] # need to slice x:x+1 rather than x
            tvalues_one[:,vars_dict[x]:(vars_dict[x] + 1)] = temp[:, np.newaxis] 
            df = pd.read_csv(os.path.join(root_dir, ind_var, "tables", "table_coefs", filename), usecols=["parameter_comp", "Pr(>|t|)"])
            df = df.loc[df['parameter_comp'] == category]
            temp = df.loc[:,"Pr(>|t|)"] # need to slice x:x+1 rather than x
            pvalues_one[:,vars_dict[x]:(vars_dict[x] + 1)] = temp[:, np.newaxis] 
            print(filename)
                
# Generate heat maps
save_dir = os.path.join(os.getcwd(),'plots', ind_var, str(ind_var) + '_heatmaps.pdf') # {ENTER SAVE DIR HERE}
heatmap(tvalues_one, pvalues_one, save_dir, xlabels)

# partition ordered parcellation ind_var
i = 0 
for value in roi_parc:
    roi_parc[i] = substr_after(roi_parc[i], "desikan_")
    i = i + 1
# output .csv files with significant ROIs
save_dir = os.path.join(os.getcwd(),'plots', ind_var, 'lists')  # {ENTER SAVE DIR HERE}
output = sig_rois(tvalues_one, pvalues_one, xlabels, roi_parc, save_dir, str(ind_var) + '.csv')

#%% Import t-scores and p-values (DTI ATLAS)

vars_dict = {"dmri_rsi.vol_fiber" : 0, "dmri_dti.fa_fiber.at" : 1, "dmri_dti.md_fiber.at" : 2, # {ENTER ANY NUMBER OF DEP VARS BELOW WITH INCREMENTAL INTEGERS}
             "dmri_dti.ld_fiber.at" : 3, "dmri_dti.td_fiber.at" : 4, "dmri_rsi.n0_fiber.at" : 5, "dmri_rsi.nd_fiber.at": 6}

xlabels = ["Vol", "FA", "MD", "LD", "TD", "N0", "ND"]
# {MAKE SURE # OF COLUMNS MATCH # OF DEP VARS BELOW}
tvalues_one = np.zeros((42,7), dtype=float)
pvalues_one = np.zeros((42,7), dtype=float)
run_once = 0
roi_parc = []
# iterate through folders and build dataframe
for file in os.listdir(os.path.join(root_dir, ind_var, "tables", "table_coefs")): # iterate through a file list in dir
    filename = os.fsdecode(file)
    for x in vars_dict:
        if "fiber.at" in filename and run_once == 0:
            roi_parc = pd.read_csv(os.path.join(root_dir, ind_var, "tables", "table_coefs", filename), usecols=["dep_var"])
            roi_parc = roi_parc["dep_var"].values # convert to np array
            run_once = 1 # only grab ROIs once
        if filename.startswith(x):
            df = pd.read_csv(os.path.join(root_dir, ind_var, "tables", "table_coefs", filename), usecols=["parameter_comp", "t value"])
            df = df.loc[df['parameter_comp'] == category]
            temp = df.loc[:,"t value"] # need to slice x:x+1 rather than x
            tvalues_one[:,vars_dict[x]:(vars_dict[x] + 1)] = temp[:, np.newaxis] 
            df = pd.read_csv(os.path.join(root_dir, ind_var, "tables", "table_coefs", filename), usecols=["parameter_comp", "Pr(>|t|)"])
            df = df.loc[df['parameter_comp'] == category]
            temp = df.loc[:,"Pr(>|t|)"] # need to slice x:x+1 rather than x
            pvalues_one[:,vars_dict[x]:(vars_dict[x] + 1)] = temp[:, np.newaxis] 
            print(filename)
        
# Generate heat maps
save_dir = os.path.join(os.getcwd(),'plots', ind_var, str(ind_var) + '_heatmaps_AT.pdf') # {ENTER SAVE DIR HERE}
heatmap(tvalues_one, pvalues_one, save_dir, xlabels)

# partition ordered parcellation ind_var
i = 0 
for value in roi_parc:
    roi_parc[i] = substr_after(roi_parc[i], "fiber.at_")
    i = i + 1
# output .csv files with significant ROIs
save_dir = os.path.join(os.getcwd(),'plots', ind_var, 'lists')  # {ENTER SAVE DIR HERE}
output = sig_rois(tvalues_one, pvalues_one, xlabels, roi_parc, save_dir, str(ind_var) + '_AT.csv')

#%% Import t-scores and p-values (ASEG)

vars_dict = {"smri_vol_subcort.aseg" : 0, "dmri_dti.fa_subcort.aseg" : 1, "dmri_dti.md_subcort.aseg" : 2, # {ENTER ANY NUMBER OF DEP VARS BELOW WITH INCREMENTAL INTEGERS}
             "dmri_dti.ld_subcort.aseg" : 3, "dmri_dti.td_subcort.aseg" : 4, "dmri_rsi.n0_subcort.aseg" : 5, "dmri_rsi.nd_subcort.aseg": 6}

xlabels = ["Vol", "FA", "MD", "LD", "TD", "N0", "ND"]
# {MAKE SURE # OF COLUMNS MATCH # OF DEP VARS BELOW}
tvalues_one = np.zeros((30,7), dtype=float)
pvalues_one = np.zeros((30,7), dtype=float)
run_once = 0
roi_parc = []
# iterate through folders and build dataframe
for file in os.listdir(os.path.join(root_dir, ind_var, "tables", "table_coefs")): # iterate through a file list in dir
    filename = os.fsdecode(file)
    for x in vars_dict:
        if "aseg" in filename and run_once == 0:
            roi_parc = pd.read_csv(os.path.join(root_dir, ind_var, "tables", "table_coefs", filename), usecols=["dep_var"])
            roi_parc = roi_parc["dep_var"].values # convert to np array
            run_once = 1 # only grab ROIs once
        if filename.startswith(x):
            df = pd.read_csv(os.path.join(root_dir, ind_var, "tables", "table_coefs", filename), usecols=["parameter_comp", "t value"])
            df = df.loc[df['parameter_comp'] == category]
            temp = df.loc[:,"t value"] # need to slice x:x+1 rather than x
            tvalues_one[:,vars_dict[x]:(vars_dict[x] + 1)] = temp[:, np.newaxis] 
            df = pd.read_csv(os.path.join(root_dir, ind_var, "tables", "table_coefs", filename), usecols=["parameter_comp", "Pr(>|t|)"])
            df = df.loc[df['parameter_comp'] == category]
            temp = df.loc[:,"Pr(>|t|)"] # need to slice x:x+1 rather than x
            pvalues_one[:,vars_dict[x]:(vars_dict[x] + 1)] = temp[:, np.newaxis] 
            print(filename)
                
# Generate heat maps
save_dir = os.path.join(os.getcwd(),'plots', ind_var, str(ind_var) + '_heatmaps_ASEG.pdf') # {ENTER SAVE DIR HERE}
heatmap(tvalues_one, pvalues_one, save_dir, xlabels)

# partition ordered parcellation ind_var
i = 0 
for value in roi_parc:
    roi_parc[i] = substr_after(roi_parc[i], "aseg_")
    i = i + 1
# output .csv files with significant ROIs
save_dir = os.path.join(os.getcwd(),'plots', ind_var, 'lists')  # {ENTER SAVE DIR HERE}
output = sig_rois(tvalues_one, pvalues_one, xlabels, roi_parc, save_dir, str(ind_var) + '_ASEG.csv')