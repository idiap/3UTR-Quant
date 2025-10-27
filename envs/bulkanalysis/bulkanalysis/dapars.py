import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations
from statannotations.Annotator import Annotator
from scipy.stats import ttest_ind
import csv



def load_DaPars_outputs(path: str):
    # Load DaPars output
    dfs = []
    for dir in os.listdir(path):
        if 'chr' in dir:
            for file in os.listdir(os.path.join(path, dir)):
                # df = pd.read_csv(os.path.join(dir))
                if file.endswith('.txt'):
                    dfs.append(pd.read_csv(os.path.join(path, dir, file), sep='\t'))

    df = pd.concat(dfs)
    df.columns = [col.split('/')[-1] for col in df.columns]
    return df

def PDUI_boxplots_distributions(df: pd.DataFrame, treatments_dict: dict, fs: int = 17, **kwargs):
    ## Mean of the DPUI for the replicates
    for treatment, ids in treatments_dict.items():    
        df[f"{treatment}_PDUI"] = df[[f"{id}_PDUI" for id in ids]].mean(axis=1)

    fig, ax = plt.subplots(1, 1, figsize=(2, 2.5*len(treatments_dict)))
    # create stats annotations
    pairs = list(combinations([f"{treatment}_PDUI" for treatment in treatments_dict], 2))

    # compute p-values
    p_values = []
    for p1, p2 in pairs:
        p_values.append(ttest_ind(df[p1], 
                        df[p2], 
                        equal_var=False, 
                        alternative='two-sided',
                        nan_policy='omit')[1])
    formatted_pvalues = [f'P={pvalue:.2e}' for pvalue in p_values]

    plotting_parameters = {'data': df[[f"{treatment}_PDUI" for treatment in treatments_dict.keys()]],
                        'linewidth': 2.3, **kwargs}

    sns.boxplot(**plotting_parameters, ax=ax,boxprops={'edgecolor': 'k'},
                medianprops={'color':'black'},
                whiskerprops={'color':'k'},
                capprops={'color':'k'})
    annotator = Annotator(ax, pairs, **plotting_parameters)
    annotator.configure(fontsize=fs)

    annotator.set_custom_annotations(formatted_pvalues)
    annotator.annotate()
    plt.xticks(labels=treatments_dict.keys(), ticks=list(range(len(treatments_dict))), rotation=90);

